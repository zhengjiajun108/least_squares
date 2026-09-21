#include "XlsxReader.h"

#include <cstdint>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<uint8_t> readFile(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f)
        throw std::runtime_error("cannot open file: " + path);
    f.seekg(0, std::ios::end);
    const std::streamoff size = f.tellg();
    f.seekg(0, std::ios::beg);
    std::vector<uint8_t> data(static_cast<std::size_t>(size));
    if (size > 0)
        f.read(reinterpret_cast<char*>(data.data()), size);
    return data;
}

uint16_t le16(const uint8_t* p) {
    return static_cast<uint16_t>(p[0] | (p[1] << 8));
}

uint32_t le32(const uint8_t* p) {
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

// ---------------------- DEFLATE (inflate) ----------------------
// 参考 zlib 的 puff.c 实现，支持 stored / fixed / dynamic 三种块

const int kMaxBits = 15;
const int kMaxLcodes = 286;
const int kMaxDcodes = 30;
const int kMaxCodes = kMaxLcodes + kMaxDcodes;
const int kFixLcodes = 288;

struct Huffman {
    short count[kMaxBits + 1];
    short symbol[kMaxCodes];
};

struct Inflater {
    const uint8_t* in;
    std::size_t inlen;
    std::size_t incnt;
    uint32_t bitbuf;
    int bitcnt;
    std::vector<uint8_t>* out;

    int bits(int need) {
        uint32_t val = bitbuf;
        while (bitcnt < need) {
            if (incnt == inlen)
                throw std::runtime_error("inflate: out of input");
            val |= static_cast<uint32_t>(in[incnt++]) << bitcnt;
            bitcnt += 8;
        }
        bitbuf = val >> need;
        bitcnt -= need;
        return static_cast<int>(val & ((1u << need) - 1));
    }
};

int decode(Inflater& s, const Huffman& h) {
    int len = 1, code = 0, first = 0, index = 0;
    for (; len <= kMaxBits; ++len) {
        code |= s.bits(1);
        const int count = h.count[len];
        if (code - count < first)
            return h.symbol[index + (code - first)];
        index += count;
        first += count;
        first <<= 1;
        code <<= 1;
    }
    return -10;
}

int construct(Huffman& h, const short* length, int n) {
    for (int len = 0; len <= kMaxBits; ++len)
        h.count[len] = 0;
    for (int symbol = 0; symbol < n; ++symbol)
        h.count[length[symbol]]++;
    if (h.count[0] == n)
        return 0;
    int left = 1;
    for (int len = 1; len <= kMaxBits; ++len) {
        left <<= 1;
        left -= h.count[len];
        if (left < 0)
            return left;
    }
    short offs[kMaxBits + 1];
    offs[1] = 0;
    for (int len = 1; len < kMaxBits; ++len)
        offs[len + 1] = static_cast<short>(offs[len] + h.count[len]);
    for (int symbol = 0; symbol < n; ++symbol)
        if (length[symbol] != 0)
            h.symbol[offs[length[symbol]]++] = static_cast<short>(symbol);
    return left;
}

const short kLens[29] = {3,  4,  5,  6,  7,  8,  9,  10, 11,  13,  15,  17,
                         19, 23, 27, 31, 35, 43, 51, 59, 67,  83,  99,  115,
                         131, 163, 195, 227, 258};
const short kLext[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
                         2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
const short kDists[30] = {1,    2,    3,    4,    5,    7,    9,    13,
                          17,   25,   33,   49,   65,   97,   129,  193,
                          257,  385,  513,  769,  1025, 1537, 2049, 3073,
                          4097, 6145, 8193, 12289, 16385, 24577};
const short kDext[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6,
                         6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

int codes(Inflater& s, const Huffman& lencode, const Huffman& distcode) {
    for (;;) {
        int symbol = decode(s, lencode);
        if (symbol < 0)
            return symbol;
        if (symbol < 256) {
            s.out->push_back(static_cast<uint8_t>(symbol));
        } else if (symbol == 256) {
            return 0;
        } else {
            symbol -= 257;
            if (symbol >= 29)
                return -9;
            const int len = kLens[symbol] + s.bits(kLext[symbol]);
            const int dsym = decode(s, distcode);
            if (dsym < 0)
                return dsym;
            const int dist = kDists[dsym] + s.bits(kDext[dsym]);
            if (static_cast<std::size_t>(dist) > s.out->size())
                return -11;
            for (int i = 0; i < len; ++i)
                s.out->push_back((*s.out)[s.out->size() - dist]);
        }
    }
}

int stored(Inflater& s) {
    s.bitbuf = 0;
    s.bitcnt = 0;
    if (s.incnt + 4 > s.inlen)
        return 2;
    const int len = s.in[s.incnt] | (s.in[s.incnt + 1] << 8);
    const int nlen = s.in[s.incnt + 2] | (s.in[s.incnt + 3] << 8);
    s.incnt += 4;
    if ((len ^ 0xffff) != nlen)
        return -4;
    if (s.incnt + len > s.inlen)
        return 2;
    for (int i = 0; i < len; ++i)
        s.out->push_back(s.in[s.incnt++]);
    return 0;
}

int fixed(Inflater& s) {
    Huffman lencode, distcode;
    short lengths[kFixLcodes];
    int symbol;
    for (symbol = 0; symbol < 144; ++symbol)
        lengths[symbol] = 8;
    for (; symbol < 256; ++symbol)
        lengths[symbol] = 9;
    for (; symbol < 280; ++symbol)
        lengths[symbol] = 7;
    for (; symbol < kFixLcodes; ++symbol)
        lengths[symbol] = 8;
    int err = construct(lencode, lengths, kFixLcodes);
    if (err)
        return err;
    for (symbol = 0; symbol < kMaxDcodes; ++symbol)
        lengths[symbol] = 5;
    err = construct(distcode, lengths, kMaxDcodes);
    if (err)
        return err;
    return codes(s, lencode, distcode);
}

int dynamic(Inflater& s) {
    static const short kOrder[19] = {16, 17, 18, 0, 8,  7, 9,  6, 10, 5,
                                     11, 4,  12, 3, 13, 2, 14, 1, 15};
    Huffman lencode, distcode;
    short lengths[kMaxCodes];
    for (int i = 0; i < kMaxCodes; ++i)
        lengths[i] = 0;

    const int nlen = s.bits(5) + 257;
    const int ndist = s.bits(5) + 1;
    const int ncode = s.bits(4) + 4;
    if (nlen > kMaxLcodes || ndist > kMaxDcodes)
        return -3;

    short lens[19];
    for (int i = 0; i < 19; ++i)
        lens[i] = 0;
    int index = 0;
    for (; index < ncode; ++index)
        lens[kOrder[index]] = static_cast<short>(s.bits(3));

    int err = construct(lencode, lens, 19);
    if (err)
        return -4;

    index = 0;
    while (index < nlen + ndist) {
        const int symbol = decode(s, lencode);
        if (symbol < 0)
            return symbol;
        if (symbol < 16) {
            lengths[index++] = static_cast<short>(symbol);
        } else {
            int len = 0;
            int count;
            if (symbol == 16) {
                if (index == 0)
                    return -5;
                len = lengths[index - 1];
                count = 3 + s.bits(2);
            } else if (symbol == 17) {
                count = 3 + s.bits(3);
            } else {
                count = 11 + s.bits(7);
            }
            if (index + count > nlen + ndist)
                return -6;
            while (count--)
                lengths[index++] = static_cast<short>(len);
        }
    }
    if (lengths[256] == 0)
        return -7;

    err = construct(lencode, lengths, nlen);
    if (err && (err < 0 || nlen != lencode.count[0] + lencode.count[1]))
        return -8;
    err = construct(distcode, lengths + nlen, ndist);
    if (err && (err < 0 || ndist != distcode.count[0] + distcode.count[1]))
        return -9;
    return codes(s, lencode, distcode);
}

std::vector<uint8_t> inflate(const std::vector<uint8_t>& in) {
    Inflater s{in.data(), in.size(), 0, 0, 0, nullptr};
    std::vector<uint8_t> out;
    s.out = &out;
    int last;
    do {
        last = s.bits(1);
        const int type = s.bits(2);
        const int err = (type == 0)   ? stored(s)
                        : (type == 1) ? fixed(s)
                        : (type == 2) ? dynamic(s)
                                      : -1;
        if (err != 0)
            throw std::runtime_error("inflate failed");
    } while (!last);
    return out;
}

// ---------------------- ZIP 解析 ----------------------

struct ZipEntry {
    std::string name;
    uint16_t method;
    uint32_t compSize;
    uint32_t localOffset;
};

std::vector<ZipEntry> readCentralDirectory(const std::vector<uint8_t>& d) {
    if (d.size() < 22)
        throw std::runtime_error("file too small to be xlsx");
    std::size_t eocd = std::string::npos;
    const std::size_t lo = d.size() >= 65557 ? d.size() - 65557 : 0;
    for (std::size_t i = d.size() - 22 + 1; i-- > lo;) {
        if (le32(&d[i]) == 0x06054b50u) {
            eocd = i;
            break;
        }
    }
    if (eocd == std::string::npos)
        throw std::runtime_error("end of central directory not found");

    const uint32_t count = le16(&d[eocd + 10]);
    const uint32_t cdOffset = le32(&d[eocd + 16]);

    std::vector<ZipEntry> entries;
    std::size_t p = cdOffset;
    for (uint32_t k = 0; k < count; ++k) {
        if (p + 46 > d.size() || le32(&d[p]) != 0x02014b50u)
            break;
        ZipEntry e;
        e.method = le16(&d[p + 10]);
        e.compSize = le32(&d[p + 20]);
        const uint16_t nameLen = le16(&d[p + 28]);
        const uint16_t extraLen = le16(&d[p + 30]);
        const uint16_t commentLen = le16(&d[p + 32]);
        e.localOffset = le32(&d[p + 42]);
        e.name.assign(reinterpret_cast<const char*>(&d[p + 46]), nameLen);
        entries.push_back(e);
        p += 46 + nameLen + extraLen + commentLen;
    }
    return entries;
}

std::vector<uint8_t> extractEntry(const std::vector<uint8_t>& d,
                                  const ZipEntry& e) {
    const std::size_t p = e.localOffset;
    if (p + 30 > d.size() || le32(&d[p]) != 0x04034b50u)
        throw std::runtime_error("bad local file header");
    const uint16_t nameLen = le16(&d[p + 26]);
    const uint16_t extraLen = le16(&d[p + 28]);
    const std::size_t dataOff = p + 30 + nameLen + extraLen;
    if (dataOff + e.compSize > d.size())
        throw std::runtime_error("zip data out of range");
    std::vector<uint8_t> raw(d.begin() + dataOff,
                             d.begin() + dataOff + e.compSize);
    if (e.method == 0)
        return raw;
    if (e.method == 8)
        return inflate(raw);
    throw std::runtime_error("unsupported zip compression method");
}

// ---------------------- XML 解析 ----------------------

std::vector<std::pair<double, double>> parseSheet(const std::string& xml) {
    std::map<int, double> xs, ys;
    std::size_t pos = 0;
    for (;;) {
        const std::size_t c = xml.find("<c ", pos);
        if (c == std::string::npos)
            break;
        const std::size_t end = xml.find("</c>", c);
        if (end == std::string::npos)
            break;
        const std::string cell = xml.substr(c, end - c);
        pos = end + 4;

        const std::size_t rp = cell.find("r=\"");
        if (rp == std::string::npos)
            continue;
        const std::size_t rq = cell.find('"', rp + 3);
        const std::string ref = cell.substr(rp + 3, rq - (rp + 3));

        std::string col;
        int row = 0;
        for (const char ch : ref) {
            if (ch >= 'A' && ch <= 'Z')
                col += ch;
            else if (ch >= '0' && ch <= '9')
                row = row * 10 + (ch - '0');
        }

        if (cell.find("t=\"s\"") != std::string::npos)
            continue;

        const std::size_t vp = cell.find("<v>");
        if (vp == std::string::npos)
            continue;
        const std::size_t vq = cell.find("</v>", vp);
        if (vq == std::string::npos)
            continue;
        const double value = std::stod(cell.substr(vp + 3, vq - (vp + 3)));

        if (col == "A")
            xs[row] = value;
        else if (col == "B")
            ys[row] = value;
    }

    std::vector<std::pair<double, double>> points;
    for (const auto& kv : xs) {
        const auto it = ys.find(kv.first);
        if (it != ys.end())
            points.emplace_back(kv.second, it->second);
    }
    return points;
}

}  // namespace

void readPointsFromXlsx(const std::string& path,
                        std::vector<double>& x,
                        std::vector<double>& y) {
    const std::vector<uint8_t> data = readFile(path);
    const std::vector<ZipEntry> entries = readCentralDirectory(data);

    const ZipEntry* sheet = nullptr;
    for (const auto& e : entries) {
        if (e.name == "xl/worksheets/sheet1.xml") {
            sheet = &e;
            break;
        }
    }
    if (sheet == nullptr)
        throw std::runtime_error("xl/worksheets/sheet1.xml not found");

    const std::vector<uint8_t> xmlBytes = extractEntry(data, *sheet);
    const std::string xml(xmlBytes.begin(), xmlBytes.end());
    const std::vector<std::pair<double, double>> points = parseSheet(xml);

    x.clear();
    y.clear();
    x.reserve(points.size());
    y.reserve(points.size());
    for (const auto& p : points) {
        x.push_back(p.first);
        y.push_back(p.second);
    }
}
