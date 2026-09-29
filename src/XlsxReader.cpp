// xlsx表格读写操作
#include "XlsxReader.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
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

const ZipEntry* findEntry(const std::vector<ZipEntry>& entries,
                          const std::string& name) {
    for (const auto& e : entries)
        if (e.name == name)
            return &e;
    return nullptr;
}

std::string toString(const std::vector<uint8_t>& bytes) {
    return std::string(bytes.begin(), bytes.end());
}

// ---------------------- XML 工具 ----------------------

std::string getAttr(const std::string& element, const std::string& name) {
    const std::string key = name + "=\"";
    const std::size_t p = element.find(key);
    if (p == std::string::npos)
        return "";
    const std::size_t begin = p + key.size();
    const std::size_t q = element.find('"', begin);
    if (q == std::string::npos)
        return "";
    return element.substr(begin, q - begin);
}

std::string unescapeXml(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '&') {
            if (s.compare(i, 5, "&amp;") == 0) {
                out += '&';
                i += 4;
            } else if (s.compare(i, 4, "&lt;") == 0) {
                out += '<';
                i += 3;
            } else if (s.compare(i, 4, "&gt;") == 0) {
                out += '>';
                i += 3;
            } else if (s.compare(i, 6, "&quot;") == 0) {
                out += '"';
                i += 5;
            } else if (s.compare(i, 6, "&apos;") == 0) {
                out += '\'';
                i += 5;
            } else {
                out += s[i];
            }
        } else {
            out += s[i];
        }
    }
    return out;
}

// 拼接元素内所有 <t>...</t> 的文本(用于共享字符串 / 内联字符串)
std::string extractText(const std::string& s) {
    std::string out;
    std::size_t p = 0;
    for (;;) {
        const std::size_t t = s.find("<t", p);
        if (t == std::string::npos)
            break;
        const std::size_t gt = s.find('>', t);
        if (gt == std::string::npos)
            break;
        const std::size_t close = s.find("</t>", gt);
        if (close == std::string::npos)
            break;
        out += s.substr(gt + 1, close - (gt + 1));
        p = close + 4;
    }
    return unescapeXml(out);
}

// 取 <v>...</v> 的原始内容
std::string extractRawV(const std::string& cell) {
    const std::size_t vp = cell.find("<v>");
    if (vp == std::string::npos)
        return "";
    const std::size_t vq = cell.find("</v>", vp);
    if (vq == std::string::npos)
        return "";
    return cell.substr(vp + 3, vq - (vp + 3));
}

// ---------------------- 工作表解析 ----------------------

struct Cell {
    bool hasValue = false;
    bool isString = false;
    double number = 0.0;
    std::string text;
};

// 解析共享字符串表 sharedStrings.xml，返回按索引排列的字符串
std::vector<std::string> parseSharedStrings(const std::string& xml) {
    std::vector<std::string> shared;
    std::size_t p = 0;
    for (;;) {
        const std::size_t s = xml.find("<si>", p);
        if (s == std::string::npos)
            break;
        const std::size_t e = xml.find("</si>", s);
        if (e == std::string::npos)
            break;
        shared.push_back(extractText(xml.substr(s, e - s)));
        p = e + 5;
    }
    return shared;
}

// 解析工作表: 返回 cells[行][列] = Cell，行列均为 1 起始
std::map<int, std::map<int, Cell>> parseSheet(
    const std::string& xml, const std::vector<std::string>& shared) {
    std::map<int, std::map<int, Cell>> cells;
    std::size_t pos = 0;
    for (;;) {
        const std::size_t c = xml.find("<c ", pos);
        if (c == std::string::npos)
            break;
        const std::size_t endClose = xml.find("</c>", c);
        const std::size_t endSelf = xml.find("/>", c);
        std::size_t end;
        std::size_t next;
        if (endSelf != std::string::npos &&
            (endClose == std::string::npos || endSelf < endClose)) {
            end = endSelf;
            next = endSelf + 2;
        } else {
            if (endClose == std::string::npos)
                break;
            end = endClose;
            next = endClose + 4;
        }
        const std::string cell = xml.substr(c, end - c);
        pos = next;

        const std::string ref = getAttr(cell, "r");
        if (ref.empty())
            continue;

        int col = 0;
        int row = 0;
        for (const char ch : ref) {
            if (ch >= 'A' && ch <= 'Z')
                col = col * 26 + (ch - 'A' + 1);
            else if (ch >= '0' && ch <= '9')
                row = row * 10 + (ch - '0');
        }

        Cell value;
        if (cell.find("t=\"inlineStr\"") != std::string::npos) {
            value.hasValue = true;
            value.isString = true;
            value.text = extractText(cell);
        } else if (cell.find("t=\"s\"") != std::string::npos) {
            const std::string v = extractRawV(cell);
            if (!v.empty()) {
                const std::size_t idx = static_cast<std::size_t>(std::stoul(v));
                value.hasValue = true;
                value.isString = true;
                if (idx < shared.size())
                    value.text = shared[idx];
            }
        } else if (cell.find("t=\"str\"") != std::string::npos) {
            value.hasValue = true;
            value.isString = true;
            value.text = unescapeXml(extractRawV(cell));
        } else {
            const std::string v = extractRawV(cell);
            if (!v.empty()) {
                value.hasValue = true;
                value.isString = false;
                value.number = std::stod(v);
            }
        }

        if (value.hasValue)
            cells[row][col] = value;
    }
    return cells;
}

// 由 workbook.xml 与 workbook.xml.rels 求工作表对应的 zip 路径
std::string resolveSheetPath(const std::vector<uint8_t>& data,
                             const std::vector<ZipEntry>& entries,
                             const std::string& sheetName) {
    const ZipEntry* wb = findEntry(entries, "xl/workbook.xml");
    if (wb == nullptr)
        throw std::runtime_error("xl/workbook.xml not found");
    const std::string wbxml = toString(extractEntry(data, *wb));

    std::string rid;
    std::size_t p = 0;
    for (;;) {
        const std::size_t s = wbxml.find("<sheet ", p);
        if (s == std::string::npos)
            break;
        const std::size_t e = wbxml.find('>', s);
        if (e == std::string::npos)
            break;
        const std::string el = wbxml.substr(s, e - s);
        if (getAttr(el, "name") == sheetName) {
            rid = getAttr(el, "r:id");
            break;
        }
        p = e + 1;
    }
    if (rid.empty())
        throw std::runtime_error("sheet not found: " + sheetName);

    const ZipEntry* rel = findEntry(entries, "xl/_rels/workbook.xml.rels");
    if (rel == nullptr)
        throw std::runtime_error("xl/_rels/workbook.xml.rels not found");
    const std::string relxml = toString(extractEntry(data, *rel));

    std::string target;
    p = 0;
    for (;;) {
        const std::size_t s = relxml.find("<Relationship ", p);
        if (s == std::string::npos)
            break;
        const std::size_t e = relxml.find('>', s);
        if (e == std::string::npos)
            break;
        const std::string el = relxml.substr(s, e - s);
        if (getAttr(el, "Id") == rid) {
            target = getAttr(el, "Target");
            break;
        }
        p = e + 1;
    }
    if (target.empty())
        throw std::runtime_error("relationship not found for sheet: " +
                                 sheetName);

    if (target[0] == '/')
        return target.substr(1);
    return "xl/" + target;
}

// ---------------------- 写出 (ZIP / worksheet) ----------------------

using PartList = std::vector<std::pair<std::string, std::vector<uint8_t>>>;

uint32_t crc32(const uint8_t* data, std::size_t len) {
    uint32_t crc = 0xFFFFFFFFu;
    for (std::size_t i = 0; i < len; ++i) {
        crc ^= data[i];
        for (int j = 0; j < 8; ++j)
            crc = (crc & 1u) ? ((crc >> 1) ^ 0xEDB88320u) : (crc >> 1);
    }
    return ~crc;
}

std::string colName(int col) {
    std::string s;
    while (col > 0) {
        const int rem = (col - 1) % 26;
        s.insert(s.begin(), static_cast<char>('A' + rem));
        col = (col - 1) / 26;
    }
    return s;
}

std::string escapeXml(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (const char ch : s) {
        switch (ch) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out += ch;
        }
    }
    return out;
}

std::string formatNumber(double v) {
    std::ostringstream os;
    os << std::setprecision(15) << v;
    return os.str();
}

std::string buildSheetXml(const std::vector<std::string>& headers,
                          const std::vector<std::vector<double>>& columns) {
    std::string xml =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/"
        "2006/main\"><sheetData>";

    xml += "<row r=\"1\">";
    for (std::size_t c = 0; c < headers.size(); ++c) {
        const std::string ref = colName(static_cast<int>(c) + 1) + "1";
        xml += "<c r=\"" + ref + "\" t=\"inlineStr\"><is><t>" +
               escapeXml(headers[c]) + "</t></is></c>";
    }
    xml += "</row>";

    std::size_t rows = 0;
    for (const auto& col : columns)
        rows = std::max(rows, col.size());
    for (std::size_t r = 0; r < rows; ++r) {
        xml += "<row r=\"" + std::to_string(r + 2) + "\">";
        for (std::size_t c = 0; c < columns.size(); ++c) {
            if (r >= columns[c].size())
                continue;
            const std::string ref =
                colName(static_cast<int>(c) + 1) + std::to_string(r + 2);
            xml += "<c r=\"" + ref + "\"><v>" + formatNumber(columns[c][r]) +
                   "</v></c>";
        }
        xml += "</row>";
    }
    xml += "</sheetData></worksheet>";
    return xml;
}

void appendLe16(std::vector<uint8_t>& b, uint16_t v) {
    b.push_back(static_cast<uint8_t>(v & 0xff));
    b.push_back(static_cast<uint8_t>((v >> 8) & 0xff));
}

void appendLe32(std::vector<uint8_t>& b, uint32_t v) {
    b.push_back(static_cast<uint8_t>(v & 0xff));
    b.push_back(static_cast<uint8_t>((v >> 8) & 0xff));
    b.push_back(static_cast<uint8_t>((v >> 16) & 0xff));
    b.push_back(static_cast<uint8_t>((v >> 24) & 0xff));
}

std::vector<uint8_t> toBytes(const std::string& s) {
    return std::vector<uint8_t>(s.begin(), s.end());
}

// 以不压缩(stored)方式写出 zip
void writeZip(const std::string& path, const PartList& parts) {
    struct Central {
        std::string name;
        uint32_t crc;
        uint32_t size;
        uint32_t offset;
    };
    std::vector<uint8_t> out;
    std::vector<Central> central;

    for (const auto& part : parts) {
        Central c;
        c.name = part.first;
        c.crc = crc32(part.second.data(), part.second.size());
        c.size = static_cast<uint32_t>(part.second.size());
        c.offset = static_cast<uint32_t>(out.size());

        appendLe32(out, 0x04034b50u);
        appendLe16(out, 20);
        appendLe16(out, 0);
        appendLe16(out, 0);     // stored
        appendLe16(out, 0);     // time
        appendLe16(out, 0x21);  // date (1980-01-01)
        appendLe32(out, c.crc);
        appendLe32(out, c.size);
        appendLe32(out, c.size);
        appendLe16(out, static_cast<uint16_t>(c.name.size()));
        appendLe16(out, 0);
        out.insert(out.end(), c.name.begin(), c.name.end());
        out.insert(out.end(), part.second.begin(), part.second.end());
        central.push_back(c);
    }

    const uint32_t cdOffset = static_cast<uint32_t>(out.size());
    for (const auto& c : central) {
        appendLe32(out, 0x02014b50u);
        appendLe16(out, 20);
        appendLe16(out, 20);
        appendLe16(out, 0);
        appendLe16(out, 0);
        appendLe16(out, 0);
        appendLe16(out, 0x21);
        appendLe32(out, c.crc);
        appendLe32(out, c.size);
        appendLe32(out, c.size);
        appendLe16(out, static_cast<uint16_t>(c.name.size()));
        appendLe16(out, 0);
        appendLe16(out, 0);
        appendLe16(out, 0);
        appendLe16(out, 0);
        appendLe32(out, 0);
        appendLe32(out, c.offset);
        out.insert(out.end(), c.name.begin(), c.name.end());
    }
    const uint32_t cdSize = static_cast<uint32_t>(out.size()) - cdOffset;

    appendLe32(out, 0x06054b50u);
    appendLe16(out, 0);
    appendLe16(out, 0);
    appendLe16(out, static_cast<uint16_t>(central.size()));
    appendLe16(out, static_cast<uint16_t>(central.size()));
    appendLe32(out, cdSize);
    appendLe32(out, cdOffset);
    appendLe16(out, 0);

    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    if (!f)
        throw std::runtime_error("cannot write file: " + path);
    f.write(reinterpret_cast<const char*>(out.data()),
            static_cast<std::streamsize>(out.size()));
}

PartList::iterator findPart(PartList& parts, const std::string& name) {
    for (auto it = parts.begin(); it != parts.end(); ++it)
        if (it->first == name)
            return it;
    return parts.end();
}

PartList buildEmptyPackage() {
    PartList parts;
    parts.emplace_back(
        "[Content_Types].xml",
        toBytes("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
                "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/"
                "content-types\"><Default Extension=\"rels\" "
                "ContentType=\"application/vnd.openxmlformats-package."
                "relationships+xml\"/><Default Extension=\"xml\" "
                "ContentType=\"application/xml\"/><Override "
                "PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd."
                "openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
                "</Types>"));
    parts.emplace_back(
        "_rels/.rels",
        toBytes("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
                "<Relationships xmlns=\"http://schemas.openxmlformats.org/"
                "package/2006/relationships\"><Relationship Id=\"rId1\" "
                "Type=\"http://schemas.openxmlformats.org/officeDocument/2006/"
                "relationships/officeDocument\" Target=\"xl/workbook.xml\"/>"
                "</Relationships>"));
    parts.emplace_back(
        "xl/workbook.xml",
        toBytes(
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
            "<workbook xmlns=\"http://schemas.openxmlformats.org/"
            "spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats."
            "org/officeDocument/2006/relationships\"><sheets></sheets>"
            "</workbook>"));
    parts.emplace_back(
        "xl/_rels/workbook.xml.rels",
        toBytes("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
                "<Relationships xmlns=\"http://schemas.openxmlformats.org/"
                "package/2006/relationships\"></Relationships>"));
    return parts;
}

// 在 parts 中新增或覆盖指定工作表，并同步 workbook / rels / content-types
void upsertSheet(PartList& parts, const std::string& sheetName,
                 const std::string& sheetXml) {
    auto wbIt = findPart(parts, "xl/workbook.xml");
    auto relIt = findPart(parts, "xl/_rels/workbook.xml.rels");
    auto ctIt = findPart(parts, "[Content_Types].xml");
    if (wbIt == parts.end() || relIt == parts.end())
        throw std::runtime_error("workbook parts missing");

    const std::string wb = toString(wbIt->second);
    const std::string rels = toString(relIt->second);

    std::string rid;
    std::size_t p = 0;
    for (;;) {
        const std::size_t s = wb.find("<sheet ", p);
        if (s == std::string::npos)
            break;
        const std::size_t e = wb.find('>', s);
        if (e == std::string::npos)
            break;
        const std::string el = wb.substr(s, e - s);
        if (getAttr(el, "name") == sheetName) {
            rid = getAttr(el, "r:id");
            break;
        }
        p = e + 1;
    }

    if (!rid.empty()) {
        std::string target;
        p = 0;
        for (;;) {
            const std::size_t s = rels.find("<Relationship ", p);
            if (s == std::string::npos)
                break;
            const std::size_t e = rels.find('>', s);
            if (e == std::string::npos)
                break;
            const std::string el = rels.substr(s, e - s);
            if (getAttr(el, "Id") == rid) {
                target = getAttr(el, "Target");
                break;
            }
            p = e + 1;
        }
        if (target.empty())
            throw std::runtime_error("relationship missing for sheet");
        const std::string path =
            (target[0] == '/') ? target.substr(1) : "xl/" + target;
        auto it = findPart(parts, path);
        if (it == parts.end())
            parts.emplace_back(path, toBytes(sheetXml));
        else
            it->second = toBytes(sheetXml);
        return;
    }

    int maxSheet = 0;
    const std::string pre = "xl/worksheets/sheet";
    const std::string suf = ".xml";
    for (const auto& part : parts) {
        if (part.first.rfind(pre, 0) == 0 &&
            part.first.size() > pre.size() + suf.size() &&
            part.first.compare(part.first.size() - suf.size(), suf.size(),
                               suf) == 0) {
            const std::string num =
                part.first.substr(pre.size(),
                                  part.first.size() - pre.size() - suf.size());
            bool digits = !num.empty();
            for (const char ch : num)
                if (!(ch >= '0' && ch <= '9'))
                    digits = false;
            if (digits)
                maxSheet = std::max(maxSheet, std::stoi(num));
        }
    }
    const int newSheet = maxSheet + 1;
    const std::string newPath =
        "xl/worksheets/sheet" + std::to_string(newSheet) + ".xml";

    int maxRid = 0;
    for (const std::string& src : {wb, rels}) {
        std::size_t q = 0;
        for (;;) {
            const std::size_t r = src.find("rId", q);
            if (r == std::string::npos)
                break;
            std::size_t j = r + 3;
            int val = 0;
            bool any = false;
            while (j < src.size() &&
                   std::isdigit(static_cast<unsigned char>(src[j]))) {
                val = val * 10 + (src[j] - '0');
                ++j;
                any = true;
            }
            if (any) {
                maxRid = std::max(maxRid, val);
                q = j;
            } else {
                q = r + 3;
            }
        }
    }
    const std::string newRid = "rId" + std::to_string(maxRid + 1);

    int maxSheetId = 0;
    {
        std::size_t q = 0;
        for (;;) {
            std::size_t s = wb.find("sheetId=\"", q);
            if (s == std::string::npos)
                break;
            s += 9;
            const std::size_t e = wb.find('"', s);
            maxSheetId = std::max(maxSheetId, std::stoi(wb.substr(s, e - s)));
            q = e + 1;
        }
    }
    const int newSheetId = maxSheetId + 1;

    std::string newWb = wb;
    const std::size_t sp = newWb.find("</sheets>");
    const std::string sheetEl = "<sheet name=\"" + escapeXml(sheetName) +
                                "\" sheetId=\"" + std::to_string(newSheetId) +
                                "\" r:id=\"" + newRid + "\"/>";
    newWb.insert(sp, sheetEl);
    wbIt->second = toBytes(newWb);

    std::string newRels = rels;
    const std::size_t rp = newRels.find("</Relationships>");
    const std::string relEl =
        "<Relationship Id=\"" + newRid +
        "\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/"
        "relationships/worksheet\" Target=\"worksheets/sheet" +
        std::to_string(newSheet) + ".xml\"/>";
    newRels.insert(rp, relEl);
    relIt->second = toBytes(newRels);

    if (ctIt != parts.end()) {
        std::string ct = toString(ctIt->second);
        const std::size_t cp = ct.find("</Types>");
        const std::string ov = "<Override PartName=\"/" + newPath +
                               "\" ContentType=\"application/vnd."
                               "openxmlformats-officedocument.spreadsheetml."
                               "worksheet+xml\"/>";
        ct.insert(cp, ov);
        ctIt->second = toBytes(ct);
    }

    parts.emplace_back(newPath, toBytes(sheetXml));
}

}  // namespace

void readPointsFromXlsx(const std::string& path,
                        const std::string& sheetName,
                        const std::string& xHeader,
                        const std::string& yHeader,
                        std::vector<double>& x,
                        std::vector<double>& y) {
    const std::vector<uint8_t> data = readFile(path);
    const std::vector<ZipEntry> entries = readCentralDirectory(data);

    std::vector<std::string> shared;
    if (const ZipEntry* ss = findEntry(entries, "xl/sharedStrings.xml"))
        shared = parseSharedStrings(toString(extractEntry(data, *ss)));

    const std::string sheetPath = resolveSheetPath(data, entries, sheetName);
    const ZipEntry* sheet = findEntry(entries, sheetPath);
    if (sheet == nullptr)
        throw std::runtime_error("worksheet not found: " + sheetPath);

    const std::string xml = toString(extractEntry(data, *sheet));
    const std::map<int, std::map<int, Cell>> cells = parseSheet(xml, shared);

    int headerRow = -1;
    int xCol = -1;
    int yCol = -1;
    for (const auto& row : cells) {
        for (const auto& col : row.second) {
            const Cell& cell = col.second;
            if (!cell.isString)
                continue;
            if (cell.text == xHeader) {
                headerRow = row.first;
                xCol = col.first;
            } else if (cell.text == yHeader) {
                headerRow = row.first;
                yCol = col.first;
            }
        }
    }
    if (headerRow < 0 || xCol < 0 || yCol < 0)
        throw std::runtime_error("header not found: " + xHeader + " / " +
                                 yHeader);

    x.clear();
    y.clear();
    for (const auto& row : cells) {
        if (row.first <= headerRow)
            continue;
        const auto itx = row.second.find(xCol);
        const auto ity = row.second.find(yCol);
        if (itx == row.second.end() || ity == row.second.end())
            continue;
        if (itx->second.isString || ity->second.isString)
            continue;
        if (!itx->second.hasValue || !ity->second.hasValue)
            continue;
        x.push_back(itx->second.number);
        y.push_back(ity->second.number);
    }
}

void writeColumnsToXlsx(const std::string& path,
                        const std::string& sheetName,
                        const std::vector<std::string>& headers,
                        const std::vector<std::vector<double>>& columns) {
    if (headers.size() != columns.size())
        throw std::runtime_error("headers and columns size mismatch");

    const std::string sheetXml = buildSheetXml(headers, columns);

    PartList parts;
    bool loaded = false;
    try {
        const std::vector<uint8_t> data = readFile(path);
        const std::vector<ZipEntry> entries = readCentralDirectory(data);
        for (const auto& e : entries)
            parts.emplace_back(e.name, extractEntry(data, e));
        loaded = true;
    } catch (...) {
        loaded = false;
    }
    if (!loaded)
        parts = buildEmptyPackage();

    upsertSheet(parts, sheetName, sheetXml);
    writeZip(path, parts);
}
