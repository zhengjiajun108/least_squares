#include "core/PlaneFit.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace {

// 对行主序 N×3 矩阵 Q 做 one-sided Jacobi SVD：
// 逐对正交化 Q 的三列，收敛后
//   - q 每列的范数即奇异值(降序前)
//   - v 每列为右奇异向量，满足 Q = U * diag(sv) * v^T
void jacobiSvd3(std::vector<double>& q, std::size_t n, double v[3][3],
                double sv[3]) {
    // 初始化右奇异向量v为单位矩阵
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            v[i][j] = (i == j) ? 1.0 : 0.0;

    // 迭代收敛(60轮)
    const int kMaxSweeps = 60;
    for (int sweep = 0; sweep < kMaxSweeps; ++sweep) {
        double off = 0.0;
        // 对每一对列做 Jacobi 旋转
        // (0, 1), (0, 2), (1, 2)
        for (int p = 0; p < 2; ++p) {
            for (int r = p + 1; r < 3; ++r) {
                // 计算列 p 与列 r 的内积
                double a = 0.0, b = 0.0, c = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    const double qp = q[i * 3 + p];
                    const double qr = q[i * 3 + r];
                    // a = p列的平方和，b = r列的平方和，c = p列与r列的内积
                    a += qp * qp;
                    b += qr * qr;
                    c += qp * qr;
                }
                // 累加非对角量||c||，用于判断是否收敛
                off += std::fabs(c);
                // if (std::fabs(c) <= 1e-15 * std::sqrt(a * b))
                if (std::fabs(c) <= 0 * std::sqrt(a * b))
                    continue;
                // 计算 Jacobi 旋转角度(cs, sn)
                const double tau = (b - a) / (2.0 * c);
                const double t = ((tau >= 0.0) ? 1.0 : -1.0) /
                                 (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
                const double cs = 1.0 / std::sqrt(1.0 + t * t);
                const double sn = cs * t;
                // 对 Q 的列 p 与列 r 做旋转
                for (std::size_t i = 0; i < n; ++i) {
                    const double qp = q[i * 3 + p];
                    const double qr = q[i * 3 + r];
                    q[i * 3 + p] = cs * qp - sn * qr;
                    q[i * 3 + r] = sn * qp + cs * qr;
                }
                // 对 V 的列 p 与列 r 做旋转
                for (int i = 0; i < 3; ++i) {
                    const double vp = v[i][p];
                    const double vr = v[i][r];
                    v[i][p] = cs * vp - sn * vr;
                    v[i][r] = sn * vp + cs * vr;
                }
            }
        }
        // 如果 off 已经足够小，则认为收敛
        // if (off <= 1e-14)
        if (off <= 0)
            break;
    }
    // 计算奇异值
    // q的三列的范数(即向量的长度)即为奇异值(降序前)
    for (int j = 0; j < 3; ++j) {
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double qq = q[i * 3 + j];
            s += qq * qq;
        }
        sv[j] = std::sqrt(s);
    }
}

}  // namespace

PlaneFitResult fitPlaneOrthogonal(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& z) {
    // 数据校验
    const std::size_t n = x.size();
    if (y.size() != n || z.size() != n)
        throw std::invalid_argument("x, y, z must have the same size");
    if (n < 3)
        throw std::invalid_argument("at least three points are required");

    // 计算质心
    PlaneFitResult r;
    r.cx = 0.0;
    r.cy = 0.0;
    r.cz = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        r.cx += x[i];
        r.cy += y[i];
        r.cz += z[i];
    }
    const double dn = static_cast<double>(n);
    r.cx /= dn;
    r.cy /= dn;
    r.cz /= dn;

    // 计算相对质心坐标
    std::vector<double> q(n * 3);
    for (std::size_t i = 0; i < n; ++i) {
        q[i * 3 + 0] = x[i] - r.cx;
        q[i * 3 + 1] = y[i] - r.cy;
        q[i * 3 + 2] = z[i] - r.cz;
    }

    // Q = U * diag(sv) * V^T
    double v[3][3];
    double sv[3];
    jacobiSvd3(q, n, v, sv);

    // 最小奇异值对应的右奇异向量即平面法向量
    int kmin = 0;
    for (int j = 1; j < 3; ++j)
        if (sv[j] < sv[kmin])
            kmin = j;
    r.nx = v[0][kmin];
    r.ny = v[1][kmin];
    r.nz = v[2][kmin];

    // 统一符号便于显示：尽量使 nz >= 0
    if (r.nz < 0.0 || (r.nz == 0.0 && r.ny < 0.0)) {
        r.nx = -r.nx;
        r.ny = -r.ny;
        r.nz = -r.nz;
    }

    r.d = -(r.nx * r.cx + r.ny * r.cy + r.nz * r.cz);

    // 计算每个点到平面的正交距离平方和
    r.residual.resize(n);
    r.sse = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double qx = x[i] - r.cx;
        const double qy = y[i] - r.cy;
        const double qz = z[i] - r.cz;
        const double dist = std::fabs(r.nx * qx + r.ny * qy + r.nz * qz);
        r.residual[i] = dist;
        r.sse += dist * dist;
    }

    // 将奇异值按降序排列
    double s[3] = {sv[0], sv[1], sv[2]};
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (s[j] > s[i]) {
                const double tmp = s[i];
                s[i] = s[j];
                s[j] = tmp;
            }
    r.sv1 = s[0];
    r.sv2 = s[1];
    r.sv3 = s[2];

    return r;
}
