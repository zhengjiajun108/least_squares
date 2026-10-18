// 分方向二次最小二乘(以时间为自变量)
#include "core/TimeQuadFit.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

// one-sided Jacobi SVD：对行主序 N×m 矩阵 V 做右正交化，使 V·J = U·Σ。
// 输出累积旋转矩阵 J(m×m，行主序)与奇异值 sv(m)。旋转后各列范数即奇异值。
void jacobiSvd(std::vector<double>& V, int N, int m, std::vector<double>& J,
               std::vector<double>& sv) {
    J.assign(static_cast<std::size_t>(m) * static_cast<std::size_t>(m), 0.0);
    for (int i = 0; i < m; ++i)
        J[static_cast<std::size_t>(i) * m + i] = 1.0;

    const int kMaxSweeps = 100;
    const double kPairTol = 1e-15;
    const double kOffTol = 1e-14;
    for (int sweep = 0; sweep < kMaxSweeps; ++sweep) {
        double off = 0.0;
        double total = 0.0;
        for (int p = 0; p < m; ++p) {
            for (int r = p + 1; r < m; ++r) {
                double a = 0.0, b = 0.0, c = 0.0;
                for (int i = 0; i < N; ++i) {
                    const double vp = V[static_cast<std::size_t>(i) * m + p];
                    const double vr = V[static_cast<std::size_t>(i) * m + r];
                    a += vp * vp;
                    b += vr * vr;
                    c += vp * vr;
                }
                off += std::fabs(c);
                total += a + b;
                if (std::fabs(c) <= kPairTol * std::sqrt(a * b))
                    continue;
                const double tau = (b - a) / (2.0 * c);
                const double t = ((tau >= 0.0) ? 1.0 : -1.0) /
                                 (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
                const double cs = 1.0 / std::sqrt(1.0 + t * t);
                const double sn = cs * t;
                for (int i = 0; i < N; ++i) {
                    const double vp = V[static_cast<std::size_t>(i) * m + p];
                    const double vr = V[static_cast<std::size_t>(i) * m + r];
                    V[static_cast<std::size_t>(i) * m + p] = cs * vp - sn * vr;
                    V[static_cast<std::size_t>(i) * m + r] = sn * vp + cs * vr;
                }
                for (int i = 0; i < m; ++i) {
                    const double jp = J[static_cast<std::size_t>(i) * m + p];
                    const double jr = J[static_cast<std::size_t>(i) * m + r];
                    J[static_cast<std::size_t>(i) * m + p] = cs * jp - sn * jr;
                    J[static_cast<std::size_t>(i) * m + r] = sn * jp + cs * jr;
                }
            }
        }
        if (off <= kOffTol * total)
            break;
    }
    sv.assign(static_cast<std::size_t>(m), 0.0);
    for (int k = 0; k < m; ++k) {
        double s2 = 0.0;
        for (int i = 0; i < N; ++i) {
            const double v = V[static_cast<std::size_t>(i) * m + k];
            s2 += v * v;
        }
        sv[static_cast<std::size_t>(k)] = std::sqrt(s2);
    }
}

// SVD 最小二乘：解 min ||A a - b||，A 行主序 N×m。近零奇异值用伪逆置 0。
std::vector<double> lstsqSVD(const std::vector<double>& A, int N, int m,
                             const std::vector<double>& b) {
    std::vector<double> V = A;
    std::vector<double> J;
    std::vector<double> sv;
    jacobiSvd(V, N, m, J, sv);

    double smax = 0.0;
    for (double s : sv)
        if (s > smax)
            smax = s;
    const double cutoff = 1e-12 * smax;

    std::vector<double> d(static_cast<std::size_t>(m), 0.0);
    for (int k = 0; k < m; ++k) {
        const double sk = sv[static_cast<std::size_t>(k)];
        if (sk <= cutoff)
            continue;
        double dot = 0.0;
        for (int i = 0; i < N; ++i)
            dot += V[static_cast<std::size_t>(i) * m + k] * b[i];
        d[static_cast<std::size_t>(k)] = dot / (sk * sk);
    }
    std::vector<double> a(static_cast<std::size_t>(m), 0.0);
    for (int i = 0; i < m; ++i) {
        double s = 0.0;
        for (int k = 0; k < m; ++k)
            s += J[static_cast<std::size_t>(i) * m + k] *
                 d[static_cast<std::size_t>(k)];
        a[static_cast<std::size_t>(i)] = s;
    }
    return a;
}

// 求实系数三次方程 c3 t^3 + c2 t^2 + c1 t + c0 = 0 的所有实根(Cardano/三角法)。
// 返回实根个数，根写入 roots[0..count-1]。c3 相对很小时自动降次为二次/一次。
int solveCubicReal(double c3, double c2, double c1, double c0, double roots[3]) {
    const double scale =
        std::max(std::max(std::fabs(c3), std::fabs(c2)),
                 std::max(std::fabs(c1), std::fabs(c0)));
    if (scale == 0.0)
        return 0;

    // 降次：二次或一次
    if (std::fabs(c3) <= 1e-14 * scale) {
        if (std::fabs(c2) <= 1e-14 * scale) {
            if (std::fabs(c1) <= 1e-14 * scale)
                return 0;
            roots[0] = -c0 / c1;
            return 1;
        }
        const double disc = c1 * c1 - 4.0 * c2 * c0;
        if (disc < 0.0)
            return 0;
        const double sq = std::sqrt(disc);
        roots[0] = (-c1 + sq) / (2.0 * c2);
        roots[1] = (-c1 - sq) / (2.0 * c2);
        return 2;
    }

    // 归一化并压缩为 y^3 + p y + q = 0
    const double a = c2 / c3;
    const double b = c1 / c3;
    const double c = c0 / c3;
    const double p = b - a * a / 3.0;
    const double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
    const double half = q / 2.0;
    const double disc = half * half + (p / 3.0) * (p / 3.0) * (p / 3.0);

    if (disc > 1e-14) {
        // 一个实根
        const double sq = std::sqrt(disc);
        roots[0] = std::cbrt(-half + sq) + std::cbrt(-half - sq) - a / 3.0;
        return 1;
    }
    if (disc >= -1e-14) {
        // 重根
        const double u = std::cbrt(-half);
        roots[0] = 2.0 * u - a / 3.0;
        roots[1] = -u - a / 3.0;
        return 2;
    }
    // 三个实根(三角解法)
    const double m = 2.0 * std::sqrt(-p / 3.0);
    double arg = (3.0 * q) / (2.0 * p) * std::sqrt(-3.0 / p);
    if (arg < -1.0)
        arg = -1.0;
    if (arg > 1.0)
        arg = 1.0;
    const double theta = std::acos(arg);
    const double pi = 3.14159265358979323846;
    for (int k = 0; k < 3; ++k)
        roots[k] = m * std::cos((theta - 2.0 * pi * k) / 3.0) - a / 3.0;
    return 3;
}

// 点 Q 到参数曲线 P(t)=A t^2+B t+C 的距离平方
double pointCurveDistSq(const std::vector<double>& A,
                        const std::vector<double>& B,
                        const std::vector<double>& C,
                        const std::vector<double>& Q, double t) {
    double s = 0.0;
    for (std::size_t k = 0; k < A.size(); ++k) {
        const double pk = A[k] * t * t + B[k] * t + C[k];
        const double d = Q[k] - pk;
        s += d * d;
    }
    return s;
}

// 点 Q 到曲线的最小距离平方，最近点参数限定在 [tmin, tmax]。
// 候选为落在区间内的三次方程实根以及两端点。
double nearestDistSq(const std::vector<double>& A, const std::vector<double>& B,
                     const std::vector<double>& C, const std::vector<double>& Q,
                     double tmin, double tmax) {
    double AA = 0.0, AB = 0.0, BB = 0.0, AD = 0.0, BD = 0.0;
    for (std::size_t k = 0; k < A.size(); ++k) {
        const double d = Q[k] - C[k];
        AA += A[k] * A[k];
        AB += A[k] * B[k];
        BB += B[k] * B[k];
        AD += A[k] * d;
        BD += B[k] * d;
    }
    const double c3 = 2.0 * AA;
    const double c2 = 3.0 * AB;
    const double c1 = BB - 2.0 * AD;
    const double c0 = -BD;

    double roots[3];
    const int cnt = solveCubicReal(c3, c2, c1, c0, roots);

    double best = std::numeric_limits<double>::infinity();
    for (int k = 0; k < cnt; ++k) {
        const double t = roots[k];
        if (t < tmin || t > tmax)
            continue;
        const double d2 = pointCurveDistSq(A, B, C, Q, t);
        if (d2 < best)
            best = d2;
    }
    // 区间端点(始终作为候选)
    const double dmin = pointCurveDistSq(A, B, C, Q, tmin);
    if (dmin < best)
        best = dmin;
    const double dmax = pointCurveDistSq(A, B, C, Q, tmax);
    if (dmax < best)
        best = dmax;
    return best;
}

}  // namespace

TimeQuadFitResult fitTimeQuadratic(
    const std::vector<std::vector<double>>& values, double dt) {
    if (values.empty())
        throw std::invalid_argument("no coordinate axes given");
    const std::size_t N = values[0].size();
    if (N < 3)
        throw std::invalid_argument("at least three points are required");
    for (const std::vector<double>& col : values)
        if (col.size() != N)
            throw std::invalid_argument(
                "all coordinate columns must have the same size");

    TimeQuadFitResult r;
    r.dt = dt;
    r.t.resize(N);
    for (std::size_t i = 0; i < N; ++i)
        r.t[i] = static_cast<double>(i) * dt;

    // 设计矩阵 V：列依次为 [t^2, t, 1]，三个坐标轴共用
    std::vector<double> V(N * 3);
    for (std::size_t i = 0; i < N; ++i) {
        const double ti = r.t[i];
        V[i * 3 + 0] = ti * ti;
        V[i * 3 + 1] = ti;
        V[i * 3 + 2] = 1.0;
    }

    const int n = static_cast<int>(N);
    r.sse = 0.0;
    r.axes.reserve(values.size());
    for (const std::vector<double>& col : values) {
        AxisQuadFit ax;
        ax.coeff = lstsqSVD(V, n, 3, col);
        ax.fit.resize(N);
        ax.residual.resize(N);
        ax.sse = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            const double ti = r.t[i];
            const double f = ax.coeff[0] * ti * ti + ax.coeff[1] * ti +
                             ax.coeff[2];
            ax.fit[i] = f;
            ax.residual[i] = col[i] - f;
            ax.sse += ax.residual[i] * ax.residual[i];
        }
        ax.rmse = std::sqrt(ax.sse / static_cast<double>(N));
        r.sse += ax.sse;
        r.axes.push_back(std::move(ax));
    }

    // 合成 RMSE：按点看待各方向残差(每点 N 个坐标)，sqrt(总 SSE / N)
    r.rmse = std::sqrt(r.sse / static_cast<double>(N));

    // 计算原始点到拟合曲线的最小几何距离平方和(最近点参数限定在观测时间区间)
    const std::size_t m = r.axes.size();
    std::vector<double> A(m), B(m), C(m);
    for (std::size_t k = 0; k < m; ++k) {
        A[k] = r.axes[k].coeff[0];
        B[k] = r.axes[k].coeff[1];
        C[k] = r.axes[k].coeff[2];
    }
    const double tmin = 0.0;
    const double tmax = (N > 0) ? static_cast<double>(N - 1) * dt : 0.0;
    r.sseGeo = 0.0;
    std::vector<double> Q(m);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t k = 0; k < m; ++k)
            Q[k] = values[k][i];
        r.sseGeo += nearestDistSq(A, B, C, Q, tmin, tmax);
    }
    r.rmseGeo = std::sqrt(r.sseGeo / static_cast<double>(N));

    return r;
}
