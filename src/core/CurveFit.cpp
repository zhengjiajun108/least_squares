// 正交距离多项式曲线拟合：交替优化 + AICc 选阶
#include "core/CurveFit.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

// 多项式求值(霍纳法)：a 为升幂系数 a0 + a1 x + ... + an x^n
double polyVal(const std::vector<double>& a, double x) {
    double r = 0.0;
    for (std::size_t k = a.size(); k-- > 0;)
        r = r * x + a[k];
    return r;
}

// 一阶导数求和：sum_{k>=1} k a_k x^(k-1)
double polyDer1(const std::vector<double>& a, double x) {
    double r = 0.0;
    for (std::size_t k = a.size(); k-- > 1;)
        r = r * x + static_cast<double>(k) * a[k];
    return r;
}

// 二阶导数求和：sum_{k>=2} k(k-1) a_k x^(k-2)
double polyDer2(const std::vector<double>& a, double x) {
    double r = 0.0;
    for (std::size_t k = a.size(); k-- > 2;)
        r = r * x + static_cast<double>(k) * static_cast<double>(k - 1) * a[k];
    return r;
}

// 点 (xi,yi) 到曲线在横坐标 t 处点的距离平方
double dist2(double xi, double yi, double t, const std::vector<double>& a) {
    const double dx = xi - t;
    const double dy = yi - polyVal(a, t);
    return dx * dx + dy * dy;
}

// one-sided Jacobi SVD：对行主序 N×m 矩阵 V 做右正交化，使 V·J = U·Σ。
// 输出累积旋转矩阵 J(m×m，行主序)与奇异值 sv(m)。求最小二乘只需右奇异方向，
// 故不显式构造 U。旋转后各列的范数即为奇异值。
void jacobiSvd(std::vector<double>& V, int N, int m, std::vector<double>& J,
               std::vector<double>& sv) {
    // 初始化 J 为单位矩阵
    J.assign(static_cast<std::size_t>(m) * static_cast<std::size_t>(m), 0.0);
    for (int i = 0; i < m; ++i)
        J[static_cast<std::size_t>(i) * m + i] = 1.0;

    const int kMaxSweeps = 100;
    const double kPairTol = 1e-15;  // 列对内积相对阈值
    const double kOffTol = 1e-14;   // 整轮非对角量相对阈值
    for (int sweep = 0; sweep < kMaxSweeps; ++sweep) {
        double off = 0.0;    // 本轮的列对非对角量累加
        double total = 0.0;  // 各列范数平方和，用于相对判据
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
                // 已正交则跳过
                if (std::fabs(c) <= kPairTol * std::sqrt(a * b))
                    continue;
                // Jacobi 旋转(cs, sn)
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
    // 旋转后各列范数即奇异值
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

// SVD 最小二乘：解 min ||A a - b||，A 行主序 N×m。用伪逆处理秩亏/病态。
// 由 V·J = U·Σ 得 a = J·d，其中 d_k = (col_k·b) / σ_k²，近零奇异值对应的 d_k 置 0。
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

    // d_k = (col_k · b) / σ_k²，近零奇异值置 0
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
    // a = J · d
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

// 构造 Vandermonde 矩阵 V[i][j] = t[i]^j，j=0..n，行主序
std::vector<double> vandermonde(const std::vector<double>& t, int n) {
    const std::size_t N = t.size();
    std::vector<double> V(N * static_cast<std::size_t>(n + 1));
    for (std::size_t i = 0; i < N; ++i) {
        double p = 1.0;
        for (int j = 0; j <= n; ++j) {
            V[i * static_cast<std::size_t>(n + 1) + j] = p;
            p *= t[i];
        }
    }
    return V;
}

// 最小二乘拟合多项式(升幂)，自变量 t、因变量 yv
std::vector<double> lsPolyFit(const std::vector<double>& t,
                              const std::vector<double>& yv, int n) {
    const std::vector<double> V =
        vandermonde(t, n);
    return lstsqSVD(V, static_cast<int>(t.size()), n + 1, yv);
}

// 牛顿法求距离函数驻点：G(t) = (t-xi) + (p(t)-yi) p'(t) = 0
double newtonNearest(const std::vector<double>& a, double xi, double yi,
                     double t0, double lo, double hi) {
    double t = t0;
    if (t < lo)
        t = lo;
    if (t > hi)
        t = hi;
    for (int k = 0; k < 60; ++k) {
        const double p = polyVal(a, t);
        const double dp = polyDer1(a, t);
        const double ddp = polyDer2(a, t);
        const double G = (t - xi) + (p - yi) * dp;
        const double Gp = 1.0 + dp * dp + (p - yi) * ddp;
        if (std::fabs(Gp) < 1e-14)
            break;
        double tn = t - G / Gp;
        if (tn < lo)
            tn = lo;
        if (tn > hi)
            tn = hi;
        if (std::fabs(tn - t) < 1e-14) {
            t = tn;
            break;
        }
        t = tn;
    }
    return t;
}

// 求点 (xi,yi) 到曲线的最近点横坐标：网格扫描 + 牛顿精修的兜底策略，
// 同时保留上一轮 t 的牛顿结果，选距离更小者，避免落入较远的局部极小。
double nearestT(const std::vector<double>& a, double xi, double yi,
                double tPrev, double lo, double hi) {
    const int M = 200;
    double bestT = lo;
    double bestD = dist2(xi, yi, lo, a);
    for (int j = 1; j <= M; ++j) {
        const double t = lo + (hi - lo) * static_cast<double>(j) / M;
        const double d = dist2(xi, yi, t, a);
        if (d < bestD) {
            bestD = d;
            bestT = t;
        }
    }
    // 候选一：从网格最优点做牛顿精修
    double tr = newtonNearest(a, xi, yi, bestT, lo, hi);
    double dr = dist2(xi, yi, tr, a);
    if (bestD < dr) {
        dr = bestD;
        tr = bestT;
    }
    // 候选二：从上一轮 t 做牛顿精修
    const double tp = newtonNearest(a, xi, yi, tPrev, lo, hi);
    const double dp = dist2(xi, yi, tp, a);
    if (dp < dr) {
        dr = dp;
        tr = tp;
    }
    return tr;
}

// 组合数 C(n,k)
double binom(int n, int k) {
    double r = 1.0;
    for (int i = 0; i < k; ++i)
        r = r * static_cast<double>(n - i) / static_cast<double>(i + 1);
    return r;
}

// 把标准化坐标下的系数 a 换算回原坐标系数：y = sum b_j x^j
std::vector<double> toOriginalCoeffs(const std::vector<double>& a, double mx,
                                     double my, double s) {
    const int n = static_cast<int>(a.size()) - 1;
    std::vector<double> b(static_cast<std::size_t>(n + 1), 0.0);
    for (int k = 0; k <= n; ++k) {
        const double ck = a[k] * std::pow(s, 1.0 - static_cast<double>(k));
        for (int j = 0; j <= k; ++j)
            b[j] += ck * binom(k, j) * std::pow(-mx, static_cast<double>(k - j));
    }
    b[0] += my;
    return b;
}

// 单个阶数的拟合结果(便于在主函数中汇总)
struct OneDegree {
    bool ok;
    double sse;
    std::vector<double> coeff;     // 原坐标系数
    std::vector<double> tOrig;     // 原坐标投影点横坐标
    std::vector<double> residual;  // 正交距离
};

// 对固定阶数 n 做一次交替优化
OneDegree fitOneDegree(const std::vector<double>& xs,
                       const std::vector<double>& ys, int n, double mx,
                       double my, double s) {
    OneDegree out;
    out.ok = false;
    const std::size_t N = xs.size();

    // 初值：普通最小二乘多项式拟合，投影参数取 x
    std::vector<double> a = lsPolyFit(xs, ys, n);
    std::vector<double> t = xs;

    // 投影参数限幅范围，防止最近点跑到外推区
    double mn = xs[0], mxx = xs[0];
    for (double v : xs) {
        if (v < mn)
            mn = v;
        if (v > mxx)
            mxx = v;
    }
    const double R = mxx - mn;
    const double lo = mn - 2.0 * R;
    const double hi = mxx + 2.0 * R;

    const int kMaxOuter = 100;
    const double kTol = 1e-10;
    double Jold = std::numeric_limits<double>::infinity();
    double J = Jold;

    for (int iter = 0; iter < kMaxOuter; ++iter) {
        // 固定 a，更新投影参数 t
        for (std::size_t i = 0; i < N; ++i)
            t[i] = nearestT(a, xs[i], ys[i], t[i], lo, hi);
        // 固定 t，用最小二乘更新 a
        a = lsPolyFit(t, ys, n);
        // 目标函数
        J = 0.0;
        for (std::size_t i = 0; i < N; ++i)
            J += dist2(xs[i], ys[i], t[i], a);
        if (std::fabs(Jold - J) < kTol * (1.0 + Jold))
            break;
        Jold = J;
    }

    // 换算回原坐标：系数与正交距离
    out.coeff = toOriginalCoeffs(a, mx, my, s);
    out.tOrig.resize(N);
    out.residual.resize(N);
    out.sse = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        out.tOrig[i] = t[i] * s + mx;
        const double dx = xs[i] - t[i];
        const double dy = ys[i] - polyVal(a, t[i]);
        const double d = s * std::sqrt(dx * dx + dy * dy);
        out.residual[i] = d;
        out.sse += d * d;
    }
    out.ok = true;
    return out;
}

}  // namespace

CurveFitResult fitCurveOrthogonal(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  int maxDegree, bool verbose) {
    if (x.size() != y.size())
        throw std::invalid_argument("x and y must have the same size");
    const std::size_t N = x.size();
    if (N < 2)
        throw std::invalid_argument("at least two points are required");

    // 有效扫描上限：不超过 N-1
    int nMax = maxDegree;
    if (nMax > static_cast<int>(N) - 1)
        nMax = static_cast<int>(N) - 1;
    if (nMax < 1)
        throw std::invalid_argument("not enough points for curve fitting");

    // 标准化：中心化 + 统一缩放(保持欧氏距离比例)
    double mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        mx += x[i];
        my += y[i];
    }
    mx /= static_cast<double>(N);
    my /= static_cast<double>(N);
    double vx = 0.0, vy = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        vx += dx * dx;
        vy += dy * dy;
    }
    vx /= static_cast<double>(N);
    vy /= static_cast<double>(N);
    double s = std::sqrt(0.5 * (vx + vy));
    if (s <= 0.0)
        throw std::invalid_argument("all points are identical; cannot fit");

    std::vector<double> xs(N), ys(N);
    for (std::size_t i = 0; i < N; ++i) {
        xs[i] = (x[i] - mx) / s;
        ys[i] = (y[i] - my) / s;
    }

    // 各阶拟合并按 AICc 选最优
    CurveFitResult best;
    best.sse = std::numeric_limits<double>::infinity();
    best.aicc = std::numeric_limits<double>::infinity();
    best.sseByDegree.reserve(static_cast<std::size_t>(nMax));
    best.aiccByDegree.reserve(static_cast<std::size_t>(nMax));

    bool hasBest = false;
    for (int n = 1; n <= nMax; ++n) {
        const OneDegree od = fitOneDegree(xs, ys, n, mx, my, s);
        const double sse = od.sse;
        const int k = n + 1;  // 有效参数个数
        const double sseSafe = (sse > 0.0) ? sse : 1e-300;
        const double aic =
            2.0 * static_cast<double>(N) * std::log(sseSafe /
                                                     static_cast<double>(N)) +
            2.0 * static_cast<double>(k);
        const double denom =
            2.0 * static_cast<double>(N) - static_cast<double>(k) - 1.0;
        const double aicc =
            aic + 2.0 * static_cast<double>(k) *
                      (static_cast<double>(k) + 1.0) / denom;

        // 每个阶数拟合完成后按需打印
        if (verbose) {
            std::cout << "[CurveFit scan] n=" << n << "  SSE=" << sse
                      << "  AICc=" << aicc << "  coeff=[";
            for (std::size_t j = 0; j < od.coeff.size(); ++j) {
                if (j > 0)
                    std::cout << ", ";
                std::cout << od.coeff[j];
            }
            std::cout << "]\n";
        }

        best.sseByDegree.push_back(sse);
        best.aiccByDegree.push_back(aicc);

        if (!hasBest || aicc < best.aicc) {
            hasBest = true;
            best.n = n;
            best.coeff = od.coeff;
            best.t = od.tOrig;
            best.residual = od.residual;
            best.sse = sse;
            best.aicc = aicc;
        }
    }

    best.rmse = std::sqrt(best.sse / static_cast<double>(N));
    return best;
}
