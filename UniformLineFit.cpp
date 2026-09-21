#include "LineFit.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>

LineFitResult uniformLineFit(const std::vector<double>& x,
                             const std::vector<double>& y) {
    //数据校验
    if (x.size() != y.size())
        throw std::invalid_argument("x and y must have the same size");
    const std::size_t n = x.size();
    if (n < 2)
        throw std::invalid_argument("at least two points are required");

    //求和，方便后续计算
    //例如：Six = i * Xi 从0到(n-1)的和
    double Si = 0.0, Sii = 0.0;
    double Sx = 0.0, Sy = 0.0, Six = 0.0, Siy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double di = static_cast<double>(i);
        Si += di;
        Sii += di * di;
        Sx += x[i];
        Sy += y[i];
        Six += di * x[i];
        Siy += di * y[i];
    }

    const double dn = static_cast<double>(n);
    const double denom = dn * Sii - Si * Si;
    if (denom == 0.0)
        throw std::invalid_argument("degenerate index set; cannot fit");

    //解未知量
    const double beta = (dn * Six - Si * Sx) / denom;
    const double alpha = (Sx - beta * Si) / dn;
    const double delta = (dn * Siy - Si * Sy) / denom;
    const double gamma = (Sy - delta * Si) / dn;

    //记录结果
    LineFitResult r;
    r.Xb = alpha;
    r.d = beta;
    if (beta == 0.0)
        throw std::invalid_argument("x spacing is zero; line is vertical");
    r.k = delta / beta;
    r.b = gamma - r.k * alpha;

    //保存去噪后的点的位置
    r.xFit.resize(n);
    r.yFit.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        r.xFit[i] = r.Xb + static_cast<double>(i) * r.d;
        r.yFit[i] = r.k * r.xFit[i] + r.b;
    }

    //计算残差与残差平方和
    r.residual.resize(n);
    r.sse = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = r.xFit[i] - x[i];
        const double dy = r.yFit[i] - y[i];
        r.residual[i] = std::sqrt(dx * dx + dy * dy);
        r.sse += dx * dx + dy * dy;
    }
    return r;
}
