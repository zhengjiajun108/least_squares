#include "LineFit.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>

LineFitResult orthogonalLineFit(const std::vector<double>& x,
                                          const std::vector<double>& y) {
    
    //数据检验                                        
    if (x.size() != y.size())
        throw std::invalid_argument("x and y must have the same size");
    const std::size_t n = x.size();
    if (n < 2)
        throw std::invalid_argument("at least two points are required");

    //计算x和y的平均值
    double meanX = 0.0, meanY = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        meanX += x[i];
        meanY += y[i];
    }
    meanX /= static_cast<double>(n);
    meanY /= static_cast<double>(n);

    //Sxx = Xi - meanX从1到n的和
    //Sxy = (Xi - meanX)(Yi - meanY)从1到n的和
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = x[i] - meanX;
        const double dy = y[i] - meanY;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if (sxx + syy == 0.0)
        throw std::invalid_argument("all points are identical; cannot fit");

    //计算k,b
    const double theta = 0.5 * std::atan2(2.0 * sxy, sxx - syy);
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    if (std::fabs(c) < 1e-12)
        throw std::invalid_argument("principal direction is vertical");
    const double k = s / c;
    const double b = meanY - k * meanX;

    //Sii = i的平方从0到(n-1)的和
    //Sim = i * m 从0到(n-1)的和
    //m是原始点在固定直线上的投影的x坐标
    const double onePlusK2 = 1.0 + k * k;
    double Si = 0.0, Sii = 0.0, Sm = 0.0, Sim = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double di = static_cast<double>(i);
        const double m = (x[i] + k * (y[i] - b)) / onePlusK2;
        Si += di;
        Sii += di * di;
        Sm += m;
        Sim += di * m;
    }

    //解得x方向间隔d和起始横坐标Xb
    const double dn = static_cast<double>(n);
    const double denom = dn * Sii - Si * Si;
    if (denom == 0.0)
        throw std::invalid_argument("degenerate index set; cannot fit");

    LineFitResult r;
    r.d = (dn * Sim - Si * Sm) / denom;
    r.Xb = (Sm - r.d * Si) / dn;
    r.k = k;
    r.b = b;

    if (r.d == 0.0)
        throw std::invalid_argument("x spacing is zero; degenerate solution");

    //记录去噪后的点位置
    r.xFit.resize(n);
    r.yFit.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        r.xFit[i] = r.Xb + static_cast<double>(i) * r.d;
        r.yFit[i] = r.k * r.xFit[i] + r.b;
    }

    //计算残差
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
