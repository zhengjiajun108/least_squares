#ifndef CURVE_FIT_H
#define CURVE_FIT_H

#include <vector>

struct CurveFitResult {
    int n;                            // 选定的多项式次数
    std::vector<double> coeff;        // 原坐标多项式系数 a0..an(升幂)
    std::vector<double> t;            // 各点最近点(投影点)的横坐标(原坐标)
    std::vector<double> residual;     // 各点到曲线的正交距离
    std::vector<double> sseByDegree;  // n=1..maxDegree 的正交距离平方和
    std::vector<double> aiccByDegree; // n=1..maxDegree 的 AICc
    double sse;                       // 选定模型的 SSE
    double rmse;                      // 选定模型的 RMSE 根均方误差
    double aicc;                      // 选定模型的 AICc
};

// 正交距离多项式曲线拟合(y=p(x))：交替优化最小化各点到曲线的正交距离平方和。
// 先在 n=1..maxDegree 内自动选阶，取 AICc 最小的次数；AICc 按样本数 2N、
// 参数个数 k=n+1 计算。返回选定模型及全部候选阶数的 SSE / AICc。
// verbose=true 时，每完成一个阶数的拟合并算出 SSE/AICc 后打印一行
// (n、SSE、AICc 及该阶的原坐标系数 a0..an)。
CurveFitResult fitCurveOrthogonal(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  int maxDegree = 5,
                                  bool verbose = false);

#endif
