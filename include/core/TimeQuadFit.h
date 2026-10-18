#ifndef TIME_QUAD_FIT_H
#define TIME_QUAD_FIT_H

#include <vector>

// 单个坐标轴的分方向二次最小二乘结果
struct AxisQuadFit {
    std::vector<double> coeff;     // [a, b, c] 对应 t^2, t, 1
    std::vector<double> fit;       // 各点拟合值
    std::vector<double> residual;  // 观测 - 拟合
    double sse;                    // 残差平方和
    double rmse;                   // sqrt(sse/N)
};

// 以时间为自变量的分方向二次最小二乘整体结果
struct TimeQuadFitResult {
    double dt;                      // 时间步长(秒)
    std::vector<double> t;          // t_i = (i-1)*dt
    std::vector<AxisQuadFit> axes;  // 每个坐标轴一个结果
    double sse;                     // 各轴 SSE 之和(固定时间残差)
    double rmse;                    // 合成 RMSE = sqrt(各轴 SSE 之和 / N)
    double sseGeo;                  // 原始点到曲线的联合最小几何距离平方和
    double rmseGeo;                 // sqrt(sseGeo / N)
};

// 分方向最小二乘：时间 t_i=(i-1)*dt，对每个坐标轴独立拟合 v(t)=a t^2+b t+c。
// values 每个内层 vector 为一个坐标轴的观测值，需长度一致且行序即时间序(第 1 行为
// t=0)。所有轴共用同一设计矩阵 [t^2, t, 1]，用 SVD 求解。要求 N>=3。
// 拟合后还会计算每个原始点到该参数曲线的最小几何距离平方和(sseGeo)，
// 最近点参数 t 限定在观测时间区间 [0, (N-1)*dt]。
TimeQuadFitResult fitTimeQuadratic(
    const std::vector<std::vector<double>>& values, double dt);

#endif
