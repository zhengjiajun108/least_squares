#ifndef SVD_PLANE_FIT_H
#define SVD_PLANE_FIT_H

#include <vector>

struct PlaneFitResult {
    double nx;   // 单位法向量 x
    double ny;   // 单位法向量 y
    double nz;   // 单位法向量 z
    double d;    // 平面方程 n·p + d = 0
    double cx;   // 质心 x
    double cy;   // 质心 y
    double cz;   // 质心 z
    double sv1;  // 奇异值(降序)
    double sv2;
    double sv3;
    double sse;                     // 正交距离平方和 = sv3^2
    std::vector<double> residual;   // 每个点到平面的绝对距离
};

// 用 SVD(one-sided Jacobi) 对三维点做正交距离平面拟合：
// 最小化各点到平面的正交距离平方和，返回法向量、平面偏移与残差。
PlaneFitResult fitPlaneOrthogonal(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& z);

#endif
