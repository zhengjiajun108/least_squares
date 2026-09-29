#ifndef UNIFORM_LINE_FIT_H
#define UNIFORM_LINE_FIT_H

#include <vector>

struct LineFitResult {
    double k;   //斜率
    double b;   //截距
    double Xb;  //起始点横坐标
    double d;   //x方向间隔
    std::vector<double> xFit;   //去噪点的x坐标
    std::vector<double> yFit;   //去噪点的y坐标
    std::vector<double> residual;  // 每个原始点到对应去噪点的欧氏距离
    double sse;                    // 残差平方和，即目标函数 D 的最小值
};

// 联合优化函数(残值均匀)
LineFitResult uniformLineFit(const std::vector<double>& x,
                             const std::vector<double>& y);

//正交距离回归(总体最小二乘法，同时考虑X和Y方向残差)
LineFitResult orthogonalLineFit(const std::vector<double>& x,
                             const std::vector<double>& y);

#endif
