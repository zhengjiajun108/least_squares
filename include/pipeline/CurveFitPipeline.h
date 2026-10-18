#ifndef CURVE_FIT_PIPELINE_H
#define CURVE_FIT_PIPELINE_H

#include "core/CurveFit.h"

#include <string>

// 从 inputPath 的 inputSheet 中按表头名读取点集(xHeader 作 x、yHeader 作 y)，
// 做正交距离多项式曲线拟合(自动选阶)，并把汇总结果打印到控制台，不写回文件。
CurveFitResult runCurveFitPipeline(const std::string& inputPath,
                                   const std::string& inputSheet,
                                   const std::string& xHeader,
                                   const std::string& yHeader,
                                   int maxDegree = 5);

#endif
