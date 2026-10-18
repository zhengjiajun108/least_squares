#ifndef TIME_QUAD_FIT_PIPELINE_H
#define TIME_QUAD_FIT_PIPELINE_H

#include "core/TimeQuadFit.h"

#include <string>

// 从 inputPath 的 inputSheet 中按表头名读取 x 列(xHeader)与 y 列(yHeader)，
// 以 dt(秒)为时间步长(第 1 行为 t=0)做分方向二次最小二乘拟合，
// 并把各轴系数与误差汇总打印到控制台，不写回文件。
TimeQuadFitResult runTimeQuadFitPipeline(const std::string& inputPath,
                                         const std::string& inputSheet,
                                         const std::string& xHeader,
                                         const std::string& yHeader,
                                         double dt);

#endif
