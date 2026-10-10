#ifndef PLANE_FIT_PIPELINE_H
#define PLANE_FIT_PIPELINE_H

#include "core/PlaneFit.h"

#include <string>

// 从 inputPath 的 inputSheet 中按表头名读取三维点：
// 用 (xHeader, yHeader) 读 x, y，用 (xHeader, zHeader) 读 x, z，拼成 (x, y, z)；
// 做正交距离平面拟合并将汇总参数打印到控制台，不写回文件。
PlaneFitResult runPlaneFitPipeline(const std::string& inputPath,
                                   const std::string& inputSheet,
                                   const std::string& xHeader,
                                   const std::string& yHeader,
                                   const std::string& zHeader);

#endif
