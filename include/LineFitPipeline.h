#ifndef LINE_FIT_PIPELINE_H
#define LINE_FIT_PIPELINE_H

#include "LineFit.h"

#include <cstddef>
#include <string>

struct LineFitSummary {
    std::size_t pointCount;
    LineFitResult uniform;
    LineFitResult orthogonal;
};

// 从 inputPath 的 inputSheet 中按表头名(第 1 列 inputXHeader、第 2 列 inputYHeader)
// 读取点数据，分别做联合优化(uniform)与正交距离(orthogonal)拟合，
// 并将点数与两组拟合参数(k, b, Xb, d, sse)打印到控制台，不写回文件。
LineFitSummary runLineFitPipeline(const std::string& inputPath,
                                  const std::string& inputSheet,
                                  const std::string& inputXHeader,
                                  const std::string& inputYHeader);

#endif
