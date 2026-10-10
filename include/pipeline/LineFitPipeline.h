#ifndef LINE_FIT_PIPELINE_H
#define LINE_FIT_PIPELINE_H

#include "core/LineFit.h"

#include <cstddef>
#include <string>
#include <vector>

struct LineFitSummary {
    std::size_t pointCount;
    LineFitResult uniform;
    LineFitResult orthogonal;
};

// 从 inputPath 的 inputSheet 中按表头名(第 1 列 inputXHeader、第 2 列 inputYHeader)读取点集，进行线性拟合，返回拟合结果
// 如果 outputSheet 和 outputHeaders 非空，则将拟合结果写入 inputPath 的 outputSheet 中，表头为 outputHeaders
LineFitSummary runLineFitPipeline(const std::string& inputPath,
                                  const std::string& inputSheet,
                                  const std::string& inputXHeader,
                                  const std::string& inputYHeader,
                                  const std::string& outputSheet = "",
                                  const std::vector<std::string>& outputHeaders = {});

#endif
