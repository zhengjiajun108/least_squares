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
// 再把两组去噪点写入 inputPath 的 outputSheet(列头 x0,y0,x1,y1)。
LineFitSummary runLineFitPipeline(const std::string& inputPath,
                                  const std::string& inputSheet,
                                  const std::string& inputXHeader,
                                  const std::string& inputYHeader,
                                  const std::string& outputSheet);

#endif
