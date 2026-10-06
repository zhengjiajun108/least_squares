// 联合优化和分布优化操作封装函数
#include "LineFitPipeline.h"

#include "XlsxReader.h"
#include "SortUtils.h"

#include <iostream>
#include <vector>

LineFitSummary runLineFitPipeline(const std::string& inputPath,
                                  const std::string& inputSheet,
                                  const std::string& inputXHeader,
                                  const std::string& inputYHeader) {
    std::vector<double> x, y;
    readPointsFromXlsx(inputPath, inputSheet, inputXHeader, inputYHeader, x, y);
    
    // 按 x 升序排序，y 列随 x 列一起重排
    sortutils::sortPoints(x, &y);

    LineFitSummary summary;
    summary.pointCount = x.size();
    // summary.uniform = uniformLineFit(x, y);
    // summary.orthogonal = orthogonalLineFit(x, y);
    summary.orthogonal = orthogonalFitResidual(x, y);

    std::cout << "=== LineFit: " << inputPath << " / " << inputSheet << " ("
              << inputXHeader << ", " << inputYHeader << ") ===\n";
    std::cout << "points = " << summary.pointCount << "\n";
    // std::cout << "[uniform]    k=" << summary.uniform.k
    //           << "  b=" << summary.uniform.b << "  Xb=" << summary.uniform.Xb
    //           << "  d=" << summary.uniform.d << "  sse=" << summary.uniform.sse
    //           << "\n";
    std::cout << "[orthogonal] k=" << summary.orthogonal.k
              << "  b=" << summary.orthogonal.b
              << "  Xb=" << summary.orthogonal.Xb
              << "  d=" << summary.orthogonal.d
              << "  sse=" << summary.orthogonal.sse << "\n";

    return summary;
}
