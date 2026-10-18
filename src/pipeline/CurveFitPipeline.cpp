// 正交距离多项式曲线拟合的操作封装函数
#include "pipeline/CurveFitPipeline.h"

#include "utils/XlsxReader.h"

#include <cstddef>
#include <iostream>
#include <vector>

CurveFitResult runCurveFitPipeline(const std::string& inputPath,
                                   const std::string& inputSheet,
                                   const std::string& xHeader,
                                   const std::string& yHeader,
                                   int maxDegree) {
    std::vector<double> x, y;
    readPointsFromXlsx(inputPath, inputSheet, xHeader, yHeader, x, y);

    const CurveFitResult r = fitCurveOrthogonal(x, y, maxDegree, true);

    // std::cout << "=== CurveFit: " << inputPath << " / " << inputSheet << " ("
    //           << xHeader << ", " << yHeader << ") ===\n";
    // std::cout << "points = " << r.residual.size() << "\n";

    // // 各候选阶数的 SSE 与 AICc
    // std::cout << "[degree scan]\n";
    // for (std::size_t i = 0; i < r.sseByDegree.size(); ++i) {
    //     std::cout << "  n=" << (i + 1) << "  SSE=" << r.sseByDegree[i]
    //               << "  AICc=" << r.aiccByDegree[i] << "\n";
    // }

    // std::cout << "[selected] n=" << r.n << "\n";
    // std::cout << "  coeff (a0..an) =";
    // for (double c : r.coeff)
    //     std::cout << " " << c;
    // std::cout << "\n";
    // std::cout << "  SSE=" << r.sse << "  RMSE=" << r.rmse
    //           << "  AICc=" << r.aicc << "\n";

    return r;
}
