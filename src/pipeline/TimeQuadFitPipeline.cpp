// 分方向二次最小二乘的操作封装函数
#include "pipeline/TimeQuadFitPipeline.h"

#include "utils/XlsxReader.h"

#include <cstddef>
#include <iostream>
#include <vector>

TimeQuadFitResult runTimeQuadFitPipeline(const std::string& inputPath,
                                         const std::string& inputSheet,
                                         const std::string& xHeader,
                                         const std::string& yHeader,
                                         double dt) {
    std::vector<double> x, y;
    readPointsFromXlsx(inputPath, inputSheet, xHeader, yHeader, x, y);

    const TimeQuadFitResult r = fitTimeQuadratic({x, y}, dt);

    const std::size_t n = r.t.size();
    std::cout << "=== TimeQuadFit: " << inputPath << " / " << inputSheet
              << " (" << xHeader << ", " << yHeader << ") ===\n";
    std::cout << "points = " << n << "\n";
    std::cout << "dt = " << r.dt << " s   t in [0, "
              << (n > 0 ? (n - 1) * r.dt : 0.0) << "] s\n";

    const std::string labels[2] = {xHeader, yHeader};
    for (std::size_t a = 0; a < r.axes.size(); ++a) {
        const AxisQuadFit& ax = r.axes[a];
        const std::string label =
            (a < 2) ? labels[a] : ("axis" + std::to_string(a));
        std::cout << "[axis " << a << "] " << label << "  a(t^2)=" << ax.coeff[0]
                  << "  b(t)=" << ax.coeff[1] << "  c=" << ax.coeff[2]
                  << "  SSE=" << ax.sse << "  RMSE=" << ax.rmse << "\n";
    }
    std::cout << "combined: SSE=" << r.sse << "  RMSE=" << r.rmse << "\n";
    std::cout << "geometric: SSE_geo=" << r.sseGeo
              << "  RMSE_geo=" << r.rmseGeo << "\n";

    return r;
}
