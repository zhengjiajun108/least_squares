// 平面拟合操作封装函数
#include "PlaneFitPipeline.h"

#include "XlsxReader.h"

#include <iostream>
#include <stdexcept>
#include <vector>

PlaneFitResult runPlaneFitPipeline(const std::string& inputPath,
                                   const std::string& inputSheet,
                                   const std::string& xHeader,
                                   const std::string& yHeader,
                                   const std::string& zHeader) {
    // 现有读取函数一次只读两列：分别读 (x, y) 与 (x, z) 组合成三维点
    std::vector<double> x, y, x2, z;
    readPointsFromXlsx(inputPath, inputSheet, xHeader, yHeader, x, y);
    readPointsFromXlsx(inputPath, inputSheet, xHeader, zHeader, x2, z);
    if (x.size() != x2.size())
        throw std::runtime_error("x 列两次读取的行数不一致");

    const PlaneFitResult r = fitPlaneOrthogonal(x, y, z);

    std::cout << "=== PlaneFit: " << inputPath << " / " << inputSheet << " ("
              << xHeader << ", " << yHeader << ", " << zHeader << ") ===\n";
    std::cout << "points = " << x.size() << "\n";
    std::cout << "centroid = (" << r.cx << ", " << r.cy << ", " << r.cz
              << ")\n";
    std::cout << "normal = (" << r.nx << ", " << r.ny << ", " << r.nz << ")\n";
    std::cout << "d = " << r.d << "\n";
    std::cout << "singular values = (" << r.sv1 << ", " << r.sv2 << ", "
              << r.sv3 << ")\n";
    std::cout << "sse = " << r.sse << "\n";

    return r;
}
