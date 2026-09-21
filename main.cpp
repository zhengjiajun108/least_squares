#include "LineFit.h"
#include "XlsxReader.h"

#include <exception>
#include <iostream>
#include <vector>

int main() {
    std::vector<double> x, y;
    try {
        readPointsFromXlsx("points.xlsx", x, y);
    } catch (const std::exception& e) {
        std::cerr << "读取 points.xlsx 失败: " << e.what() << "\n";
        return 1;
    }

    LineFitResult r = uniformLineFit(x, y);
    std::cout << "[Uniform] k = " << r.k << ", b = " << r.b
              << ", Xb = " << r.Xb << ", d = " << r.d << "\n";
    // for(std::size_t i = 0; i < r.xFit.size(); ++i) {
    //     std::cout << "xFit[" << i << "] = " << r.xFit[i]
    //               << ", yFit[" << i << "] = " << r.yFit[i]
    //               << ", residual[" << i << "] = " << r.residual[i] << "\n";
    // }
    std::cout << "sse = " << r.sse << "\n";

    r = orthogonalLineFit(x, y);
    std::cout << "[Orthogonal] k = " << r.k << ", b = " << r.b
              << ", Xb = " << r.Xb << ", d = " << r.d << "\n";
    std::cout << "sse = " << r.sse << "\n";
    return 0;
}
