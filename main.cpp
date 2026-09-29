#include "LineFit.h"
#include "XlsxReader.h"

#include <exception>
#include <iostream>
#include <vector>

int main() {
    std::vector<double> x, y;
    try {
        readPointsFromXlsx("label.xlsx", "L1", "X0", "Y0", x, y);
    } catch (const std::exception& e) {
        std::cerr << "读取 label.xlsx 失败: " << e.what() << "\n";
        return 1;
    }
    std::cout << "points = " << x.size() << "\n";

    LineFitResult ur = uniformLineFit(x, y);
    std::cout << "[Uniform] k = " << ur.k << ", b = " << ur.b
              << ", Xb = " << ur.Xb << ", d = " << ur.d << "\n";
    std::cout << "sse = " << ur.sse << "\n";

    LineFitResult orr = orthogonalLineFit(x, y);
    std::cout << "[Orthogonal] k = " << orr.k << ", b = " << orr.b
              << ", Xb = " << orr.Xb << ", d = " << orr.d << "\n";
    std::cout << "sse = " << orr.sse << "\n";

    const std::vector<std::string> headers = {"x0", "y0", "x1", "y1"};
    const std::vector<std::vector<double>> cols = {ur.xFit, ur.yFit,
                                                   orr.xFit, orr.yFit};
    // try {
    //     writeColumnsToXlsx("label.xlsx", "最小二乘法优化", headers, cols);
    //     std::cout << "已写入 sheet: 最小二乘法优化\n";
    // } catch (const std::exception& e) {
    //     std::cerr << "写入 label.xlsx 失败: " << e.what() << "\n";
    //     return 1;
    // }
    return 0;
}
