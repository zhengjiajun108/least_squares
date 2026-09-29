#include "LineFitPipeline.h"

#include <exception>
#include <iostream>

int main() {
    try {
        const LineFitSummary s =
            runLineFitPipeline("label.xlsx", "L1", "X0", "Y0", "球心降噪");
        std::cout << "points = " << s.pointCount << "\n";
        std::cout << "[Uniform] k = " << s.uniform.k << ", b = " << s.uniform.b
                  << ", Xb = " << s.uniform.Xb << ", d = " << s.uniform.d
                  << ", sse = " << s.uniform.sse << "\n";
        std::cout << "[Orthogonal] k = " << s.orthogonal.k
                  << ", b = " << s.orthogonal.b << ", Xb = " << s.orthogonal.Xb
                  << ", d = " << s.orthogonal.d
                  << ", sse = " << s.orthogonal.sse << "\n";
        std::cout << "已写入 sheet: 球心降噪\n";
    } catch (const std::exception& e) {
        std::cerr << "失败: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
