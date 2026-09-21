#include "LineFit.h"

#include <iostream>
#include <vector>

int main() {
    std::vector<double> x = {0.1, 1.2, 1.9, 3.1, 3.8};
    std::vector<double> y = {1.2, 2.9, 5.1, 6.8, 9.2};

    LineFitResult r = uniformLineFit(x, y);
    std::cout << "[Uniform] k = " << r.k << ", b = " << r.b
              << ", Xb = " << r.Xb << ", d = " << r.d << "\n";

    // std::cout << "denoised points:\n";
    // for (std::size_t i = 0; i < r.xFit.size(); ++i)
    //     std::cout << "  (" << r.xFit[i] << ", " << r.yFit[i] << ")"
    //               << "  residual = " << r.residual[i] << "\n";
    std::cout << "sse = " << r.sse << "\n";

    r = orthogonalLineFit(x, y);
    std::cout << "[Orthogonal] k = " << r.k << ", b = " << r.b
              << ", Xb = " << r.Xb << ", d = " << r.d << "\n";

    // std::cout << "denoised points:\n";
    // for (std::size_t i = 0; i < r.xFit.size(); ++i)
    //     std::cout << "  (" << r.xFit[i] << ", " << r.yFit[i] << ")"
    //               << "  residual = " << r.residual[i] << "\n";
    std::cout << "sse = " << r.sse << "\n";
    return 0;
}
