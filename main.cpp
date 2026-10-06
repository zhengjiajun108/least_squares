#include "SvdPlaneFit.h"
#include "XlsxReader.h"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    try {
        // 现有读取函数一次只读两列：分别读 (X6,Y6) 与 (X6,Z6) 组合成三维点
        std::vector<double> x, y, x2, z;
        readPointsFromXlsx("label.xlsx", "L2", "X6", "Y6", x, y);
        readPointsFromXlsx("label.xlsx", "L2", "X6", "Z6", x2, z);
        if (x.size() != x2.size())
            throw std::runtime_error("x 列两次读取的行数不一致");

        const PlaneFitResult r = fitPlaneOrthogonal(x, y, z);

        std::cout << "points = " << x.size() << "\n";
        std::cout << "centroid = (" << r.cx << ", " << r.cy << ", " << r.cz
                  << ")\n";
        std::cout << "normal = (" << r.nx << ", " << r.ny << ", " << r.nz
                  << ")\n";
        std::cout << "d = " << r.d << "\n";
        std::cout << "singular values = (" << r.sv1 << ", " << r.sv2 << ", "
                  << r.sv3 << ")\n";
        std::cout << "sse = " << r.sse << "\n";
        for (std::size_t i = 0; i < r.residual.size(); ++i)
            std::cout << "residual[" << i << "] = " << r.residual[i] << "\n";
    } catch (const std::exception& e) {
        std::cerr << "失败: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
