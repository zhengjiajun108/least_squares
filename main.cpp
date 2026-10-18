#include "pipeline/CurveFitPipeline.h"
#include "pipeline/TimeQuadFitPipeline.h"

#include <exception>
#include <iostream>

int main() {
    try {
        // runCurveFitPipeline("lab110.xlsx", "L1", "X0", "Y0");
        runTimeQuadFitPipeline("lab110.xlsx", "L1", "X0", "Y0", 0.0055);
    } catch (const std::exception& e) {
        std::cerr << "失败: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
