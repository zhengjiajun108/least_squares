#include "LineFitPipeline.h"
#include "PlaneFitPipeline.h"

#include <exception>
#include <iostream>

int main() {
    try {
        // runPlaneFitPipeline("label.xlsx", "L2", "X6(降噪)", "Y6(降噪)", "Z6(降噪)");
        runLineFitPipeline("label.xlsx", "旋转轴与投影", "X6(降噪)", "Y6(降噪)");
        runLineFitPipeline("label.xlsx", "旋转轴与投影", "X6(未降噪)", "Y6(未降噪)");
        runLineFitPipeline("label.xlsx", "旋转轴与投影", "X6(对比)", "Y6(对比)");
    } catch (const std::exception& e) {
        std::cerr << "失败: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
