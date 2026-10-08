#include "LineFitPipeline.h"
#include "PlaneFitPipeline.h"

#include <exception>
#include <iostream>

int main() {
    try {
        // runPlaneFitPipeline("label.xlsx", "L2", "X6(降噪)", "Y6(降噪)", "Z6(降噪)");
        // runLineFitPipeline("label.xlsx", "旋转轴与角速度", "方位角θ(球心降噪)", "俯度角Φ(球心降噪)", "测试", {"方位角θ降噪(球心降噪)", "俯度角Φ降噪(球心降噪)"});
        // runLineFitPipeline("label.xlsx", "旋转轴与角速度", "方位角θ(球心未降噪)", "俯度角Φ(球心未降噪)", "测试", {"方位角θ降噪(球心未降噪)", "俯度角Φ降噪(球心未降噪)"});
        // runLineFitPipeline("label.xlsx", "旋转轴与角速度", "方位角θ(参考)", "俯度角Φ(参考)", "测试", {"方位角θ降噪(参考)", "俯度角Φ降噪(参考)"});
    } catch (const std::exception& e) {
        std::cerr << "失败: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
