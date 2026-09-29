// 联合优化和分布优化操作封装函数
#include "LineFitPipeline.h"

#include "XlsxReader.h"

#include <vector>

LineFitSummary runLineFitPipeline(const std::string& inputPath,
                                  const std::string& inputSheet,
                                  const std::string& inputXHeader,
                                  const std::string& inputYHeader,
                                  const std::string& outputSheet) {
    std::vector<double> x, y;
    readPointsFromXlsx(inputPath, inputSheet, inputXHeader, inputYHeader, x, y);

    LineFitSummary summary;
    summary.pointCount = x.size();
    summary.uniform = uniformLineFit(x, y);
    summary.orthogonal = orthogonalLineFit(x, y);

    const std::vector<std::string> headers = {"x0", "y0", "x1", "y1"};
    const std::vector<std::vector<double>> columns = {
        summary.uniform.xFit, summary.uniform.yFit,
        summary.orthogonal.xFit, summary.orthogonal.yFit};
    writeColumnsToXlsx(inputPath, outputSheet, headers, columns);

    return summary;
}
