// 排序工具
#include "utils/SortUtils.h"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>

namespace sortutils {

void sortPoints(std::vector<double>& x,
                std::vector<double>* y,
                std::vector<double>* z,
                SortAxis axis,
                bool descending) {
    // 收集非空列，顺序固定为 x, y, z
    std::vector<std::vector<double>*> cols;
    cols.push_back(&x);
    if (y != nullptr)
        cols.push_back(y);
    if (z != nullptr)
        cols.push_back(z);

    // 校验各非空列长度一致
    const std::size_t n = x.size();
    for (const std::vector<double>* col : cols)
        if (col->size() != n)
            throw std::invalid_argument(
                "sortPoints: all provided columns must have the same size");

    // 定位排序依据列
    std::size_t axisIndex = 0;
    if (axis == SortAxis::Y)
        axisIndex = 1;
    else if (axis == SortAxis::Z)
        axisIndex = 2;
    if (axisIndex >= cols.size())
        throw std::invalid_argument(
            "sortPoints: sort axis column was not provided");
    if (n == 0)
        return;

    const std::vector<double>& keyCol = *cols[axisIndex];

    // 生成下标并按依据列排序
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::stable_sort(idx.begin(), idx.end(),
                     [&keyCol, descending](std::size_t a, std::size_t b) {
                         return descending ? keyCol[a] > keyCol[b]
                                           : keyCol[a] < keyCol[b];
                     });

    // 按同一置换原地重排每一列
    std::vector<double> tmp(n);
    for (std::vector<double>* col : cols) {
        for (std::size_t i = 0; i < n; ++i)
            tmp[i] = (*col)[idx[i]];
        col->swap(tmp);
    }
}

}  // namespace sortutils
