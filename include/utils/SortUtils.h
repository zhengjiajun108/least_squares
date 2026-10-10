#ifndef SORT_UTILS_H
#define SORT_UTILS_H

#include <vector>

namespace sortutils {

enum class SortAxis { X, Y, Z };

// 按指定轴就地排序：所有非空列按同一顺序重排，保持点对应关系。
// y、z 非必须，缺省传 nullptr；axis 指定依据列(该列必须存在且非空)；
// descending=true 时降序。各非空列长度必须一致，否则抛 std::invalid_argument。
void sortPoints(std::vector<double>& x,
                std::vector<double>* y = nullptr,
                std::vector<double>* z = nullptr,
                SortAxis axis = SortAxis::X,
                bool descending = false);

}  // namespace sortutils

#endif
