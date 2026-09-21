#ifndef XLSX_READER_H
#define XLSX_READER_H

#include <string>
#include <vector>

// 读取 xlsx 中第一个工作表(xl/worksheets/sheet1.xml)的点数据
// 第一列(A)为 x 坐标，第二列(B)为 y 坐标
void readPointsFromXlsx(const std::string& path,
                        std::vector<double>& x,
                        std::vector<double>& y);

#endif
