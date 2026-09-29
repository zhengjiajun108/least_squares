#ifndef XLSX_READER_H
#define XLSX_READER_H

#include <string>
#include <vector>

// 读取 xlsx 中指定工作表的点数据。
// 先在该工作表中按表头名 xHeader / yHeader 定位两列，再读取其下方所有数据行，
// 分别作为 x / y 坐标。
void readPointsFromXlsx(const std::string& path,
                        const std::string& sheetName,
                        const std::string& xHeader,
                        const std::string& yHeader,
                        std::vector<double>& x,
                        std::vector<double>& y);

// 将若干列数据写入(或覆盖)指定工作表：第 1 行为列头(内联字符串)，
// 之后每一行为各列对应的数据。若工作表不存在则自动创建。
// 原地修改 xlsx 文件，原有其它工作表/图表等保持不变。
void writeColumnsToXlsx(const std::string& path,
                        const std::string& sheetName,
                        const std::vector<std::string>& headers,
                        const std::vector<std::vector<double>>& columns);

#endif
