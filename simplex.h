#ifndef SIMPLEX_H_
#define SIMPLEX_H_
#include <vector>
#include <string>
void show_task(std::vector<double> c, std::vector<double> b, double** a);
void show_new_task(std::vector<double> c, std::vector<double> b, double** a);
void set_data(std::vector<double> &c, std::vector<double> &b, double** a, unsigned int &n, unsigned int &m);
double * simplex_method(std::vector<double> c, std::vector<double> b, double ** a);
bool not_reference(double ** arr, unsigned int m);
bool not_optimal(double ** arr, unsigned int m, unsigned int n);
std::vector<std::string> print(double ** arr, unsigned int m, unsigned int n);
void set_new_data(std::vector<double> &c, std::vector<double> &b, double** a);
static unsigned int ri = {0}; // Разрешающая строка
static unsigned int rj = {0}; // Разрешающий столбец
static std::string min_max; // Ищем min/max
static unsigned int chk = {0};
static unsigned int count = {0};
#endif