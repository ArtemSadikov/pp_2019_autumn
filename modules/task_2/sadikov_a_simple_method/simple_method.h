// Copyright 2019 Sadikov Artem
#ifndef MODULES_TASK_2_SADIKOV_A_SIMPLE_METHOD_SIMPLE_METHOD_H_
#define MODULES_TASK_2_SADIKOV_A_SIMPLE_METHOD_SIMPLE_METHOD_H_

#include <vector>

double* get_rand_matrix(int size);
double* solve_simple(std::vector<double> delta_a, double* x,
                 double error, int size, int rank, int row_count,
                 int size_proc);
bool is_equal(double* x, double* y);
double* get_res(double* matrix, int size, double error);

#endif  // MODULES_TASK_2_SADIKOV_A_SIMPLE_METHOD_SIMPLE_METHOD_H_
