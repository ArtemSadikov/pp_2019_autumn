// Copyright 2019 Sadikov Artem
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <cmath>
#include <limits>
#include "../../../modules/task_2/sadikov_a_simple_method/simple_method.h"

std::vector<double> get_rand_matrix(int size) {
    if (size < 2)
        throw "ERR";

    std::mt19937 gen;
    gen.seed(static_cast<double>(time(NULL)));

    std::vector<double> matrix(size * (size + 1));
    for (int i = 0; i < size * (size + 1); i++) {
        matrix[i] = gen() % 20;
    }

    return matrix;
}

bool is_equal(std::vector<double> x, std::vector<double> y) {
    for (int i = 0; i < static_cast<int>(x.size()); i++) {
        if (!(std::fabs(x[i] - y[i]) < 1e-4))
            return false;
    }
    return true;
}

std::vector<double> solve_simple(std::vector<double> delta_a, std::vector<double> x,
                 double error, int size, int rank, int row_count,
                 int size_proc) {
    std::vector<double> x_old;
    int iter = 0, core;
    double norm = 0;
    std::vector<int> sendcounts, displs;

    sendcounts.resize(size_proc);
    displs.resize(size_proc);

    MPI_Scan(&row_count, &core, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    core -= row_count;

    MPI_Allgather(&row_count, 1, MPI_INT, &sendcounts[0], 1, MPI_INT,
                 MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    displs[0] = 0;
    for (int i = 1; i < size_proc; i++) {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    x_old.resize(size);
    do {
        iter++;
        norm = x[0];

        for (int i = 0; i < size; i++) { x_old[i] = x[i]; }

        double sum;

        for (int i = 0; i < row_count; i++) {
            sum = 0.0;
            for (int j = 0; j < i + core; j++) {
                sum += delta_a[i * (size + 1) + j] * x_old[j];
            }

            for (int j = 1 + core + i; j < size; j++) {
                sum += delta_a[i * (size + 1) + j] * x_old[j];
            }

            x[i + core] = static_cast<double>(delta_a[i * (size + 1) + size]
                         - sum) /
                         static_cast<double>(delta_a[i * (size + 1) + i + core]);
        }

        MPI_Allgatherv(&x[0] + core, row_count, MPI_DOUBLE,
                        &x[0], &sendcounts[0], &displs[0], MPI_DOUBLE,
                        MPI_COMM_WORLD);

        if (rank == 0) {
            std::vector<double> val(size);
            for (int i = 0; i < size; i++) {
                val[i] = fabs(x[i] - x_old[i]);
            }
            double mx = val[0];
            for (int i = 0; i < size; i++) {
                if (mx < val[i])
                    mx = val[i];
            }
            norm = mx;
        }

        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } while (error < norm);

    return x;
}

std::vector<double> get_res(std::vector<double> matrix, int size, double error) {
    if (size * (size + 1) != static_cast<int>(matrix.size()))
        throw "WRONG";

    int row_count, SIZE, size_proc, rank;
    std::vector<double> x;
    std::vector<double> delta_a;
    std::vector<int> sendcounts, displs;
    MPI_Comm_size(MPI_COMM_WORLD, &size_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /*if (size_proc > size)
        throw "ERR";*/

    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    x.resize(size);
    if (rank == 0) {
        for (int i = 0; i < size; i++) { x[i] = 0.0; }
    }
    MPI_Bcast(&x[0], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    row_count = (size / size_proc) +
                ((size % size_proc) > rank ? 1 : 0);

    delta_a.resize((size + 1) * row_count);
    displs.resize(size_proc);
    sendcounts.resize(size_proc);

    SIZE = (size + 1) * row_count;

    MPI_Gather(&SIZE, 1, MPI_INT, &sendcounts[0], 1, MPI_INT, 0,
                MPI_COMM_WORLD);

    displs[0] = 0;
    for (int i = 1; i < size_proc; i++) {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    MPI_Scatterv(&matrix[0], &sendcounts[0], &displs[0], MPI_DOUBLE,
                 &delta_a[0], SIZE, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    x = solve_simple(delta_a, x, error, size, rank, row_count, size_proc);

    return x;
}
