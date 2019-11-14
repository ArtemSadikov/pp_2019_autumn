// Copyright 2019 Sadikov Artem
#include <mpi.h>
#include <gtest/gtest.h>
#include <gtest-mpi-listener.hpp>

#include <vector>

#include "./simple_method.h"

TEST(Simple_Method_Slae, Test_On_Matrix_Size_3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double* matrix = new double[12];
    if (rank == 0) {
        matrix[0] = 4;
        matrix[1] = 1;
        matrix[2] = 1;
        matrix[3] = 9;
        matrix[4] = 1;
        matrix[5] = 6;
        matrix[6] = -1;
        matrix[7] = 10;
        matrix[8] = 1;
        matrix[9] = 2;
        matrix[10] = 5;
        matrix[11] = 20;
    }

    double* actual = new double[3];

        if (rank == 0) {
            double* expected = new double[3];
            expected[0] = 1.03958;
            expected[1] = 2.00833;
            expected[2] = 3.05;
            EXPECT_TRUE(is_equal(actual, expected));
        }
}

/*TEST(Simple_Method_Slae, Test_On_Matrix_Size_10_Test) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double* matrix = new double[10];

    if (rank == 0) {
        matrix = get_rand_matrix(10);
    }

    double* actual = new double[10];

    actual = get_res(matrix, 10, 0.0001);
    if (rank == 0) {
        EXPECT_FALSE(0);
    }
}*/

/* TEST(DISABLED_Simple_Method_Slae, Cant_Get_Res_With_Wrong_Size) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> matrix;

    if (rank == 0) {
        matrix = get_rand_matrix(4);
    }
    if (rank == 0) {
        std::vector<double> act(4);
        ASSERT_ANY_THROW(act = get_res(matrix, 2, 0.001));
    }
}

TEST(DISABLED_Simple_Method_Slae, Test_On_Matrix_Size_4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> matrix(20);
    if (rank == 0) {
        matrix[0] = 10;
        matrix[1] = -1;
        matrix[2] = 2;
        matrix[3] = 0;
        matrix[4] = 6;
        matrix[5] = -1;
        matrix[6] = 11;
        matrix[7] = -1;
        matrix[8] = 3;
        matrix[9] = 25;
        matrix[10] = 2;
        matrix[11] = -1;
        matrix[12] = 10;
        matrix[13] = -1;
        matrix[14] = -11;
        matrix[15] = 0;
        matrix[16] = 3;
        matrix[17] = -1;
        matrix[18] = 8;
        matrix[19] = 15;
    }

    std::vector<double> actual(4);

    actual = get_res(matrix, 4, 0.99);

    std::vector<double> expected(4);
    if (rank == 0) {
        expected[0] = 1.04727;
        expected[1] = 1.7159;
        expected[2] = -0.80522;
        expected[3] = 0.88522;

        EXPECT_TRUE(is_equal(actual, expected));
    }
}

TEST(DISABLED_Simple_Method_Slae, Not_Create_Matrix_If_Size_Less_Then_2) {
    ASSERT_ANY_THROW(std::vector<double> matrix = get_rand_matrix(1));
}

TEST(DISABLED_Simple_Method_Slae, Can_Create_Matrix_On_5_Size) {
    ASSERT_NO_THROW(std::vector<double> matrix = get_rand_matrix(5));
}

TEST(DISABLED_Simple_Method_Slae, Can_Create_Matrix_On_10_Size) {
    ASSERT_NO_THROW(std::vector<double> matrix = get_rand_matrix(10));
} */

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
