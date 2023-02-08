#include "../src//SquareMatrix.h"
#include "doctest.h"

// Checks if two vectors have the same size and same elements
void compare_vectors(vector<vector<double>> a,
                     vector<vector<double>> b,
                     const double tolerance = 0.001) {
  // allow for a 0.1% error
  for (int i = 0; i < a.size(); i++) {
    CHECK_EQ(a[i].size(), b[i].size());
    for (int j = 0; j < a[i].size(); j++) {
      REQUIRE(a[i][j] - b[i][j] == doctest::Approx(0).epsilon(tolerance));
    }
  }
}

TEST_CASE("Test constructor") {
  SUBCASE("2D matrix as parameter") {
    CHECK_THROWS_AS(SquareMatrix A({{1, 1}}), std::exception);
    CHECK_THROWS_AS(SquareMatrix A({{1, 1}, {1}}), std::exception);
    CHECK_THROWS_AS(SquareMatrix A({}), std::exception);
    CHECK_NOTHROW(SquareMatrix A({{1, 2}, {3, 4}}));
    CHECK_NOTHROW(SquareMatrix A({{1}}));
  }
  SUBCASE("Augmented matrix as parameter") {
    CHECK_NOTHROW(SquareMatrix A({{1, 2, 3}, {3, 4, 5}}, 1));
    CHECK_NOTHROW(SquareMatrix A({{1, 2}, {3, 4}}, 1));
  }
}

TEST_CASE("Test determinant") {
  SUBCASE("1x1 matrix") {
    SquareMatrix A({{1}});
    CHECK_EQ(A.det(), 1);
  }
  SUBCASE("2x2 identity matrix") {
    SquareMatrix A({{1, 0}, {0, 1}});
    CHECK_EQ(A.det(), 1);
  }
  SUBCASE("2x2 zero matrix") {
    SquareMatrix A({{0, 0}, {0, 0}});
    CHECK_EQ(A.det(), 0);
  }
  SUBCASE("3x3 non-singular matrix") {
    SquareMatrix A({{25, 125, 35}, {3, 4, 1}, {0, 1, 6}});
    CHECK_EQ(A.det(), -1570);
  }
  SUBCASE("3x3 singular matrix") {
    SquareMatrix A({{25, 125, 35}, {50, 250, 70}, {0, 1, 6}});
    CHECK_EQ(A.det(), 0);
  }
  SUBCASE("5x5 non-singular matrix") {
    SquareMatrix A({{5, 2, 1, 4, 6},
                    {9, 4, 2, 5, 2},
                    {11, 5, 7, 3, 9},
                    {5, 6, 6, 7, 2},
                    {7, 5, 9, 3, 3}});
    CHECK_EQ(A.det(), -2004);
  }
}

TEST_CASE("Check for diagonal dominance") {
  SUBCASE("3x4 DD augmented matrix") {
    SquareMatrix A({{-22, 3, 8, 27}, {3, -5, 2, -14}, {3, 3, 21, 48}}, true);
    CHECK_EQ(A.is_diagonally_dominant(), true);
  }

  SUBCASE("3x4 non-DD augmented matrix") {
    SquareMatrix A({{-2, 3, 8, 27}, {3, -5, 2, -14}, {3, 3, 21, 48}}, true);
    CHECK_EQ(A.is_diagonally_dominant(), 0);
  }

  SUBCASE("3x3 non-DD matrix") {
    SquareMatrix A({{-22, 3, 8}, {3, -5, 2}, {3, 3, 21}}, true);
    CHECK_EQ(A.is_diagonally_dominant(), true);
  }
}

TEST_CASE("Test Gauss-Jacobi and Gauss-Seidel") {
  SUBCASE("Test input parameters") {
    SquareMatrix A({{4, -1, -1, 3}, {-2, 6, 1, 9}, {-1, 1, 7, -6}}, true);
    SUBCASE("Invalid initial_approx") {
      CHECK_THROWS_WITH(
          A.solve_approx(1, {0, 0, 0, 0}),
          "Size of approximation array must be equal to number of variables");
    }

    SUBCASE("Invalid number of iterations places") {
      CHECK_THROWS_WITH(A.solve_approx(1, {0, 0, 0}, -1),
                        "Number of iterations must be a positive value");
    }

    SUBCASE("Invalid decimal places") {
      CHECK_THROWS_WITH(A.solve_approx(1, {0, 0, 0}, 5, 0),
                        "Number of decimal places must be a positive value");
    }

    SUBCASE("Matrix is not augmented") {
      SquareMatrix A({{4, -1, -1}, {-2, 6, 1}, {-1, 1, 7}}, false);
      CHECK_THROWS_WITH(A.solve_approx(1, {0, 0, 0}), "Matrix must augmented.");
    }

    SUBCASE("Matrix has 0s in its leading diagonal") {
      SquareMatrix A({{0, -1, -1, 3}, {-2, 6, 1, 9}, {-1, 1, 7, -6}}, true);
      CHECK_THROWS_WITH(
          A.solve_approx(1, {0, 0, 0}),
          "Leading diagonal elements of matrix must be non-zero.");
    }
  }

  SUBCASE("Test values calculated") {
    SUBCASE("3x3 using Seidel") {
      // question 3 notion
      SquareMatrix A({{6, -1, -5, -19}, {-1, -7, -1, -22}, {-2, 3, 8, 27}},
                     true);
      vector<vector<double>> table = A.solve_approx(1, {0, 0, 0}, 3);
      const vector<vector<double>> expectedAnswer = {{0, 0, 0},
                                                     {-3.1667, 3.5952, 1.2351},
                                                     {-1.5382, 3.1862, 1.7956},
                                                     {-1.1393, 3.0491, 1.9468}};
      compare_vectors(table, expectedAnswer);
    }
    SUBCASE("3x3 using Jacobi") {
      SquareMatrix A({{4, -1, -1, 3}, {-2, 6, 1, 9}, {-1, 1, 7, -6}}, true);

      vector<vector<double>> table = A.solve_approx(0, {0, 0, 0}, 3);
      const vector<vector<double>> expectedAnswer = {{0, 0, 0},
                                                     {0.75, 1.5, -0.857},
                                                     {0.911, 1.893, -0.964},
                                                     {0.982, 1.964, -0.997}};
      // A.calc_cout();
      compare_vectors(table, expectedAnswer);
    }
    SUBCASE("4x4 using Seidel") {
      SquareMatrix A({{10, 1, 2, 3, 30},
                      {1, 15, 2, -5, 17},
                      {0, 1, 20, 3, 74},
                      {3, -10, -1, 25, 80}},
                     true);
      vector<vector<double>> table = A.solve_approx(1, {0, 0, 0, 0}, 3);
      const vector<vector<double>> expectedAnswer = {
          {0, 0, 0, 0},
          {3.0000, 0.9333, 3.6533, 3.3595},
          {1.1682, 1.6882, 3.1117, 3.8596},
          {1.0510, 1.9349, 3.0243, 3.9688}};
      // A.calc_cout();
      compare_vectors(expectedAnswer, table);
    }
  }
}

TEST_CASE("Testing") {
  SUBCASE("2015-1-sp qu 1b") {
    SquareMatrix A({{52, 0, 38, 24}, {-1, 19, 0, 58}, {-2, 3, 8, 27}}, true);
    A.solve_approx(1, {0, 0, 0}, 10, 2);
    // A.calc_cout();
  }
  SUBCASE("2015-1-sp qu 1b") {
    SquareMatrix A(
        {
            {24, 14, 0, 38},
            {0, -45, -10, -55},
            {0, -5, -10, -15},
        },
        true);
    A.solve_approx(0, {0, 0, 0}, 10, 4);
    // A.calc_cout();
  }
}

TEST_CASE("Tranpose matrix") {
  SUBCASE("4x4 matrix") {
    SquareMatrix A({{5, 6, -1}, {1, 4, 2}, {1, -2, 5}});
    A.transpose();
    vector<vector<double>> expected = {{5, 1, 1}, {6, 4, -2}, {-1, 2, 5}};
    // std::cout << A.stringify();
    compare_vectors(expected, A.get_vec());
  }
}

TEST_CASE("Test matrix inverse using Leibniz") {
  SUBCASE("3x3 invertible matrix") {
    SquareMatrix A({{5, 6, -1}, {1, 4, 2}, {1, -2, 5}});
    A.leb_inv();
    // A.calc_cout();
    vector<vector<double>> expected = {
        {double(24) / 108, double(-28) / 108, double(16) / 108},
        {double(-3) / 108, double(26) / 108, double(-11) / 108},
        {double(-6) / 108, double(16) / 108, double(14) / 108}};
    compare_vectors(expected, A.get_vec());
  }

  SUBCASE("3x3 non-invertible matrix") {
    SquareMatrix A({{5, 6, -1}, {5, 6, -1}, {1, -2, 5}});
    A.leb_inv();
    A.calc_cout();
    compare_vectors({{5, 6, -1}, {5, 6, -1}, {1, -2, 5}}, A.get_vec());
  }

  SUBCASE("4x4 invertible matrix") {
    SquareMatrix A(
        {{1, 2, 3, 4}, {1, 4, 2, 0}, {1, -2, 5, -10}, {1, -2, 7, -10}});
    A.leb_inv();
    // A.calc_cout();
    compare_vectors(
        {{double(10) / 11, double(-3) / 11, double(26) / 11, -2},
         {double(-5) / 22, double(7) / 22, double(-15) / 44, double(1) / 4},
         {0, 0, double(-1) / 2, double(1) / 2},
         {double(3) / 22, double(-1) / 11, double(-1) / 22, 0}},
        A.get_vec());
  }
}