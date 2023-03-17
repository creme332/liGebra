#include "../src//SquareMatrix.h"
#include "doctest.h"

void compare_1D_vector(vector<double> a,
                       vector<double> b,
                       const double tolerance = 0.001) {
  // Checks if two vectors have the same size and same elements
  CHECK_EQ(a.size(), b.size());

  for (int j = 0; j < a.size(); j++) {
    REQUIRE(a[j] - b[j] == doctest::Approx(0).epsilon(tolerance));
  }
}
void compare_2D_vectors(vector<vector<double>> a,
                        vector<vector<double>> b,
                        const double tolerance = 0.001) {
  // allow for a 0.1% error
  for (int i = 0; i < a.size(); i++) {
    compare_1D_vector(a[i], b[i]);
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

TEST_CASE("Test diagonal dominance") {
  SUBCASE("3x4 DD augmented matrix") {
    SquareMatrix A({{-22, 3, 8, 27}, {3, -6, 2, -14}, {3, 3, 21, 48}}, true);
    CHECK_EQ(A.is_diag_dominant(1), 1);
    CHECK_EQ(A.is_diag_dominant(), 1);
  }

  SUBCASE("3x4 non-DD augmented matrix") {
    SquareMatrix A({{-2, 3, 8, 27}, {3, -5, 2, -14}, {3, 3, 21, 48}}, true);
    CHECK_EQ(A.is_diag_dominant(1), 0);
    CHECK_EQ(A.is_diag_dominant(), 0);
  }

  SUBCASE("3x3 non-DD matrix") {
    SquareMatrix A({{-22, 3, 8}, {3, -10, 2}, {3, 3, 21}}, true);
    CHECK_EQ(A.is_diag_dominant(1), 1);
    CHECK_EQ(A.is_diag_dominant(), 1);
  }

  SUBCASE("3x3 non-strict DD matrix") {
    SquareMatrix A({{-22, 3, 8}, {3, -5, 2}, {3, 3, 21}}, true);
    CHECK_EQ(A.is_diag_dominant(1), 0);
    CHECK_EQ(A.is_diag_dominant(), 1);
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
      compare_2D_vectors(table, expectedAnswer);
    }
    SUBCASE("3x3 using Jacobi") {
      SquareMatrix A({{4, -1, -1, 3}, {-2, 6, 1, 9}, {-1, 1, 7, -6}}, true);

      vector<vector<double>> table = A.solve_approx(0, {0, 0, 0}, 3);
      const vector<vector<double>> expectedAnswer = {{0, 0, 0},
                                                     {0.75, 1.5, -0.857},
                                                     {0.911, 1.893, -0.964},
                                                     {0.982, 1.964, -0.997}};
      // A.calc_cout();
      compare_2D_vectors(table, expectedAnswer);
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
      compare_2D_vectors(expectedAnswer, table);
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
    compare_2D_vectors(expected, A.get_vec());
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
    compare_2D_vectors(expected, A.get_vec());
  }

  SUBCASE("3x3 non-invertible matrix") {
    SquareMatrix A({{5, 6, -1}, {5, 6, -1}, {1, -2, 5}});
    A.leb_inv();
    // A.calc_cout();
    compare_2D_vectors({{5, 6, -1}, {5, 6, -1}, {1, -2, 5}}, A.get_vec());
  }

  SUBCASE("4x4 invertible matrix") {
    SquareMatrix A(
        {{1, 2, 3, 4}, {1, 4, 2, 0}, {1, -2, 5, -10}, {1, -2, 7, -10}});
    A.leb_inv();
    // A.calc_cout();
    compare_2D_vectors(
        {{double(10) / 11, double(-3) / 11, double(26) / 11, -2},
         {double(-5) / 22, double(7) / 22, double(-15) / 44, double(1) / 4},
         {0, 0, double(-1) / 2, double(1) / 2},
         {double(3) / 22, double(-1) / 11, double(-1) / 22, 0}},
        A.get_vec());
  }
}

TEST_CASE("Test matrix inverse using Gauss") {
  SUBCASE("3x3 invertible matrix") {
    SquareMatrix A({{5, 6, -1}, {1, 4, 2}, {1, -2, 5}});
    vector<vector<double>> result = A.gauss_inv();
    // A.calc_cout();
    vector<vector<double>> expected = {
        {double(24) / 108, double(-28) / 108, double(16) / 108},
        {double(-3) / 108, double(26) / 108, double(-11) / 108},
        {double(-6) / 108, double(16) / 108, double(14) / 108}};
    compare_2D_vectors(expected, result);
    compare_2D_vectors(A.get_vec(), result);
  }

  SUBCASE("3x3 non-invertible matrix") {
    SquareMatrix A({{5, 6, -1}, {5, 6, -1}, {1, -2, 5}});
    vector<vector<double>> result = A.gauss_inv();
    // A.calc_cout();
    compare_2D_vectors({{}}, result);
  }

  SUBCASE("4x4 invertible matrix") {
    SquareMatrix A(
        {{1, 2, 3, 4}, {1, 4, 2, 0}, {1, -2, 5, -10}, {1, -2, 7, -10}});
    vector<vector<double>> result = A.gauss_inv();
    // A.calc_cout();
    compare_2D_vectors(
        {{double(10) / 11, double(-3) / 11, double(26) / 11, -2},
         {double(-5) / 22, double(7) / 22, double(-15) / 44, double(1) / 4},
         {0, 0, double(-1) / 2, double(1) / 2},
         {double(3) / 22, double(-1) / 11, double(-1) / 22, 0}},
        result);
    compare_2D_vectors(A.get_vec(), result);
  }
}

TEST_CASE("Test row echelon form") {
  SUBCASE("3x3 matrix") {
    SquareMatrix A({{5, 6, -1}, {1, 4, 2}, {1, -2, 5}});
    A.to_ref();
    // A.calc_cout();
    vector<vector<double>> expected = {
        {1, 1.2, -0.2}, {0, 1, 0.786}, {0, 0, 1}};
    compare_2D_vectors(expected, A.get_vec());
  }

  SUBCASE("3x3 matrix in REF with a zero row") {
    SquareMatrix A({{1, 5, 0}, {0, 0, 1}, {0, 0, 0}});
    A.to_ref();
    // A.calc_cout();
    compare_2D_vectors({{1, 5, 0}, {0, 0, 1}, {0, 0, 0}}, A.get_vec());
  }

  SUBCASE("3x3 matrix already in REF") {
    SquareMatrix A({{1, 5, 0}, {0, 0, 1}, {0, 0, 0}});
    A.to_ref();
    compare_2D_vectors({{1, 5, 0}, {0, 0, 1}, {0, 0, 0}}, A.get_vec());
  }

  SUBCASE("3-var [A | I] augmented matrix") {
    SquareMatrix A(
        {{3, -2, 4, 1, 0, 0}, {0, 2, -3, 0, 1, 0}, {4, 15, 11, 0, 0, 1}}, true);

    A.to_ref();
    compare_2D_vectors({{1, -2.0 / 3, 4.0 / 3, 1.0 / 3, 0, 0},
                        {0, 1, -3.0 / 2, 0, 0.5, 0},
                        {0, 0, 1, -8.0 / 193, -53.0 / 193, 6.0 / 193}},
                       A.get_vec());
  }
}

TEST_CASE("Test reduced row echelon form") {
  SUBCASE("3x3 matrix") {
    SquareMatrix A({{5, 6, -1}, {1, 4, 2}, {1, -2, 5}});
    A.to_rref();
    // A.calc_cout();
    vector<vector<double>> expected = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    compare_2D_vectors(expected, A.get_vec());
  }

  SUBCASE("3x3 matrix in REF with a zero row") {
    SquareMatrix A({{1, 5, 0}, {0, 0, 1}, {0, 0, 0}});
    A.to_rref();
    // A.calc_cout();
    compare_2D_vectors({{1, 5, 0}, {0, 0, 1}, {0, 0, 0}}, A.get_vec());
  }

  SUBCASE("4x4 system of equations with unique solution") {
    SquareMatrix A({{10, 1, 2, 3, 30},
                    {1, 15, 2, -5, 17},
                    {0, 1, 20, 3, 74},
                    {3, -10, -1, 25, 80}},
                   true);

    A.to_rref();
    // A.calc_cout();
    compare_2D_vectors(
        {{1, 0, 0, 0, 1}, {0, 1, 0, 0, 2}, {0, 0, 1, 0, 3}, {0, 0, 0, 1, 4}},
        A.get_vec());
  }

  SUBCASE("4x4 system of equations with infinite solutions") {
    SquareMatrix A({{10, 1, 2, 3, 30},
                    {10, 1, 2, 3, 30},
                    {10, 1, 2, 3, 30},
                    {10, 1, 2, 3, 30}},
                   true);

    A.to_rref();
    // A.calc_cout();
    compare_2D_vectors({{1, 0.1, 0.2, 0.3, 3},
                        {0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0}},
                       A.get_vec());
  }
}

TEST_CASE("Test to_diag()") {
  // TODO: update valid to make it work for non-augmented matrices

  // Returns true if matrix can be used with Gauss-seidel/Gauss-Jacobi
  auto valid = [](vector<vector<double>> augmentedMatrix) {
    for (int row = 0; row < augmentedMatrix.size(); row++) {
      // count number of non-zero elements in this row excluding augmented
      // column
      int c = 0;
      for (int col = 0; col < augmentedMatrix[row].size() - 1; col++) {
        if (augmentedMatrix[row][col] != 0) {
          c++;
        };
        if (c > 1)
          break;
      }
      if (c <= 1)
        return false;
    }
    return true;
  };
  SUBCASE("3-var system with circles on each row (Example 1)") {
    SquareMatrix A({{2, 6, -1, 85}, {6, 15, 2, 72}, {1, 1, 54, 110}}, true);
    A.to_diag(true);
    // A.calc_cout();
    CHECK_EQ(A.is_diag_dominant(), 1);
    compare_1D_vector(A.solve_cramer(), {-157.1661, 67.1726, 3.7036});
  }

  SUBCASE("3-var system with 2 initial circles") {
    // TODO: DOES NOT WORK AS EXPECTED
    SquareMatrix A({{1, 1, 1, 4}, {3, -10, 0, -17}, {2, -1, -1, -1}}, true);
    A.to_diag(false);
    // A.calc_cout();
    //   TODO: Check if each row has at least 2 non-zero elements

    CHECK_EQ(A.is_diag_dominant(), 1);
    compare_1D_vector(A.solve_cramer(), {1, 2, 1});
  }
}

TEST_CASE("Test Cramer's Rule") {
  SUBCASE("3-var system with infinite solutions") {
    SquareMatrix A({{2, 6, -1, 85}, {2, 6, -1, 15}, {1, 1, 54, 110}}, true);
    CHECK_THROWS_AS(A.solve_cramer(), std::exception);
  }

  SUBCASE("non-augmented matrix") {
    SquareMatrix A({{1, 1, 1, 1}, {0, 2, -6, 2}, {3, 6, -5, 4}, {0, 2, -6, 2}},
                   false);
    CHECK_THROWS_AS(A.solve_cramer(), std::exception);
  }

  SUBCASE("2-var system with unique solutions") {
    SquareMatrix A({{1, 1, 2}, {0, 2, 2}}, true);
    vector<double> solutions = A.solve_cramer();
    compare_1D_vector(solutions, {1, 1});
  }

  SUBCASE("3-var system with unique solutions") {
    SquareMatrix A({{1, 1, 1, 1}, {0, 2, -6, 2}, {3, 6, -5, 4}}, true);
    vector<double> solutions = A.solve_cramer();
    compare_1D_vector(solutions, {8, -5, -2});
  }

  SUBCASE("4-var system with unique solutions") {
    SquareMatrix A(
        {
            {1, 1, 1, 1, 2},
            {0, 2, -6, 2, 4},
            {3, 6, -5, 4, 10},
            {5, 1, -5, 4, 5},
        },
        true);
    vector<double> solutions = A.solve_cramer();
    compare_1D_vector(solutions, {0, 1, 0, 1});
    // A.calc_cout();
  }
}

TEST_CASE("Test operator") {
  SUBCASE("add 2x2") {
    SquareMatrix A({{2, 10}, {0, 0}});
    SquareMatrix B({{0, -1}, {1, 0}});
    compare_2D_vectors((A + B).get_vec(), {{2, 9}, {1, 0}});
  }

  SUBCASE("matrix operations on 2x2 and 3x3") {
    SquareMatrix A({{2, 10}, {0, 0}});
    SquareMatrix B({
        {0, -1, 5},
        {0, -1, 5},
        {0, -1, 5},
    });
    CHECK_THROWS_AS((A + B), std::exception);
    CHECK_THROWS_AS((A - B), std::exception);
    CHECK_THROWS_AS((A * B), std::exception);
  }

  SUBCASE("matrix operation on augmented matrices") {
    SquareMatrix A(
        {
            {0, -1, 5},
            {5, -1, 5},
        },
        true);
    CHECK_THROWS_AS((A * A), std::exception);
    CHECK_THROWS_AS((A - A), std::exception);
    CHECK_THROWS_AS((A + A), std::exception);
  }

  SUBCASE("subtract 2x2") {
    SquareMatrix A({{2, 10}, {0, 0}});
    SquareMatrix B({{0, -1}, {1, 0}});
    compare_2D_vectors((A - B).get_vec(), {{2, 11}, {-1, 0}});
  }
  SUBCASE("multiply 2x2") {
    SquareMatrix A({{2, 10}, {0, 0}});
    SquareMatrix B({{0, -1}, {1, 0}});
    compare_2D_vectors((A * B).get_vec(), {{10, -2}, {0, 0}});
  }
}