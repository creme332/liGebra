#include "../src//SquareMatrix.h"
#include "doctest.h"

TEST_CASE("Test constructor") {
  SUBCASE("2D matrix") {
    CHECK_THROWS_AS(SquareMatrix A({{1, 1}}), std::exception);
    CHECK_THROWS_AS(SquareMatrix A({{1, 1}, {1}}), std::exception);
    CHECK_THROWS_AS(SquareMatrix A({}), std::exception);
    CHECK_NOTHROW(SquareMatrix A({{1, 2}, {3, 4}}));
    CHECK_NOTHROW(SquareMatrix A({{1}}));
  }
  SUBCASE("Augmented matrix") {
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
    SquareMatrix A({{-22, 3, 8, 27}, {3, -5, 2, -14}, {3, 3, 21, 48}}, true);
    CHECK_EQ(A.is_diagonally_dominant(), true);
  }

    SUBCASE("3x4 non-DD augmented matrix") {
    SquareMatrix A({{-2, 3, 8, 27}, {3, -5, 2, -14}, {3, 3, 21, 48}}, true);
    CHECK_EQ(A.is_diagonally_dominant(), 0);
  }
}

//TEST_CASE("Test Gauss-Jacobi") {
//  SquareMatrix A({{-2, 3, 8, 27}, {3, -5, 2, -14}, {3, 3, 21, 48}}, true);
//  const vector<double> initial_approx = {0, 0, 0};
//  A.solve_approx(initial_approx, 4);
//}