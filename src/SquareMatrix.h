#pragma once
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using std::cout, std::cin, std::endl, std::vector, std::string;

// A class for 2D square matrix.
class SquareMatrix {
 private:
  vector<vector<double>> myMatrix;
  bool isAugmented;          // is matrix an augmented matrix
  std::string calculations;  // stores step-by-step explanations on last
                             // operation performed on matrix

  // Checks if matrix is a 2D square matrix or valid augmented matrix.
  bool isValid(vector<vector<double>> initMatrix, bool _isAugmented);

  // Calculates determinant using Laplace expansion.
  double getDeterminant(vector<vector<double>> A);

  // Compares two floating point numbers and returns true if they are
  // approximately equal
  bool approxEqual(double a, double b);

  // Returns an (n x n) identity matrix
  vector<vector<double>> get_identity(int n);

  // Returns the i-th row of matrix
  vector<double> get_row(int i);

  // Returns the row index of a row > startRow having a non-zero entry in a
  // specified column. If no such row found, return startRow.
  int getNextPivotRow(int startRow, int col);

 public:
  SquareMatrix(vector<vector<double>> initMatrix, bool _isAugmented = 0);

  static vector<vector<double>> merge_matrices(vector<vector<double>> A,
                                               vector<vector<double>> B);

  // Returns a stringified version of matrix
  std::string stringify(int dp = 3);

  // Performs row1 + k*row2 and saves result to row1
  void add_rows(int row1, int row2, double k);

  // Scales a row of matrix by dividing each element by k.
  void scale_row(int row, double k);

  // returns element at i-th row and j-th column where i, j are zero-based
  // indices.
  double at(int i, int j);

  // Swaps rows of matrix
  void swap_row(int row1, int row2);

  void setMatrix(vector<vector<double>> initMatrix);

  // Returns rank of matrix.
  double rank();

  // Returns determinant of matrix.
  double det();

  double getEigenValues();

  // Inverse matrix using Gauss-Jordan Elimination
  // method and returns the result.
  void gauss_inv();

  vector<vector<double>> getMinor(vector<vector<double>> A, int row, int col);

  // Returns coefficient matrix from an augmented matrix
  vector<vector<double>> getCoefficientMatrix(vector<vector<double>> augmented);

  bool is_diagonally_dominant();

  // Outputs approximations of a system of equations using either Gauss-Jacobi
  // or Gauss-Seidel method
  vector<vector<double>> solve_approx(const bool useSeidelMethod,
                                      vector<double> initial_approx,
                                      int iterations = 5,
                                      int dp = 4);

  // getLUF();  // get LU factorization double getMinor();
  // vector<vector<double>> getInverseGauss(bool printSteps);
  // vector<vector<double>> getInverseLeibniz(bool printSteps);  // normal
  // method

  // vector<vector<double>> getRREF();
};
