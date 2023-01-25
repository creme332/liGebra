#pragma once
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <string>
#include <vector>

using std::cout, std::cin, std::endl, std::vector;

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

  // void addRows(int i, int j);  // add
  // void swapRows(int i, int j);
  // void swapCols (int i, int j);
 public:
  SquareMatrix(vector<vector<double>> initMatrix, bool _isAugmented = 0);

  // Returns a stringified version of matrix
  std::string stringify();
  void setMatrix(vector<vector<double>> initMatrix);

  // Returns rank of matrix.
  double rank();

  // Returns determinant of matrix.
  double det();

  double getEigenValues();

  vector<vector<double>> getMinor(vector<vector<double>> A, int row, int col);

  // Returns coefficient matrix from an augmented matrix
  vector<vector<double>> getCoefficientMatrix(vector<vector<double>> augmented);

  bool is_diagonally_dominant();

  // Outputs approximations of a system of equations using either Gauss-Jacobi
  // or Gauss-Seidel method
  void solve_approx(const vector<double> initial_approx, const int iterations);

  // getLUF();  // get LU factorization double getMinor();
  // vector<vector<double>> getInverseGauss(bool printSteps);
  // vector<vector<double>> getInverseLeibniz(bool printSteps);  // normal
  // method

  // vector<vector<double>> getRREF();
};
