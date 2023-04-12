#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using std::endl, std::vector, std::string;

// A class for 2D square matrix.
class SquareMatrix {
 private:
  vector<vector<double>> myMatrix;
  bool isAugmented;          // is matrix an augmented matrix
  std::string calculations;  // stores step-by-step explanations of last
                             // operation performed on matrix

  // Checks if matrix is a 2D square matrix or valid augmented matrix.
  bool isValid(vector<vector<double>> initMatrix, bool _isAugmented);

  // Compares two floating point numbers and returns true if they are
  // approximately equal
  bool approxEqual(double a, double b);

  // Returns an (n x n) identity matrix
  vector<vector<double>> get_identity(int n);

  // Returns the row index of a row > startRow having a non-zero entry in a
  // specified column. If no such row found, return startRow.
  int get_next_pivot_row(int col);

  // validates a given index. If invalid, output an error message
  void validate(int i, string msg);

 public:
  SquareMatrix(vector<vector<double>> initMatrix, bool _isAugmented = 0);
  SquareMatrix();

  SquareMatrix operator+(SquareMatrix A);
  SquareMatrix operator-(SquareMatrix A);
  SquareMatrix operator-();
  SquareMatrix operator*(SquareMatrix A);

  // Rotates matrix 90deg clockwise.
  void transpose();

  // Returns square matrix as a 2D vector
  vector<vector<double>> get_vec();

  // Outputs most recent calculation performed on matrix to  console
  void calc_cout();

  // Merges two N x N square matrices and returns an augmented matrix [A | B]
  static vector<vector<double>> merge_matrices(vector<vector<double>> A,
                                               vector<vector<double>> B);

  // Inverse matrix using Leibniz method with cofactor formula
  vector<vector<double>> leb_inv();

  // Returns a stringified version of matrix
  std::string stringify(int dp = 3);

  // Returns calculations in the form of a string
  std::string get_calc();

  // Performs row1 + k * row2 and saves result to row1
  void add_rows(int row1, int row2, double k);

  // Scales a row of matrix by dividing each element by k.
  void scale_row(int row, double k);

  // Returns element at i-th row and j-th column where i, j are zero-based
  // indices.
  double at(int i, int j);

  // Sets element at i-th row and j-th column to x where i, j are zero-based
  void set_val(int i, int j, double x);

  // Swaps rows of matrix
  void swap_row(int row1, int row2);

  // Swaps columns of matrix
  void swap_col(int col1, int col2);

  // Returns rank of matrix.
  int rank();

  // Returns the trace. If matrix augmented, return trace of coefficient matrix
  double trace();

  // Returns determinant of matrix. If matrix augmented, return determinant
  // of coefficient matrix.
  double det();

  // Inverse matrix using Gauss-Jordan Elimination
  // method and returns the result.
  vector<vector<double>> gauss_inv();

  // Outputs approximations of a system of equations Ax = B using either
  // Gauss-Jacobi or Gauss-Seidel method
  vector<vector<double>> solve_approx(const bool useSeidelMethod,
                                      vector<double> initial_approx,
                                      int iterations = 5,
                                      int dp = 4);

  // Returns cofactor matrix
  SquareMatrix get_cofactor();

  // Returns minor of matrix
  double get_minor(int row, int col);

  // Get LU factorization. Does not modify original matrix.
  std::unordered_map<char, SquareMatrix> get_PLU();

  // Converts matrix to row echelon form using elementary row operations
  void to_ref();

  // Converts matrix to reduced row echelon form using row operations
  void to_rref();

  // Solves a system of equation Ax = B using cramers rule.
  vector<double> solve_cramer();

  // Solves a system of equation Ax = B using LU decomposition.
  vector<double> solve_plu();

  // Returns true if matrix is diagonally dominant.
  bool is_diag_dominant(bool strict = false);

  // Makes matrix diagonally dominant
  void to_diag(bool strict);
};
