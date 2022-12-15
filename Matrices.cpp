#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

// Prints a matrix or an augmented matrix
void printVector(vector<vector<double>> v, bool augmented = false) {
  const int dp = 10;
  if (augmented) {
    for (int row = 0; row < v.size(); row++) {
      for (int col = 0; col < 2 * v.size(); col++) {
        if (col == v.size() - 1) {
          cout << v[row][col] << " | ";
        } else {
          cout << v[row][col] << " ";
        }
      }
      cout << endl;
    }
  } else {
    cout << "{ ";
    for (auto a : v) {
      cout << "{";
      for (auto b : a) {
        cout << setprecision(dp) << b << ", ";
      }
      cout << "}, ";
    }
    cout << "}" << endl;
  }
  cout << endl;
}

// Returns the minor of the element at a given coordinates where A is a square
// matrix
vector<vector<double>> getMinor(vector<vector<double>> A, int row, int col) {
  vector<vector<double>> minor;
  for (int i = 0; i < A.size(); i++) {
    vector<double> b;
    for (int j = 0; j < A.size(); j++) {
      if (i != row && j != col) {
        b.push_back(A[i][j]);
      }
    }
    if (b.size() > 0)
      minor.push_back(b);
  }
  return minor;
}

// Returns the determinant of a square matrix using Laplace expansion
double getDeterminant(vector<vector<double>> A) {
  if (A.size() == 1) {
    return A[0][0];
  }

  int determinant = 0;
  // loop through each element in row 0
  for (int i = 0; i < A[0].size(); i++) {
    // apply cofactor formula
    double el = i % 2 == 0 ? A[0][i] : -A[0][i];

    // get minor of current element
    vector<vector<double>> minor = getMinor(A, 0, i);
    determinant += getDeterminant(minor) * el;
  }
  return determinant;
}

// Returns an (n x n) identity matrix
vector<vector<double>> getIdentityMatrix(int n) {
  if (n <= 0) {
    throw std::invalid_argument("Received an invalid matrix size");
    return {{}};
  }
  vector<vector<double>> identity(n, vector<double>(n, 0));
  for (int i = 0; i < n; i++) {
    identity[i][i] = 1;
  }
  return identity;
}

// Compares two floating point numbers and returns true if they are
// approximately equal
bool approxEqual(double a, double b) {
  return abs(a - b) < 1e-9;
}
// Returns an augmented matrix A|I where I is an identity matrix and A is a
// square matrix.

// Returns true if two 2D matrices are equal
bool isEqualMatrices(vector<vector<double>> A, vector<vector<double>> B) {
  // compare row sizes
  if (A.size() != B.size())
    return false;

  // compare column sizes
  for (int i = 0; i < A.size(); i++) {
    if (A[i].size() != B[i].size())
      return false;
  }
  // if both matrices are empty
  if (A.size() == 1 && A[0].size() == 0)
    return true;

  // compare elements
  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < A.size(); j++) {
      // if (A[i][j] != B[i][j]) {
      //   return false;
      // }
      if (!approxEqual(A[i][j], B[i][j])) {
        return false;
      }
    }
  }

  return true;
}

vector<vector<double>> getAugmentedMatrix(vector<vector<double>> A) {
  const int dimension = A.size();
  vector<vector<double>> augmentedMatrix(dimension,
                                         vector<double>(dimension * 2, 0));

  vector<vector<double>> identity = getIdentityMatrix(dimension);

  // Merge A and identity matrix into a single matrix
  for (int row = 0; row < dimension; row++) {
    for (int col = 0; col < 2 * dimension; col++) {
      augmentedMatrix[row][col] =
          col < dimension ? A[row][col] : identity[row][col - dimension];
    }
  }
  return augmentedMatrix;
}

// Performs row1 + k*row2 and saves result to row1
vector<vector<double>> addMatrixRows(vector<vector<double>> A,
                                     int row1,
                                     int row2,
                                     double k) {
  const int dimension = A[A.size() - 1].size();

  for (int col = 0; col < dimension; col++) {
    A[row1][col] += A[row2][col] * k;
  }

  return A;
}

// Swaps two rows of matrix.
vector<vector<double>> swapMatrixRows(vector<vector<double>> A,
                                      int row1,
                                      int row2) {
  vector<double> copyRow1 = A[row1];
  A[row1] = A[row2];
  A[row2] = copyRow1;
  return A;
}

// Scales a row of matrix by dividing each element by k.
vector<vector<double>> scaleMatrixRow(vector<vector<double>> A,
                                      int row,
                                      double k) {
  for (int i = 0; i < A[0].size(); i++) {
    if (!approxEqual(A[row][i], 0))
      A[row][i] /= k;
  }
  return A;
}

// Returns row index of a row after startRow where A[i][col] != 0.
int getSpecialRow(vector<vector<double>> A, int startRow, int col) {
  for (int i = startRow; i < A.size(); i++) {
    if (!approxEqual(A[i][col], 0))
      return i;
  }
  return -1;
}
// Calculates the inverse of a square matrix using Gauss-Jordan Elimination
// method and returns the result.
vector<vector<double>> getInverseMatrix(vector<vector<double>> v) {
  double determinant = getDeterminant(v);
  cout << "Determinant: " << determinant << endl;
  // if determinant is 0, no inverse matrix
  if (approxEqual(0, determinant)) {
    return {{}};
  }
  vector<vector<double>> AugmentedMatrix = getAugmentedMatrix(v);
  printVector(AugmentedMatrix, true);

  // form an upper triangular matrix with leading diagonal 1 in LHS of augmented
  // matrix
  cout << "Create an upper triangular matrix in LHS" << endl;
  for (int i = 0; i < v.size(); i++) {
    if (approxEqual(AugmentedMatrix[i][i], 0)) {
      // swap rows
      // get row index of row where i-th element is not a 0
      int row = getSpecialRow(AugmentedMatrix, i + 1, i);

      // swap rows
      AugmentedMatrix = swapMatrixRows(AugmentedMatrix, i, row);

      cout << "Swap rows " << row << "and " << i << endl;
    }

    // output step
    if (!approxEqual(AugmentedMatrix[i][i], 1)) {
      cout << "Divide R" << i << " by " << AugmentedMatrix[i][i] << endl;
      AugmentedMatrix =
          scaleMatrixRow(AugmentedMatrix, i, AugmentedMatrix[i][i]);
      printVector(AugmentedMatrix, true);
    }

    for (int j = i + 1; j < v.size(); j++) {
      // output step
      if (approxEqual(AugmentedMatrix[j][i], 1)) {
        cout << "R" << j << "  - R" << i << endl;
      } else {
        cout << "R" << j << "  - R" << i << " * " << AugmentedMatrix[j][i]
             << endl;
      }
      AugmentedMatrix =
          addMatrixRows(AugmentedMatrix, j, i, -AugmentedMatrix[j][i]);

      printVector(AugmentedMatrix, true);
    }
  }

  // convert upper triangular matrix in LHS to an identity matrix
  cout << "Convert upper triangular matrix in LHS to an identity matrix"
       << endl;
  for (int i = v.size() - 2; i >= 0; i--) {
    for (int j = i; j >= 0; j--) {
      if (approxEqual(AugmentedMatrix[j][i + 1], 0)) {
        // this IF statement is optional
        continue;
      }

      // output step
      if (approxEqual(AugmentedMatrix[j][i + 1], 1)) {
        cout << "R" << j << "  - R" << i + 1 << endl;
      } else {
        cout << "R" << j << "  - R" << i + 1 << " * "
             << AugmentedMatrix[j][i + 1] << endl;
      }

      AugmentedMatrix =
          addMatrixRows(AugmentedMatrix, j, i + 1, -AugmentedMatrix[j][i + 1]);

      printVector(AugmentedMatrix, true);
    }
  }

  // extract inverse matrix from augmented matrix
  vector<vector<double>> inverseMatrix(v.size(), vector<double>(v.size(), 0));
  for (int row = 0; row < v.size(); row++) {
    for (int col = v.size(); col < AugmentedMatrix[0].size(); col++) {
      inverseMatrix[row][col - v.size()] = AugmentedMatrix[row][col];
    }
  }

  return inverseMatrix;
}

// Run tests for functions
void runTests() {
  vector<vector<double>> A, expected, received;

  // 4x4 matrix with inverse
  A = {{1, 1, 2, 5}, {1, 3, 2, 5}, {2, 1, 0, 8}, {5, 4, 2, 9}};
  expected = {{
                  0.08333333333,
                  -0.4166666667,
                  -0.1666666667,
                  0.3333333333,
              },
              {
                  -0.5,
                  0.5,
                  0,
                  0,
              },
              {
                  0.6041666667,
                  -0.1458333333,
                  -0.3333333333,
                  0.04166666667,
              },
              {
                  0.04166666667,
                  0.04166666667,
                  0.1666666667,
                  -0.08333333333,
              }};
  received = getInverseMatrix(A);
  printVector(received);
  cout << "Test 1 : "
       << ((isEqualMatrices(received, expected) == 1) ? "Passed" : "Failed")
       << endl
       << endl;

  // 4x4 matrix with no inverse
  A = {{1, 3, 2, 5}, {1, 3, 2, 5}, {2, 1, 0, 8}, {5, 4, 2, 9}};
  expected = {{}};
  received = getInverseMatrix(A);
  printVector(received);
  cout << "Test 2 : "
       << ((isEqualMatrices(received, expected) == 1) ? "Passed" : "Failed")
       << endl
       << endl;

  // 4x4 zero matrix with no inverse
  A = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  expected = {{}};
  received = getInverseMatrix(A);
  printVector(received);
  cout << "Test 3 : "
       << ((isEqualMatrices(received, expected) == 1) ? "Passed" : "Failed")
       << endl
       << endl;

  // 4x4 matrix with inverse, R11=0
  A = {{0, 1, 2, 5}, {1, -1, 2, 5}, {-10, 1, 0, 8}, {5, 4, 2, 9}};
  expected = {{
                  -0.2,
                  0.1111111111,
                  -0.04444444444,
                  0.08888888889,
              },
              {
                  0.4,
                  -0.4444444444,
                  -0.02222222222,
                  0.04444444444,
              },
              {
                  1.05,
                  -0.2638888889,
                  -0.1694444444,
                  -0.2861111111,
              },
              {
                  -0.3,
                  0.1944444444,
                  0.07222222222,
                  0.1055555556,
              }};
  received = getInverseMatrix(A);
  printVector(received);
  cout << "Test 4 : "
       << ((isEqualMatrices(received, expected) == 1) ? "Passed" : "Failed")
       << endl
       << endl;

  // 3x3 matrix with inverse
  A = {{0, 1, 2}, {5, 2, 5}, {-10, 1, -1}};
  expected = {
      {
          -1.4,
          0.6,
          0.2,
      },
      {
          -9,
          4,
          2,
      },
      {
          5,
          -2,
          -1,
      },
  };
  received = getInverseMatrix(A);
  printVector(received);
  cout << "Test 5 : "
       << ((isEqualMatrices(received, expected) == 1) ? "Passed" : "Failed")
       << endl
       << endl;

  // 2x2 matrix
  A = {{5, 1}, {-1, 7}};
  expected = {
      {
          0.1944444444,
          -0.02777777778,
      },
      {
          0.02777777778,
          0.1388888889,
      },
  };
  received = getInverseMatrix(A);
  printVector(received);
  cout << "Test 6 : "
       << ((isEqualMatrices(received, expected) == 1) ? "Passed" : "Failed")
       << endl
       << endl;

  // 1x1 matrix
  A = {{5}};
  expected = {{0.2}};
  received = getInverseMatrix(A);
  printVector(received);
  cout << "Test 7 : "
       << ((isEqualMatrices(received, expected) == 1) ? "Passed" : "Failed")
       << endl
       << endl;
}
int main() {
  runTests();
  return 0;
}
