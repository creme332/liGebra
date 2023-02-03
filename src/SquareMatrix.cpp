#include "SquareMatrix.h"

SquareMatrix::SquareMatrix(vector<vector<double>> initMatrix,
                           bool _isAugmented) {
  calculations = "";
  if (isValid(initMatrix, _isAugmented)) {
    myMatrix = initMatrix;
    isAugmented = _isAugmented;
  } else {
    throw std::invalid_argument("Matrix is not a square matrix.");
  }
}

bool SquareMatrix::isValid(vector<vector<double>> initMatrix,
                           bool _isAugmented) {
  const int rowCount = initMatrix.size();
  if (rowCount == 0)
    return 0;

  const int colCount = initMatrix[0].size();
  if (colCount == 0)
    return 0;

  for (int i = 0; i < rowCount; i++) {
    if (_isAugmented) {
      if (initMatrix[i].size() - 1 != colCount - 1)
        return 0;
    } else if (initMatrix[i].size() != rowCount) {
      return 0;
    }
  }
  return 1;
}

double SquareMatrix::at(int row, int col) {
  if (row < 0 || col < 0 || row >= myMatrix.size() || col >= myMatrix[0].size())
    throw std::invalid_argument("Matrix must augmented.");

  return myMatrix[row][col];
}

double SquareMatrix::det() {
  if (myMatrix.size() == 1) {
    return myMatrix[0][0];
  }

  double determinant = 0;
  // loop through each element in row 0
  for (int i = 0; i < myMatrix[0].size(); i++) {
    // apply cofactor formula
    double el = i % 2 == 0 ? myMatrix[0][i] : -myMatrix[0][i];

    // get minor of current element
    determinant += SquareMatrix(myMatrix, isAugmented).get_minor(0, i) * el;
  }
  return determinant;
}

bool SquareMatrix::is_diagonally_dominant() {
  const int lastCol = isAugmented ? myMatrix[0].size() - 1 : myMatrix[0].size();
  for (int row = 0; row < myMatrix.size(); row++) {
    // get absolute sum of all elements on current row except diagonal element
    double sum = 0;
    for (int col = 0; col < lastCol; col++) {
      if (row != col) {
        sum += abs(myMatrix[row][col]);
      }
    }
    if (abs(myMatrix[row][row]) < sum)
      return 0;
  }
  return 1;
}

vector<vector<double>> SquareMatrix::solve_approx(
    const bool useSeidelMethod,
    const vector<double> initial_approx,
    const int iterations,
    const int dp) {
  std::stringstream stringRep;

  if (!isAugmented) {
    throw std::invalid_argument("Matrix must augmented.");
  }

  if (iterations < 1) {
    throw std::invalid_argument(
        "Number of iterations must be a positive value");
  }

  if (dp < 1) {
    throw std::invalid_argument(
        "Number of decimal places must be a positive value");
  }

  const int vars = myMatrix.size();

  if (initial_approx.size() != vars) {
    throw std::invalid_argument(
        "Size of approximation array must be equal to number of variables");
  }

  // calculate approximations
  vector<vector<double>> table(iterations + 1, initial_approx);
  for (int row = 1; row < iterations + 1; row++) {
    for (int col = 0; col < vars; col++) {
      table[row][col] = myMatrix[col][vars];

      // calculate value of table[row][col]
      for (int i = 0; i < vars; i++) {
        if (i == col)
          continue;
        if (useSeidelMethod) {
          // for seidel method, use the most recently calculated values to
          // calculate new value
          table[row][col] -=
              myMatrix[col][i] * (i < col ? table[row][i] : table[row - 1][i]);

        } else {
          // for jacobi method, use values from previous iteration to
          // calculate new value
          table[row][col] -= myMatrix[col][i] * table[row - 1][i];
        }
      }
      if (myMatrix[col][col] == 0) {
        throw std::invalid_argument(
            "Leading diagonal elements of matrix must be non-zero.");
      }
      table[row][col] /= myMatrix[col][col];
    }
  }

  // output table
  stringRep << (useSeidelMethod ? "Gauss-seidel method" : "Gauss-jacobi method")
            << endl;
  if (!is_diagonally_dominant()) {
    stringRep << ("Solutions may not converge as matrix is NOT diagonally "
                  "dominant.")
              << endl;
  }
  const char separator = ' ';
  const int nameWidth = 15;
  const int numWidth = 10;

  // output table header
  stringRep << std::left << std::setw(nameWidth) << std::setfill(separator)
            << "Iteration";
  for (int col = 0; col < vars; col++) {
    std::string a = "x" + std::to_string(col + 1);
    stringRep << std::left << std::setw(nameWidth) << std::setfill(separator)
              << a;
  }
  stringRep << endl;

  // output table values
  for (int row = 0; row < iterations + 1; row++) {
    stringRep << std::left << std::setw(nameWidth) << std::setfill(separator)
              << row;
    for (int col = 0; col < vars; col++) {
      stringRep << std::left << std::setw(nameWidth) << std::setfill(separator)
                << std::fixed << std::setprecision(dp) << table[row][col];
    }
    stringRep << endl;
  }

  calculations = stringRep.str();
  return table;
}

string SquareMatrix::stringify(int dp) {
  std::stringstream stringRep;
  const char separator = '|';
  const char spaceChar = ' ';
  const int nameWidth = 10;
  const int numWidth = 6;
  const int rowCount = myMatrix.size();

  for (int row = 0; row < rowCount; row++) {
    for (int col = 0; col < myMatrix[row].size(); col++) {
      stringRep << std::left << std::setw(nameWidth) << std::setfill(spaceChar)
                << std::fixed << std::setprecision(dp) << myMatrix[row][col];
      if (col == rowCount - 1) {
        stringRep << std::left << std::setw(5) << std::setfill(spaceChar)
                  << separator;
      }
    }
    stringRep << endl;
  }
  return stringRep.str();
}

void SquareMatrix::add_rows(int row1, int row2, double k) {
  const int dimension = myMatrix[myMatrix.size() - 1].size();

  for (int col = 0; col < dimension; col++) {
    myMatrix[row1][col] += myMatrix[row2][col] * k;
  }
}

void SquareMatrix::scale_row(int row, double k) {
  for (int i = 0; i < myMatrix[row].size(); i++) {
    if (!approxEqual(myMatrix[row][i], 0))
      myMatrix[row][i] /= k;
  }
}

bool SquareMatrix::approxEqual(double a, double b) {
  return abs(a - b) < 1e-9;
}

void SquareMatrix::swap_row(int row1, int row2) {
  if (row1 < 0 || row2 < 0 || row1 >= myMatrix.size() ||
      row2 >= myMatrix.size()) {
    throw std::invalid_argument("Invalid row indices");
  }
  std::swap(myMatrix[row1], myMatrix[row2]);
}

void SquareMatrix::gauss_inv() {
  if (isAugmented) {
    throw std::invalid_argument("Cannot inverse an augmented matrix");
  }
  const int rowCount = myMatrix.size();

  // string containing all the steps to be printed.
  std::stringstream stringRep;

  // create an augmented matrix = [A | I]
  SquareMatrix AugmentedMatrix(
      merge_matrices(myMatrix, get_identity(myMatrix.size())), true);
  stringRep << AugmentedMatrix.stringify();

  stringRep << endl;

  // Convert left matrix of augmented matrix to an upper triangular matrix where
  // each leading diagonal element is 0 or 1.
  stringRep << "Create an upper triangular matrix in LHS" << endl;
  for (int i = 0; i < rowCount; i++) {
    // Make AugmentedMatrix[i][i] a pivot if possible

    // perform row swapping if required
    if (approxEqual(AugmentedMatrix.at(i, i), 0) && i != rowCount - 1) {
      // get row index of row where i-th element is not a 0
      int newPivotRow = AugmentedMatrix.getNextPivotRow(i + 1, i);
      if (newPivotRow != i) {
        // swap rows
        AugmentedMatrix.swap_row(i, newPivotRow);
        stringRep << "Swap rows " << newPivotRow + 1 << " and " << i + 1
                  << endl;
        stringRep << AugmentedMatrix.stringify() << endl;
      }
    }

    // Scale current row to make pivot a 1
    if (!approxEqual(AugmentedMatrix.at(i, i), 1) &&
        !approxEqual(AugmentedMatrix.at(i, i), 0)) {
      stringRep << "Divide R" << i + 1 << " by " << AugmentedMatrix.at(i, i)
                << endl;
      AugmentedMatrix.scale_row(i, AugmentedMatrix.at(i, i));
      stringRep << AugmentedMatrix.stringify() << endl;
    }

    // make current column a pivot column
    for (int j = i + 1; j < rowCount; j++) {
      // output step
      if (approxEqual(AugmentedMatrix.at(j, i), 1)) {
        stringRep << "R" << j + 1 << "  - R" << i + 1 << endl;
      } else {
        stringRep << "R" << j + 1 << "  - R" << i + 1 << " * "
                  << AugmentedMatrix.at(j, i) << endl;
      }
      AugmentedMatrix.add_rows(j, i, -AugmentedMatrix.at(j, i));

      stringRep << AugmentedMatrix.stringify() << endl;
    }
  }

  // check if inverse matrix exists. For inverse matrix to exist, product of
  // leading diagonal must be 1.
  bool isInvertible = 1;
  for (int i = 0; i < rowCount; i++) {
    if (approxEqual(AugmentedMatrix.at(i, i), 0)) {
      isInvertible = 0;
      break;
    }
  }

  if (!isInvertible) {
    stringRep << "\nNo Inverse \n";
    stringRep << stringRep.str();
    return;
  }

  // convert upper triangular matrix in LHS to an identity matrix
  stringRep << "Convert upper triangular matrix in LHS to an identity matrix"
            << endl;
  for (int i = rowCount - 2; i >= 0; i--) {
    for (int j = i; j >= 0; j--) {
      if (approxEqual(AugmentedMatrix.at(j, i + 1), 0)) {
        // this IF statement is for optimisation only and is optional
        continue;
      }

      // output step
      if (approxEqual(AugmentedMatrix.at(j, i + 1), 1)) {
        stringRep << "R" << j + 1 << "  - R" << (i + 1) + 1 << endl;
      } else {
        stringRep << "R" << j + 1 << "  - R" << (i + 1) + 1 << " * "
                  << AugmentedMatrix.at(j, i + 1) << endl;
      }
      AugmentedMatrix.add_rows(j, i + 1, -AugmentedMatrix.at(j, i + 1));
      stringRep << AugmentedMatrix.stringify() << endl;
    }
  }

  // extract inverse matrix from augmented matrix
  vector<vector<double>> inverseMatrix(rowCount, vector<double>(rowCount, 0));
  for (int row = 0; row < rowCount; row++) {
    for (int col = rowCount; col < 2 * rowCount; col++) {
      inverseMatrix[row][col - rowCount] = AugmentedMatrix.at(row, col);
    }
  }
  calculations = stringRep.str();
}

// Merges two N x N square matrices and returns an augmented matrix [A | B]
vector<vector<double>> SquareMatrix::merge_matrices(vector<vector<double>> A,
                                                    vector<vector<double>> B) {
  if (A.size() != B.size()) {
    throw std::invalid_argument(
        "Cannot merge matrices having different number of rows");
  }
  if (A.size() == 0 || B.size() == 0) {
    throw std::invalid_argument("Cannot merge empty matrices");
  }
  const int N = A.size();
  vector<vector<double>> augmentedMatrix(N, vector<double>(N * 2, 0));

  // Merge matrices
  for (int row = 0; row < N; row++) {
    for (int col = 0; col < 2 * N; col++) {
      augmentedMatrix[row][col] = col < N ? A[row][col] : B[row][col - N];
    }
  }
  return augmentedMatrix;
}

vector<vector<double>> SquareMatrix::get_identity(int n) {
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

// Returns the row index of a row > startRow having a non-zero entry in a
// specified column. If no such row found, return startRow.
int SquareMatrix::getNextPivotRow(int startRow, int col) {
  for (int i = startRow; i < myMatrix.size(); i++) {
    if (!approxEqual(myMatrix[i][col], 0))
      return i;
  }
  return startRow;
}

void SquareMatrix::calc_cout() {
  std::cout << calculations;
}

void SquareMatrix::validate(int i, string msg) {
  if (i < 0 || i >= myMatrix.size()) {
    throw std::invalid_argument(msg);
  }
}

double SquareMatrix::get_minor(int row, int col) {
  validate(row, "Cannot calculate minor. Invalid value of row");
  validate(col, "Cannot calculate minor. Invalid value of col");
  vector<vector<double>> minor_matrix;
  for (int i = 0; i < myMatrix.size(); i++) {
    vector<double> b;
    for (int j = 0; j < myMatrix.size(); j++) {
      if (i != row && j != col) {
        b.push_back(myMatrix[i][j]);
      }
    }
    if (b.size() > 0)
      minor_matrix.push_back(b);
  }
  return SquareMatrix(minor_matrix).det();
}

vector<vector<double>> SquareMatrix::get_cofactor() {
  if (isAugmented) {
    throw std::invalid_argument("Cannot find cofactor of augmented matrix.");
  }
  SquareMatrix x(myMatrix, isAugmented);
  const int row_count = myMatrix.size();
  vector<vector<double>> minor_matrix(row_count, vector<double>(row_count, 0));

  for (int i = 0; i < row_count; i++) {
    for (int j = 0; j < row_count; j++) {
      minor_matrix[i][j] = ((i + j) % 2 ? 1 : -1) + x.get_minor(i, j);
    }
  }

  return minor_matrix;
}
