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

vector<vector<double>> SquareMatrix::getMinor(vector<vector<double>> A,
                                              int row,
                                              int col) {
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

double SquareMatrix::getDeterminant(vector<vector<double>> A) {
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

double SquareMatrix::det() {
  return getDeterminant(myMatrix);
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
  cout << (useSeidelMethod ? "Gauss-seidel method" : "Gauss-jacobi method")
       << endl;
  if (!is_diagonally_dominant()) {
    cout << ("Solutions may not converge as matrix is NOT diagonally "
             "dominant.")
         << endl;
  }
  const char separator = ' ';
  const int nameWidth = 15;
  const int numWidth = 10;

  // output table header
  cout << std::left << std::setw(nameWidth) << std::setfill(separator)
       << "Iteration";
  for (int col = 0; col < vars; col++) {
    std::string a = "x" + std::to_string(col + 1);
    cout << std::left << std::setw(nameWidth) << std::setfill(separator) << a;
  }
  cout << endl;

  // output table values
  for (int row = 0; row < iterations + 1; row++) {
    cout << std::left << std::setw(nameWidth) << std::setfill(separator) << row;
    for (int col = 0; col < vars; col++) {
      cout << std::left << std::setw(nameWidth) << std::setfill(separator)
           << std::fixed << std::setprecision(dp) << table[row][col];
    }
    cout << endl;
  }

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
