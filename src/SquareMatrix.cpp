#include "SquareMatrix.h"

SquareMatrix::SquareMatrix(vector<vector<double>> initMatrix,
                           bool _isAugmented) {
  calculations = "";
  if (isValid(initMatrix, _isAugmented)) {
    myMatrix = initMatrix;
    isAugmented = _isAugmented;
  } else {
    throw std::invalid_argument("Invalid matrix or augmented matrix.");
  }
}
SquareMatrix::SquareMatrix() {
  calculations = "";
  myMatrix = {{}};
  isAugmented = false;
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
      if (int(initMatrix[i].size()) - 1 != colCount - 1)
        return 0;
    } else if (int(initMatrix[i].size()) != rowCount) {
      return 0;
    }
  }
  return 1;
}

double SquareMatrix::at(int row, int col) {
  if (row < 0 || col < 0 || row >= int(myMatrix.size()) ||
      col >= int(myMatrix[0].size()))
    throw std::invalid_argument("Matrix must augmented.");

  return myMatrix[row][col];
}

double SquareMatrix::det() {
  if (myMatrix.size() == 1) {
    return myMatrix[0][0];
  }

  double determinant = 0;
  const int rowCount = myMatrix.size();

  // loop through each element in first row
  for (int i = 0; i < rowCount; i++) {
    // apply cofactor formula
    double el = i % 2 == 0 ? myMatrix[0][i] : -myMatrix[0][i];

    // get minor of current element
    determinant += get_minor(0, i) * el;
  }
  return determinant;
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

  if (int(initial_approx.size()) != vars) {
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
  if (!is_diag_dominant()) {
    stringRep << ("Solutions may not converge as matrix is NOT diagonally "
                  "dominant.")
              << endl;
  }
  const char separator = ' ';
  const int nameWidth = 15;

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

  calculations += stringRep.str();
  return table;
}

string SquareMatrix::stringify(int dp) {
  std::stringstream stringRep;
  //! Corner characters consume 2 spaces
  const string TopLeftCorner = "┌";
  const string TopRightCorner = "┐";
  const string BottomLeftCorner = "└";
  const string BottomRightCorner = "┘";
  const char spaceChar = ' ';
  const char pipe = '|';

  const int columnWidth = 10;  // width of column in which values are placed
  const int paddingWidth = 5;  // left and right padding
  const int rowCount = myMatrix.size();
  const int columnCount = int(myMatrix[0].size());

  std::stringstream topRow;
  std::stringstream bottomRow;

  if (isAugmented) {
    topRow << std::left << std::setw(2 * paddingWidth + 2)
           << std::setfill(spaceChar) << TopLeftCorner;
    bottomRow << std::left << std::setw(2 * paddingWidth + 2)
              << std::setfill(spaceChar) << BottomLeftCorner;
  } else {
    topRow << std::left << std::setw(paddingWidth + 2)
           << std::setfill(spaceChar) << TopLeftCorner;
    bottomRow << std::left << std::setw(paddingWidth + 2)
              << std::setfill(spaceChar) << BottomLeftCorner;
  }

  for (int row = 0; row < columnCount; row++) {
    topRow << std::left << std::setw(columnWidth) << std::setfill(spaceChar)
           << spaceChar;
    bottomRow << std::left << std::setw(columnWidth) << std::setfill(spaceChar)
              << spaceChar;
  }

  topRow << std::right << std::setw(paddingWidth - 2) << std::setfill(spaceChar)
         << TopRightCorner;
  bottomRow << std::right << std::setw(paddingWidth - 2)
            << std::setfill(spaceChar) << BottomRightCorner;

  stringRep << topRow.str() << endl;

  // print values in matrix
  for (int row = 0; row < rowCount; row++) {
    stringRep << std::left << std::setw(paddingWidth) << std::setfill(spaceChar)
              << pipe;
    for (int col = 0; col < columnCount; col++) {
      stringRep << std::left << std::setw(columnWidth)
                << std::setfill(spaceChar) << std::fixed
                << std::setprecision(dp) << myMatrix[row][col];
      if (col == rowCount - 1 && isAugmented) {
        stringRep << std::left << std::setw(paddingWidth)
                  << std::setfill(spaceChar) << pipe;
      }
    }
    stringRep << std::left << std::setw(paddingWidth + 2)
              << std::setfill(spaceChar) << pipe;
    stringRep << endl;
  }
  stringRep << bottomRow.str() << endl;

  return stringRep.str();
}

void SquareMatrix::add_rows(int row1, int row2, double k) {
  const int dimension = myMatrix[myMatrix.size() - 1].size();

  for (int col = 0; col < dimension; col++) {
    myMatrix[row1][col] += myMatrix[row2][col] * k;
  }
}

void SquareMatrix::scale_row(int row, double k) {
  if (approxEqual(k, 0))
    throw std::invalid_argument("Cannot divide by 0");
  for (int i = 0; i < int(myMatrix[row].size()); i++) {
    myMatrix[row][i] /= k;
  }
}

bool SquareMatrix::approxEqual(double a, double b) {
  const double EPSILON = 1e-9;
  return fabs(a - b) < EPSILON;
}

void SquareMatrix::swap_row(int row1, int row2) {
  if (row1 < 0 || row2 < 0 || row1 >= int(myMatrix.size()) ||
      row2 >= int(myMatrix.size())) {
    throw std::invalid_argument("Invalid row indices");
  }
  std::swap(myMatrix[row1], myMatrix[row2]);
}

vector<vector<double>> SquareMatrix::gauss_inv() {
  if (isAugmented) {
    throw std::invalid_argument("Cannot inverse an augmented matrix");
  }
  const int rowCount = myMatrix.size();

  // string containing all the steps to be printed.
  std::stringstream stringRep;

  // create an augmented matrix = [A | I]
  SquareMatrix AugmentedMatrix(
      merge_matrices(myMatrix, get_identity(myMatrix.size())), true);

  AugmentedMatrix.to_ref();

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
    calculations += AugmentedMatrix.get_calc();
    calculations += "No Inverse\n";
    return {{}};
  }

  AugmentedMatrix.to_rref();
  stringRep << AugmentedMatrix.get_calc();

  // extract inverse matrix from augmented matrix
  vector<vector<double>> inverseMatrix(rowCount, vector<double>(rowCount, 0));
  for (int row = 0; row < rowCount; row++) {
    for (int col = rowCount; col < 2 * rowCount; col++) {
      inverseMatrix[row][col - rowCount] = AugmentedMatrix.at(row, col);
    }
  }
  stringRep << "Inverse matrix:\n"
            << SquareMatrix(inverseMatrix).stringify() << endl;
  calculations += stringRep.str();
  myMatrix = inverseMatrix;
  return inverseMatrix;
}

vector<vector<double>> SquareMatrix::leb_inv() {
  if (isAugmented) {
    throw std::invalid_argument("Cannot inverse an augmented matrix");
  }
  std::stringstream stringRep;
  const double determinant = det();
  SquareMatrix cof = get_cofactor();

  stringRep << "Determinant of matrix: \n";
  stringRep << determinant;
  stringRep << "\n\n";

  stringRep << "Cofactor matrix: \n";
  stringRep << cof.stringify();
  stringRep << "\n";

  cof.transpose();
  stringRep << "Adjoint matrix: \n";
  stringRep << cof.stringify();
  stringRep << "\n";

  if (determinant == 0) {
    stringRep << "Matrix has no inverse\n";
    calculations += stringRep.str();
    return {{}};
  }

  stringRep << "Inverse matrix: \n";
  for (int i = 0; i < int(myMatrix.size()); i++)
    cof.scale_row(i, determinant);
  stringRep << cof.stringify();
  myMatrix = cof.get_vec();
  calculations += stringRep.str();

  return myMatrix;
}

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

int SquareMatrix::get_next_pivot_row(int col) {
  for (int i = col; i < int(myMatrix.size()); i++) {
    if (!approxEqual(myMatrix[i][col], 0))
      return i;
  }
  return col;
}

void SquareMatrix::calc_cout() {
  std::cout << calculations;
}

void SquareMatrix::validate(int i, string msg) {
  if (i < 0 || i >= int(myMatrix.size())) {
    throw std::invalid_argument(msg);
  }
}

double SquareMatrix::get_minor(int row, int col) {
  validate(row, "Cannot calculate minor. Invalid value of row");
  validate(col, "Cannot calculate minor. Invalid value of col");
  vector<vector<double>> minor_matrix;
  for (int i = 0; i < int(myMatrix.size()); i++) {
    vector<double> b;
    for (int j = 0; j < int(myMatrix.size()); j++) {
      if (i != row && j != col) {
        b.push_back(myMatrix[i][j]);
      }
    }
    if (b.size() > 0)
      minor_matrix.push_back(b);
  }
  return SquareMatrix(minor_matrix).det();
}

SquareMatrix SquareMatrix::get_cofactor() {
  if (isAugmented) {
    throw std::invalid_argument("Cannot find cofactor of augmented matrix.");
  }
  const int row_count = myMatrix.size();
  vector<vector<double>> minor_matrix(row_count, vector<double>(row_count, 0));

  for (int i = 0; i < row_count; i++) {
    for (int j = 0; j < row_count; j++) {
      minor_matrix[i][j] = ((i + j) & 1 ? -1 : 1) * get_minor(i, j);
    }
  }

  return SquareMatrix(minor_matrix);
}

void SquareMatrix::transpose() {
  for (int i = 0; i < int(myMatrix.size()); i++) {
    for (int j = 0; j < i; j++) {
      std::swap(myMatrix[i][j], myMatrix[j][i]);
    }
  }
}

vector<vector<double>> SquareMatrix::get_vec() {
  return myMatrix;
}

void SquareMatrix::to_ref() {
  /*
    A matrix is in echelon form if:
    1. Any rows consisting entirely of zeros are grouped at the bottom of the
    matrix.
    2. The first nonzero element of each row is 1.
    3. The leading 1 of each row after the first is positioned to the right of
    the leading 1 of the previous row.(This implies that all the elements below
    a leading 1 are zero.)
  */

  const int rowCount = myMatrix.size();
  bool REF = 1;  // is matrix already in REF?

  // string containing all the steps to be printed.
  std::stringstream stringRep;

  stringRep << "Converting matrix below to row echelon form" << endl;

  stringRep << stringify() << endl;

  for (int i = 0; i < rowCount; i++) {
    // Make A[i][i] a pivot if possible

    // perform row swapping if required
    if (approxEqual(myMatrix[i][i], 0) && i != rowCount - 1) {
      // get row index of row where i-th element is not a 0
      int newPivotRow = get_next_pivot_row(i);
      if (newPivotRow != i) {
        // swap rows
        swap_row(i, newPivotRow);
        stringRep << "Swap rows " << newPivotRow + 1 << " and " << i + 1
                  << endl;
        stringRep << stringify() << endl;
        REF = 0;
      }
    }

    // Scale current row to make pivot a 1
    if (!approxEqual(myMatrix[i][i], 1) && !approxEqual(myMatrix[i][i], 0)) {
      stringRep << "R" << i + 1 << " / " << myMatrix[i][i] << endl;
      scale_row(i, myMatrix[i][i]);
      stringRep << stringify() << endl;
      REF = 0;
    }

    // make current column a pivot column
    for (int j = i + 1; j < rowCount; j++) {
      // if element is already 0, move to next element
      if (approxEqual(myMatrix[j][i], 0)) {
        continue;
      }
      // output step
      stringRep << "R" << j + 1 << "  - R" << i + 1;
      if (!approxEqual(myMatrix[j][i], 1)) {
        stringRep << " * " << myMatrix[j][i];
      }

      add_rows(j, i, -myMatrix[j][i]);
      REF = 0;
      stringRep << endl << stringify() << endl;
    }
  }
  // if matrix was already in REF
  if (REF) {
    stringRep.str("");
    stringRep << "Matrix is in Row Echelon Form" << endl << stringify() << endl;
  }

  calculations += stringRep.str();
}

void SquareMatrix::to_rref() {
  /*
  A matrix is said to be in RREF if:
  1. matrix is in REF.
  2. The elements above and below a leading 1 are zero.
  */
  const int rowCount = myMatrix.size();

  // string containing all the steps to be printed.
  std::stringstream stringRep;

  bool RREF = 1;  // is matrix already in RREF?

  // convert matrix to row echelon form
  to_ref();

  for (int i = rowCount - 2; i >= 0; i--) {
    for (int row = i; row >= 0; row--) {
      if (approxEqual(myMatrix[row][i + 1], 0) ||
          approxEqual(myMatrix[i + 1][i + 1], 0)) {
        continue;
      }

      // output step
      stringRep << "R" << row + 1 << "  - R" << (i + 1) + 1;
      if (!approxEqual(myMatrix[row][i + 1], 1)) {
        stringRep << " * " << myMatrix[row][i + 1] << endl;
      }
      stringRep << endl;

      if (!approxEqual(myMatrix[row][i + 1], 0)) {
        add_rows(row, i + 1, -myMatrix[row][i + 1]);
        RREF = 0;
      }
      stringRep << stringify() << endl;
    }
  }
  // if matrix was already in RREF
  if (RREF) {
    stringRep.str("");
    stringRep << "Matrix is already in Reduced Row Echelon Form" << endl;
  }
  calculations += stringRep.str();
}

void SquareMatrix::swap_col(int col1, int col2) {
  const int lower_bound = 0;
  const int row_count = myMatrix.size();
  const int upper_bound = row_count + (isAugmented ? 0 : -1);

  // validate user input
  if (col1 < lower_bound || col2 < lower_bound || col1 > upper_bound ||
      col2 > upper_bound) {
    throw std::invalid_argument("Invalid column indices");
  }

  // swap columns
  for (int i = 0; i < row_count; i++) {
    std::swap(myMatrix[i][col1], myMatrix[i][col2]);
  }
}

vector<double> SquareMatrix::solve_cramer() {
  if (!isAugmented) {
    throw std::invalid_argument(
        "Matrix is must be an augmented matrix to use Cramer's rule.");
  }
  std::stringstream stringRep;
  const double determinant = det();
  stringRep << "Original matrix: \n\n" << stringify() << endl;

  stringRep << "Coefficient matrix: \n\n" << (get_coef()).stringify() << endl;

  stringRep << "Determinant = " << determinant << endl << endl;

  if (approxEqual(determinant, 0)) {
    throw std::invalid_argument(
        "Determinant must be non-zero to use Cramer's Rule");
  }

  const int rowCount = myMatrix.size();
  vector<double> solutions(rowCount, 0);

  for (int i = 0; i < rowCount; i++) {
    stringRep << "Swap column " << i + 1 << " and column " << rowCount + 1
              << " in original matrix" << endl
              << endl;
    swap_col(i, rowCount);
    stringRep << stringify() << endl;
    const double new_determinant = det();
    stringRep << "New determinant: " << new_determinant << endl;

    solutions[i] = new_determinant / determinant;
    stringRep << "x" << i + 1 << " = " << new_determinant << " / "
              << determinant << " = " << solutions[i] << endl
              << endl;
    swap_col(i, rowCount);
  }

  calculations += stringRep.str();

  return solutions;
}

bool SquareMatrix::is_diag_dominant(bool strict) {
  // Strict row diagonal dominance means that for each row, the absolute value
  // of the diagonal term is greater than the result of absolute values of other
  // terms

  // if matrix is augmented, check coefficient matrix only
  const int row_count = myMatrix.size();
  // if matrix is augmented, ignore columns after row_count.
  for (int i = 0; i < row_count; i++) {
    double row_sum = 0;
    for (int j = 0; j < row_count; j++) {
      row_sum += abs(myMatrix[i][j]);
    }
    if (strict) {
      // check for strict diagonal dominance
      if (abs(myMatrix[i][i]) <= row_sum - abs(myMatrix[i][i]))
        return 0;
    } else {
      if (abs(myMatrix[i][i]) < row_sum - abs(myMatrix[i][i]))
        return 0;
    }
  }
  return 1;
}

void SquareMatrix::to_diag() {
  // https://math.stackexchange.com/a/4697266/996485

  if (is_diag_dominant(true))
    return;
  if (approxEqual(get_coef().det(),
                  0))  // if coefficient matrix is singular, not possible.
    throw(
        "Singular matrix cannot be converted to strict diagonal dominance "
        "form");

  std::stringstream stringRep;
  const int row_count = myMatrix.size();

  // convert matrix to reduced row echlon form
  to_rref();

  stringRep << "Making matrix diagonally dominant" << endl;
  const double e = 0.5;

  // Starting from top row, add each row to the e*[row directly below it]
  for (int row = 0; row < row_count - 1; row++) {
    stringRep << "R" << row << " + " << e << " * R" << row + 1 << endl;
    add_rows(row, row + 1, e);
    stringRep << stringify() << endl;
  }

  // Starting from bottom row, add each row to the e*[row directly above it]
  for (int row = row_count - 1; row > 0; row--) {
    stringRep << "R" << row << " + " << e << " * R" << row - 1 << endl;
    add_rows(row, row - 1, e);
    stringRep << stringify() << endl;
  }

  calculations += stringRep.str();
}

SquareMatrix SquareMatrix::operator+(SquareMatrix otherMatrix) {
  if (isAugmented || otherMatrix.isAugmented) {
    throw std::invalid_argument("Cannot add augmented matrices");
  }
  // add matrices
  const int row_size = myMatrix.size();
  const int col_size = myMatrix[0].size();
  vector<vector<double>> result(row_size, vector<double>(col_size, 0));
  vector<vector<double>> A = otherMatrix.get_vec();

  // check if matrices have the same dimensions
  if (myMatrix.size() != A.size()) {
    throw std::invalid_argument("Matrices must have same dimensions");
  }

  for (int i = 0; i < row_size; i++) {
    for (int j = 0; j < col_size; j++) {
      result[i][j] = myMatrix[i][j] + A[i][j];
    }
  }
  return SquareMatrix(result, isAugmented);
}

SquareMatrix SquareMatrix::operator-(SquareMatrix otherMatrix) {
  if (isAugmented || otherMatrix.isAugmented) {
    throw std::invalid_argument("Cannot subtract augmented matrices");
  }
  const int row_size = myMatrix.size();
  const int col_size = myMatrix[0].size();
  vector<vector<double>> result(row_size, vector<double>(col_size, 0));
  vector<vector<double>> A = otherMatrix.get_vec();

  // check if matrices have the same dimensions
  if (myMatrix.size() != A.size()) {
    throw std::invalid_argument("Matrices must have same dimensions");
  }

  for (int i = 0; i < row_size; i++) {
    for (int j = 0; j < col_size; j++) {
      result[i][j] = myMatrix[i][j] - A[i][j];
    }
  }
  return SquareMatrix(result, isAugmented);
}

SquareMatrix SquareMatrix::operator-() {
  const int row_size = myMatrix.size();
  const int col_size = myMatrix[0].size();
  vector<vector<double>> result(row_size, vector<double>(col_size, 0));

  for (int i = 0; i < row_size; i++) {
    for (int j = 0; j < col_size; j++) {
      result[i][j] = -myMatrix[i][j];
    }
  }
  return SquareMatrix(result, isAugmented);
}

SquareMatrix SquareMatrix::operator*(SquareMatrix otherMatrix) {
  if (isAugmented || otherMatrix.isAugmented) {
    throw std::invalid_argument("Cannot multiply augmented matrices");
  }
  // multiply matrices
  const int row_size = myMatrix.size();
  vector<vector<double>> result(row_size, vector<double>(row_size, 0));
  vector<vector<double>> A = otherMatrix.get_vec();

  // check if matrices have the same dimensions
  if (myMatrix.size() != A.size()) {
    throw std::invalid_argument("Matrices must have same dimensions");
  }

  for (int i = 0; i < row_size; i++) {
    for (int j = 0; j < row_size; j++) {
      for (int k = 0; k < row_size; k++) {
        result[i][j] += myMatrix[i][k] * A[k][j];
      }
    }
  }
  return SquareMatrix(result, isAugmented);
}

std::string SquareMatrix::get_calc() {
  return calculations;
}

void SquareMatrix::set_val(int i, int j, double x) {
  if (i < 0 || j < 0 || i >= int(myMatrix.size())) {
    throw std::invalid_argument("Invalid row/column indices");
  }
  myMatrix[i][j] = x;
}

std::unordered_map<char, SquareMatrix> SquareMatrix::get_PLU() {
  if (isAugmented) {
    throw std::invalid_argument(
        "Cannot find LU factorization for augmented matrix");
  }

  const int rowCount = myMatrix.size();
  SquareMatrix P(get_identity(rowCount));
  SquareMatrix L(get_identity(rowCount));
  SquareMatrix U(myMatrix);

  // string containing all the steps to be printed.
  std::stringstream stringRep;

  stringRep << "Initialise matrix U:" << endl;
  stringRep << U.stringify() << endl;

  for (int i = 0; i < rowCount; i++) {
    // Make A[i][i] a pivot if possible

    // perform row swapping if required
    if (approxEqual(U.at(i, i), 0) && i != rowCount - 1) {
      // get row index of row where i-th element is not a 0
      int newPivotRow = U.get_next_pivot_row(i);
      if (newPivotRow != i) {
        // swap rows
        U.swap_row(i, newPivotRow);
        P.swap_row(i, newPivotRow);

        assert(i >= 0 && newPivotRow >= 0);

        // perform swapping in L
        for (int k = 0; k < std::min(i, newPivotRow); k++) {
          double const a = L.at(i, k);
          L.set_val(i, k, L.at(newPivotRow, k));
          L.set_val(newPivotRow, k, a);
        }

        // output
        stringRep << "Matrix U: Swap rows " << newPivotRow + 1 << " and "
                  << i + 1 << endl;
        stringRep << U.stringify() << endl;
        stringRep << "Update Matrix L" << endl;
        stringRep << L.stringify() << endl;
        stringRep << "Update Matrix P" << endl;
        stringRep << P.stringify() << endl;
      }
    }

    // Scale current row to make pivot a 1
    if (!approxEqual(U.at(i, i), 1) && !approxEqual(U.at(i, i), 0)) {
      const double x = U.at(i, i);

      stringRep << "Matrix U: R" << i + 1 << " / " << x << endl;
      L.set_val(i, i, x);
      U.scale_row(i, x);

      stringRep << U.stringify() << endl;
      stringRep << "Update Matrix L" << endl;
      stringRep << L.stringify() << endl;
    }

    // make current column a pivot column
    for (int j = i + 1; j < rowCount; j++) {
      // if element is already 0, move to next element
      if (approxEqual(U.at(j, i), 0)) {
        continue;
      }
      // output step
      stringRep << "Matrix U: R" << j + 1 << "  - R" << i + 1;
      if (!approxEqual(U.at(j, i), 1)) {
        stringRep << " * " << U.at(j, i);
      }
      stringRep << endl;

      L.set_val(j, i, U.at(j, i));
      U.add_rows(j, i, -U.at(j, i));

      stringRep << U.stringify() << endl;
      stringRep << "Update Matrix L" << endl;
      stringRep << L.stringify() << endl;
    }
  }

  stringRep << "Final Matrix P:" << endl;
  stringRep << P.stringify() << endl;

  stringRep << "Final Matrix L:" << endl;
  stringRep << L.stringify() << endl;

  stringRep << "Final Matrix U:" << endl;
  stringRep << U.stringify() << endl;

  calculations += stringRep.str();

  std::unordered_map<char, SquareMatrix> answer = {{
                                                       'p',
                                                       P,
                                                   },
                                                   {
                                                       'l',
                                                       L,
                                                   },
                                                   {
                                                       'u',
                                                       U,
                                                   }};
  return answer;
}

SquareMatrix SquareMatrix::get_coef() {
  const int rowCount = myMatrix.size();
  SquareMatrix A(vector<vector<double>>(rowCount, vector<double>(rowCount, 0)));

  for (int i = 0; i < rowCount; i++) {
    for (int j = 0; j < rowCount; j++) {
      A.set_val(i, j, myMatrix[i][j]);
    }
  }
  return A;
}

vector<double> SquareMatrix::solve_plu() {
  if (!isAugmented) {
    throw std::invalid_argument("Matrix must be in augmented form");
  }

  std::stringstream stringRep;
  const int rowCount = myMatrix.size();
  stringRep << "Solving AX = B" << endl;

  // extract coefficient matrix A from myMatrix
  SquareMatrix A = get_coef();
  stringRep << "Matrix A" << endl << A.stringify() << endl;

  // extract column vector B from myMatrix
  stringRep << "Matrix B = [";
  vector<double> B(rowCount, 0);
  for (int i = 0; i < rowCount; i++) {
    B[i] = myMatrix[i][rowCount];
    stringRep << B[i];
    if (i != rowCount - 1)
      stringRep << ", ";
  }
  stringRep << "]" << endl << endl;

  // get plu decomposition of A
  std::unordered_map<char, SquareMatrix> result = A.get_PLU();
  vector<vector<double>> P = result['p'].get_vec();
  vector<vector<double>> L = result['l'].get_vec();
  vector<vector<double>> U = result['u'].get_vec();
  stringRep << "Calculate PLU decomposition of A" << endl;
  stringRep << "Matrix P" << endl;
  stringRep << SquareMatrix(P).stringify() << endl;
  stringRep << "Matrix L" << endl;
  stringRep << SquareMatrix(L).stringify() << endl;
  stringRep << "Matrix U" << endl;
  stringRep << SquareMatrix(U).stringify() << endl;

  // calculate PB
  stringRep << "Matrix PB = [";
  vector<double> PB(rowCount, 0);

  for (int i = 0; i < rowCount; i++) {
    for (int j = 0; j < rowCount; j++) {
      PB[i] += P[i][j] * B[j];
    }
    stringRep << PB[i];
    if (i != rowCount - 1)
      stringRep << ", ";
  }
  stringRep << "]" << endl << endl;

  //  forward substitution : solve LZ = PB for Z
  stringRep << "Solve LZ = PB for Z" << endl;
  stringRep << "Matrix Z = [";
  vector<double> Z(rowCount, 0);

  for (int i = 0; i < rowCount; i++) {
    Z[i] = PB[i];
    for (int j = 0; j < i; j++) {
      Z[i] -= L[i][j] * Z[j];
    }
    Z[i] /= L[i][i];
    stringRep << Z[i];
    if (i != rowCount - 1)
      stringRep << ", ";
  }
  stringRep << "]" << endl << endl;

  //  backward substitution : solve UX = Z for X
  stringRep << "Solve UX = Z for X" << endl;
  stringRep << "Matrix X = [";
  vector<double> X(rowCount, 0);
  for (int i = rowCount - 1; i >= 0; i--) {
    X[i] = Z[i];
    for (int j = rowCount - 1; j > i; j--) {
      X[i] -= U[i][j] * X[j];
    }
    X[i] /= U[i][i];
    stringRep << X[i];
    if (i != 0)
      stringRep << ", ";
  }
  stringRep << "]" << endl;

  calculations += stringRep.str();
  return X;
}

double SquareMatrix::trace() {
  std::stringstream stringRep;
  stringRep << "Trace = ";
  double trace = 0;

  for (int i = 0; i < int(myMatrix.size()); i++) {
    stringRep << myMatrix[i][i];
    trace += myMatrix[i][i];
    if (i != int(myMatrix.size()) - 1)
      stringRep << " + ";
  }

  stringRep << " = " << trace << endl;
  calculations += stringRep.str();
  return trace;
}

int SquareMatrix::rank() {
  std::stringstream stringRep;
  SquareMatrix copy(myMatrix, isAugmented);
  copy.to_ref();
  stringRep << copy.get_calc();

  const int row_size = myMatrix.size();
  const int col_size = myMatrix[0].size();
  int rank = row_size;

  // count number of non-zero rows
  for (int i = 0; i < row_size; i++) {
    bool zero_row = 1;
    for (int j = 0; j < col_size; j++) {
      if (!approxEqual(copy.at(i, j), 0)) {
        zero_row = 0;
        break;
      }
    }
    if (zero_row)
      rank--;
  }
  stringRep << "Rank = number of non-zero rows = " << rank << endl;
  calculations += stringRep.str();

  return rank;
}