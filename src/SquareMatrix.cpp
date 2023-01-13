#include "SquareMatrix.h"
#include <stdexcept>

SquareMatrix::SquareMatrix(vector<vector<double>> initMatrix,
                           bool _isAugmented) {
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