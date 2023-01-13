#pragma once
#include <iostream>
#include <vector>
using std::vector;

// A class for 2D square matrix.
class SquareMatrix {
 private:
  vector<vector<double>> myMatrix;
  bool isAugmented;  // is matrix an augmented matrix

  // Checks if matrix is a 2D square matrix or valid augmented matrix.
  bool isValid(vector<vector<double>> initMatrix, bool _isAugmented);

  // void addRows(int i, int j);  // add
  // void swapRows(int i, int j);

 public:
  SquareMatrix(vector<vector<double>> initMatrix, bool _isAugmented = 0);

  /// <summary>MyMethod is a method in the MyClass class.
  /// <para>Here's how you could make a second paragraph in a description. <see
  /// cref="System::Console::WriteLine"/> for information about output
  /// statements.</para> <seealso cref="MyClass::MyMethod2"/>
  /// </summary>
  // void stringify();  // converts matrix into a string ready to printed to
  // console void setMatrix(vector<vector<double>> initMatrix); double
  // getDeterminant(); double getRank(); double getEigenValues(); double
  // getLUF();  // get LU factorization double getMinor();
  // vector<vector<double>> getInverseGauss(bool printSteps);
  // vector<vector<double>> getInverseLeibniz(bool printSteps);  // normal
  // method

  // vector<vector<double>> getRREF();
};
