# liGebra

![Badge for test workflow](https://github.com/creme332/liGebra/actions/workflows/test.yml/badge.svg)[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/creme332/liGebra/495e7dc07fb51714e52ba4b2138332ca3e04c2a8?urlpath=lab%2Ftree%2Fdocumentation.ipynb)

A basic $C\texttt{++17}$ linear algebra library built for educational purposes. It uses the command line interface to output step-by-step calculations.

```cpp
SquareMatrix A({{5, 6, -1}, {1, 4, 2}, {1, -2, 5}});
A.to_ref(); // convert matrix to row echelon form
A.calc_cout(); // output calculations
```

[Try it online 🎉](https://mybinder.org/v2/gh/creme332/liGebra/495e7dc07fb51714e52ba4b2138332ca3e04c2a8?urlpath=lab%2Ftree%2Fdocumentation.ipynb)

# Features
* Basic matrix operations (addition, subtraction, multiplication, transpose).
* Inverse of a matrix.
  + Gauss-Jordan elimination method.
  + Leibniz method.
* Calculates matrix properties.
  + Determinant.
    - Laplace expansion method.
    - Gauss-Jordan elimination method.
  + Rank of matrix.
  + Trace of matrix.
* Convert matrix to row echelon form and reduced row echelon form.
* LU/PLU factorization (Crout's method) with partial pivoting.
* Solving system of linear equations.
  + Gauss-Jordan elimination method.
  + Cramers rule method.
  + LU/PLU decomposition method with partial pivoting.
  + Iterative methods
    - Gauss-Jacobi method.
    - Gauss-Seidel method.
* Pre-conditioning for iterative methods by making matrix strictly diagonally dominant.
* Chaining of operations is possible.
* Unit-testing with `doctest` library. 
* Automatic update to shared library using Github Actions.
# Installation

To use this library locally, you will need a compiler that supports $C\texttt{++17}$ installed on your machine. (MSVC or g++ compiler should work fine) 

Once you have a compiler installed, you can download the `src` folder from the GitHub repository. Import the `SquareMatrix.h` file in the `.cpp` where you want to use the library. For example, if your file structure is as follows:

```
src/
├─ SquareMatrix.h
├─ SquareMatrix.cpp
main.cpp
```

then your `main.cpp` file should have `#include "src/SquareMatrix.h"` at the top.

> ⚠ Do not include `test_runner.cpp` and the `test` folder in your project. 

## Shared library

If you want to use a shared library instead, you can generate the shared library with:

```bash
g++ -o libligebra.so -fpic -shared src/SquareMatrix.cpp -std=c++17
```

Then follow this [tutorial](https://betterprogramming.pub/how-to-build-a-linux-shared-library-f5b574b0c08e) on how to integrate the library in your code.

If you are using VS Code, assuming you have my `.vscode` folder, simply navigate to `test_runner.cpp` and pres `F5` .

## Testing

All files required for testing are found in the `tests` folder. Tests are using `doctest` library.
To run the tests (assuming you have a g++ compiler):

```bash
g++ -std=c++17 test_runner.cpp tests/tests.cpp src/SquareMatrix.cpp -W
./a.out
```

> ⚠ There should not be any other `.cpp` files with a `main()` function in your project. 

Comment the following lines in `test_runner.cpp` :

```cpp
  if (test_result == 1)
    throw std::runtime_error("Test failed");
```

The above lines causes the github workflow in `test.yml` to fail whenever a test fails. They are not necessary for testing locally.

If you are using VS Code, assuming you have my `.vscode` folder, simply navigate to `test_runner.cpp` and pres `F5` .

# To-do
* [x] Add binder link to documentation
* [x] Update task.json : add task to update so.
* [ ] [Reduce binder setup time](https://discourse.jupyter.org/t/how-to-reduce-mybinder-org-repository-startup-time/4956)
* [ ] Split tests into several files
* [ ] Calculate eigenvalues + vectors using QR algorithm
* [ ] Add tests for console output (Test string spit out by `calc_cout`)
* [ ] Add bodmas support for matrix operations
* [ ] Add option to output in improper fraction form instead of decimals.
# Limitations
* No support for complex numbers.
* Can only convert matrices to non-strict diagonal dominant form.
* Data type of numbers is limited to `double`.
* Supports only square matrices for coefficient matrix.
* Matrix operations do not follow BODMAS.
