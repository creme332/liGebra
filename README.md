# liGebra
![Badge for test workflow](https://github.com/creme332/liGebra/actions/workflows/test.yml/badge.svg)

A basic C++17 linear algebra calculator built for educational purposes. It uses the command line interface to output step-by-step calculations.
```cpp
SquareMatrix A({{5, 6, -1}, {1, 4, 2}, {1, -2, 5}});
A.to_ref(); // convert matrix to row echelon form
A.calc_cout(); // output calculations
```

[View documentation](documentation.ipynb)

# Features
- Basic matrix operations (addition, subtraction, multiplication, transpose).
- Inverse of a matrix.
	- Gauss-Jordan elimination method.
	- Leibniz method.
- Calculates matrix properties.
	- Determinant.
		- Laplace expansion method.
		- Gauss-Jordan elimination method.
	- Rank of matrix.
	- Trace of matrix.
- Convert matrix to row echelon form and reduced row echelon form.
- Make a matrix diagonally dominant using elementary row operations.
- LU/PLU factorization (Crout's method) with partial pivoting.
- Solving system of linear equations.
	- Gauss-Jordan elimination method.
	- Cramers rule method.
	- LU/PLU decomposition method with partial pivoting.
	- Iterative methods
		- Gauss-Jacobi method.
		- Gauss-Seidel method.
- Chaining of operations is possible.
- Tested with `doctest` library. 
- Github workflow to automatically update shared library.

# Installation
To run this program, you will need a compiler that supports $C\texttt{++17}$ installed on your machine. (MSVC or g++ compiler should work fine) 

Once you have a compiler installed, you can download the `src` folder from the GitHub repository. Add `#include "src/SquareMatrix.h"` to your  where you want to use the library. 

> ⚠ Do not include `test_runner.cpp` and the `test` folder in your project. 

## Shared library
If you want to use a shared library instead, you can generate the shared library with:
```bash
g++ -o libligebra.so -fpic -shared src/SquareMatrix.cpp -std=c++17
```
Then follow this [tutorial](https://betterprogramming.pub/how-to-build-a-linux-shared-library-f5b574b0c08e) on how to integrate the library in your code.

## Testing
All files required for testing are found in the `tests` folder. Tests are using `doctest` library.
To run your code as together with the tests (assuming you have a g++ compiler):
```linux
g++ test_runner.cpp tests/tests.cpp src/SquareMatrix.cpp -W
```
> ⚠ There should not be any other `.cpp` files with a `main()` function in your project. 

# To-do
- [ ] Review ddm algorithm
- [ ] Add binder link to documentation
- [ ] Convert into a header-only library instead of shared library.
- [ ] Update task.json : add task to update so.
- [ ] Split tests into several files
- [ ] Calculate eigenvalues + vectors using QR algorithm
- [ ] Add tests for console output (Test string spit out by `calc_cout`)
- [ ] Add bodmas support for matrix operations
- [ ] Add support for complex numbers
- [ ] Add option to output in improper fraction form instead of decimals.

Similar project: https://github.com/JNygard/MatrixCalculator

# Limitations
- No support for complex numbers.
- Data type of numbers is limited to `double`.
- Supports only square matrices for coefficient matrix.
- Matrix operations do not follow BODMAS.