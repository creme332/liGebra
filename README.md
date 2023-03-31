# liGebra
![Badge for test workflow](https://github.com/creme332/liGebra/actions/workflows/test.yml/badge.svg)

A basic C++17 linear algebra library built for educational purposes. It uses the command line interface to output step-by-step calculations.

I created this to experiment with OOP and `doctest`.

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
	- Eigenvectors and eigenvalues.
- Convert matrix to row echelon form and reduced row echelon form.
- Make a matrix diagonally dominant using elementary row operations.
- LU/PLU factorization (crouts method) with partial pivoting.
- Solving system of linear equations.
	- Gauss-Jordan elimination method.
	- Cramers rule method.
	- LU/PLU decomposition method with partial pivoting.
	- Iterative methods
		- Gauss-Jacobi method.
		- Gauss-Seidel method.

# Installation
To run this program, you will need to have a C++ compiler that supports C++17 installed on your machine.

Once you have a compiler installed, you can download the source code for this program from the GitHub repository.

# Usage

### `solve_cramer`

# To-do
add endl in ref row echelong
gauss-jacoi and seidel resetting calculation string
review ddm algorithm

- [ ] Add examples folder 
- [ ] Add bodmas support for matrix operations
- [ ] Add unary minus as operator
- [ ] Calculate determinant using row operations
- [ ] Allow chaining of functions in output. First transpose then find inverse in a single step.
- [ ] Add support for complex numbers
- [ ] Add GUI with qmake
- [ ] Add feature to output in decimals improper fraction form

Similar project: https://github.com/JNygard/MatrixCalculator

# Usage
Navigate to the `main.cpp` file and write your code.
Run code and see the console output.

## Testing
Run `tests.cpp` file in the `tests` folder.

# Limitations
- No support for complex numbers.
- Data type of numbers is limited to `double`.
- Supports only square matrices.
- Matrix operations do not follow BODMAS.