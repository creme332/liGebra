# liGebra
A C++17 linear algebra library built for educational purposes. It uses the command line interface to output step-by-step calculations.

I created this to experiment with OOP and `doctest`.

# Features
- Basic matrix operations (addition, subtraction, multiplication, transpose).
- Inverse a matrix.
	- Gauss-Jordan elimination method.
	- Leibniz method.
- Determinant calculator.
	- Laplace expansion method.
	- Gauss-Jordan elimination method.
- Rank of matrix calculator.
- Convert matrix to row echelon form and reduced row echelon form.
- Make a matrix diagonally dominant through row/column exchanges.
- LU/PLU factorization (crouts method) with partial pivoting.
- Solving a system of linear equations.
	- Gauss-Jordan elimination method.
	- Cramers rule method.
	- LU/PLU decomposition method with partial pivoting.
- Iterative methods of solving a system of linear equations.
	- Gauss-Jacobi method.
	- Gauss-Seidel method.
- Eigenvectors and eigenvalues calculator. (*)

# Installation
To run this program, you will need to have a C++ compiler that supports C++17 installed on your machine.

Once you have a compiler installed, you can download the source code for this program from the GitHub repository.

# Usage

# To-do
add endl in ref row echelong
add tests lu factorization

gauss-jacoi and seidel resetting calculation string
review ddm algorithm
strassen algorithm

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