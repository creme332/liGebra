# liGebra
A single header library which carries out some basic linear algebra calculations. 
Step-by-step explanations are included.

I created this project to practice object-oriented programming in C++ and test-driven development with `doctest`.

# Features
- Basic matrix operations (addition, subtraction, multiplication, transpose).
- Strassen algorithm for matrix multiplication.
- Inverse a matrix.
	- Gauss-Jordan elimination method.
	- Leibniz method.
- Determinant calculator.
- Rank of matrix calculator.
- Convert matrix to reduced row echelon form.
- LU factorization.
- Solving system of linear equations.
	- Gauss elimination method
	- Cramers rule method.
	- LU decomposition method.
	- Gauss-Jacobi method.
	- Gauss-Seidel method.
- Eigenvectors and eigenvalues

# Code features
- OOP 
- No external dependencies
- C++14
- Test driven development

# To-do
- [ ] Use OOP
- [ ] create a class ligebra : ligebra.add(A, B, C)

```
-- src
-- tests
```
- [ ] Convert to Reduced row echelon form ([pseudocode on wikipedia](https://en.wikipedia.org/wiki/Row_echelon_form))
- [ ] prettyprint matrix to console : https://stackoverflow.com/questions/70991246/how-to-format-output-like-this
- [ ] Calculate determinant using row operations
- [ ] Use Catch library for testing like so: https://github.com/JNygard/MatrixCalculator
- [ ] Solving system of equations : Gauss, cramers, gauss jacobi, gauss seidel, eigenvectors
- [ ] Add GUI with qmake?

# Installation

- Download header from release.

> ** Warning **
>
> C++17 required.

# Usage

## Testing
Run tests

# Limitations
- Supports only numbers.
- Non-optimal for large matrices.
- Supports only square matrices.