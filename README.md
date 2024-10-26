# <span style="font-size: 18px; font-weight: bold;">Console Application Development Using Numerical Method</span>

## <span style="font-size: 18px; font-weight: bold;">Table of Contents</span>
- [Introduction](#introduction)  
- [Algorithms](#algorithms)  
  - [1. Solutions for Linear Equations](#1-solutions-for-linear-equations)  
    - [a. Jacobi Iterative Method](#a-jacobi-iterative-method)  
    - [b. Gauss-Seidel Iterative Method](#b-gauss-seidel-iterative-method)  
    - [c. Gauss Elimination](#c-gauss-elimination)  
    - [d. Gauss-Jordan Elimination](#d-gauss-jordan-elimination)  
    - [e. LU Factorization](#e-lu-factorization)  
  - [2. Solutions for Non-linear Equations](#2-solutions-for-non-linear-equations)  
    - [a. Bisection Method](#a-bisection-method)  
    - [b. False Position Method](#b-false-position-method)
    - [c. Newton-Raphson Method](#c-newton-raphson-method) 
    - [d. Secant Method](#d-secant-method)  
     
  - [3. Solutions for Differential Equations](#3-solutions-for-differential-equations)  
    - [a. Runge-Kutta Method](#a-runge-kutta-method)  
  - [4. Matrix Inversion](#4-matrix-inversion)  
- [Conclusion](#conclusion)  

## <span style="font-size: 18px; font-weight: bold;">Introduction</span>  
The Numerical Methods Solver is an application designed to provide solutions for various mathematical problems, including linear equations, non-linear equations, differential equations, and matrix inversion. This README file outlines the algorithms used in the application, providing detailed explanations for each method, including their workings, applications, advantages, and disadvantages.  


# Algorithms

## 1. Solutions for Linear Equations

### a. Jacobi Iterative Method

The Jacobi method, named after German mathematician Carl Gustav Jacob Jacobi (1804 – 1851), is used for solving a system of linear equations:


a_11 * x_1 + a_12 * x_2 + ... + a_1n * x_n = b_1 

a_21 * x_1 + a_22 * x_2 + ... + a_2n * x_n = b_2 

a_n1 * x_1 + a_n2 * x_2 + ... + (a_nn) * x_n=b_n

To find the solution to this system Ax = B, we assume that the system has a unique solution and that there are no zero entries among the diagonal or pivot elements of the coefficient matrix A.

The equations can be rearranged as follows:


x_1 = (1/a_11)* (b_1 - a_12 * x_2 - a_13 * x_3 - ... - a_1n * x_n) 

x_2 = (1/a_22) *(b_2 - a_21 * x_1 - a_23 * x_3 - ...  - a_2n * x_n) 

x_n =(1/a_nn) * (b_n - a_n1 * x_1 - a_n2 *x_2 - ... - a_{n(n-1)} * x_{n-1})


By making an initial guess for the solution \(x^{(0)} = (x_1^{(0)}, x_2^{(0)}, \ldots, x_n^{(0)})\) and substituting these values into the right-hand side of the above equations, we get the first approximations \(x^{(1)} = (x_1^{(1)}, x_2^{(1)}, \ldots, x_n^{(1)})\). Continuing this process iteratively, we obtain a sequence of approximations \(\{x^{(k)}\}\) that converges to the exact solution of the system up to a given error tolerance as \(k \to \infty\).

**Applications:**
- Solving large systems of linear equations in fields such as engineering and physics.
- Particularly useful for sparse matrices.

**Advantages:**
- Simple to implement and understand.
- Can be parallelized, making it efficient for large systems.

**Disadvantages:**
- Convergence can be slow, especially for poorly conditioned matrices.
- May not converge for some types of matrices unless they are diagonally dominant.



### b. Gauss-Seidel Iterative Method

In the Gauss-Seidel method, the value of the variables is modified as soon as a new value is evaluated. For instance, in the Jacobi method, the value of \(x_i^{(k)}\) is not modified until the \((k + 1)\)th iteration, whereas in the Gauss-Seidel method, the value of \(x_i^{(k)}\) changes in the \(k\)th iteration.

**Applications:**
- Often used in engineering applications where large systems arise.
- Effective for systems modeled by finite difference equations.

**Advantages:**
- Typically converges faster than the Jacobi method.
- Efficient for large, sparse matrices.

**Disadvantages:**
- Convergence is not guaranteed for all types of matrices.
- Requires initial values to be reasonably close to the solution for better convergence.



### c. Gauss Elimination

**Form the Augmented Matrix:**
- Combine the coefficient matrix \(A\) and the constant vector \(b\) into an augmented matrix \([A | b]\).

**Forward Elimination:**
- Convert the augmented matrix into an upper triangular form using elementary row operations. This involves:
  - Selecting a pivot element (the first non-zero element in the current column).
  - Performing row operations to eliminate all entries below the pivot in that column.
  - Repeating for each column, moving down the rows.

**Row Operations:**
- Swap two rows.
- Multiply a row by a non-zero scalar.
- Add or subtract a multiple of one row to another row.

**Back Substitution:**
- Once the matrix is in upper triangular form, start from the last row and solve for the unknowns:
  - For the last row, you can directly solve for the corresponding variable.
  - Substitute back into the previous rows to find the remaining variables step by step.

**Applications:**
- Used extensively in engineering, physics, and computer science for solving systems of linear equations.
- Also applicable in operations research and numerical simulations.

**Advantages:**
- **Exact Solutions:** Provides a direct method to obtain exact solutions for systems of linear equations.
- **Wide Applicability:** Can be used for any size of linear systems.
- **Structured Approach:** Systematic method that is easy to follow and implement.

**Disadvantages:**
- **Computationally Intensive:** Requires \(O(n^3)\) operations for an \(n \times n\) system, making it inefficient for very large systems.
- **Numerical Stability:** Can lead to numerical instability and rounding errors, particularly for ill-conditioned matrices.
- **Requires Full Rank:** The method requires that the system has a unique solution (i.e., the matrix must be of full rank).



### d. Gauss-Jordan Elimination

**Form the Augmented Matrix:**
- Start with the augmented matrix \([A | b]\) that combines the coefficient matrix \(A\) and the constant vector \(b\).

**Forward Elimination to Row Echelon Form:**
- Convert the matrix to row echelon form (upper triangular form) using elementary row operations:
  - Select a pivot element in each column.
  - Eliminate all entries below the pivot by performing appropriate row operations.

**Back Elimination to Reduced Row Echelon Form:**
- Continue applying row operations to eliminate all entries above each pivot, transforming the matrix into reduced row echelon form (RREF):
  - The pivot elements in each column should be 1.
  - Each pivot column should have zeros in all other entries.

**Extract Solutions:**
- Once in RREF, read the solutions directly from the matrix. If the system has a unique solution, the matrix will have a pivot in every column corresponding to a variable.

**Applications:**
- Widely used in linear algebra for solving systems of linear equations.
- Applicable for finding the inverse of a matrix (if it exists).

**Advantages:**
- **Direct Solutions:** Provides a straightforward method to obtain the solutions of linear systems.
- **Unique and Consistent Solutions:** The RREF form indicates whether the system has a unique solution or is inconsistent.


#### <span style="font-size: 12px; font-weight: bold;">e. LU Factorization</span>  
**Overview:** LU factorization decomposes a square matrix A into a product of a lower triangular matrix L and an upper triangular matrix U.  

**Steps:**  
1. Start with a square matrix A.  
2. Initialize L as an identity matrix and U as a copy of A.  
3. Perform Gaussian elimination to convert U to upper triangular form while updating L with multipliers.  
4. The result is A = LU, where L contains the multipliers and U is upper triangular.  

**Applications:**  
- Solving linear systems with the same coefficient matrix.  
- Matrix inversion.  
- Calculating determinants.  

**Advantages:**  
- Efficient for multiple systems with the same matrix.  
- Simplifies the solution process into triangular systems.  
- More stable than some direct methods.  

**Disadvantages:**  
- Requires additional storage for L and U.  
- Can be unstable for singular or nearly singular matrices.  
- Partial pivoting may be needed for better stability, complicating the process.  
# Nonlinear Equations

Nonlinear equations are equations where the relationship between variables is not linear. These equations include higher powers of variables (e.g., X^2 or X^3, as well as trigonometric (e.g., sin(x), cos(x)), exponential, and logarithmic functions. In nonlinear equations, the dependent variable does not change in direct proportion to the independent variable.

For example, in the equation y = x^2 + 1, as x increases, y changes in a non-linear manner. The highest degree (power) of any term in the equation is called the degree of the equation, and the total number of solutions corresponds to this degree.

## General Form of non linear equation

A common representation of a nonlinear equation of degree 3 is:


ax^3 + bx^2 + cx + d = 0


where x and y are variables and a, b, c, and d are constants.

## Methods for Finding Roots of Nonlinear Equations

There are several numerical methods to find the roots of nonlinear equations. Some of the most common methods include:

### a. Bisection Method

The Bisection Method, also known as the binary chopping or half-interval method, works under the condition that if a function \(f(x)\) is real and continuous within an interval \(a < x < b\) and \(f(a)\) and \(f(b)\) are opposite sign, that is \(f(a) \* f(b) < 0\) then there is at least one real root in the interval between \(a\) and \(b\).

1. Start with an interval \(temp\_a, temp\_b\) where \(temp\_a = -Xmax\) & \(temp\_b = -Xmax + 1\) if \(f(temp\_a)\) and \(f(temp\_b)\) have opposite signs. Then these values are taken as \(a = temp\_a\) and \(b = temp\_b\) and the first root is calculated by following steps 2 to 5.
2. Compute the midpoint \(x_0\) as:

   X_0 = \(a + b\)/2
 
3. Check if \(f(X_0) = 0\). If yes, then \(X_0\) is the root.
4. If \(f(X_0) \* f(a) < 0\), replace \(b\) with \(X_0\) (narrowing the interval to \([a, x_0]\)).
5. If \(f(X_0) \* f(b) < 0\), replace \(a\) with \(x_0\) (narrowing the interval to \([X_0, b]\)).
6. Repeat this process until all roots are found. After finding the first root using step 1, calculate a new interval to find the remaining roots.

### b. False Position Method

The False Position Method improves upon the Bisection Method by introducing a more efficient formula to find the root. Instead of using the midpoint, it calculates an improved estimate using the function values at the interval endpoints. Here \(a\) and \(b\) is calculated as same way as bisection for every root.

The formula is:


X_0 = \(a \* f(b) - b \* f(a) \) \/ \( f(b) - f(a) \)


Steps are similar to the bisection method, but the formula for the root is recalculated using the above formula. This method tends to converge faster.

### 3. Newton-Raphson Method

The Newton-Raphson Method is a more efficient and faster method that starts with an initial guess and iteratively improves it. This method uses the derivative of the function to find successively better approximations to the root.

The formula is:


X_i = X_i-1 \- \( f(X_i-1) \/ f'(X_i-1) \)


Steps:

1. Start with an initial guess \(X_0\) where \(X_0 = \(a + b\) \/ 2\).
2. Compute the next approximation using the formula.
3. Repeat until the function value \(f(X_i)\) is sufficiently close to zero.

### d. Secant Method

The Secant Method is similar to the Newton-Raphson method but does not require the computation of the derivative. Instead, it approximates the derivative by using two previous points. The formula is:


X_n+1 = X_n - \( f(X_n) \* (X_n - X_n-1) \) \/ f(X_n) - f(X_n-1)

Steps:

1. Start with an initial guess X_n and X_n-1 where X_n = a and X_n-1 = b.
2. Compute the next approximation using the formula.
3. Repeat until the function value f(X_n+1) is sufficiently close to zero.

## Summary of Methods

| Method            | Efficiency | Approach                                |
|-------------------|------------|-----------------------------------------|
| Bisection Method  | Moderate   | Finds root by halving intervals         |
| False Position    | Faster     | Improves bisection with weighted midpoint |
| Newton-Raphson    | Fast       | Uses derivative for rapid convergence   |
| Secant Method     | Fast       | Similar to Newton-Raphson; no derivative |

Each method has its strengths and is suited for different types of equations and requirements for speed and accuracy.

### <span style="font-size: 18px; font-weight: bold;">3. Solutions for Differential Equations</span>  
#### <span style="font-size: 12px; font-weight: bold;">a. Runge-Kutta Method</span>  
1. **Define Initial Conditions:** Start with the initial condition y(x0) = y0 and step size h.  
2. **Calculate the Four Slopes:** For each step xn, calculate the following slopes:  
   - k1=h*f(x(n),y(n))
   - k2=h*f(x(n)+(h/2),y(n)+(k1/2)) 
   - k3=h*f(x(n)+(h/2),y(n)+(k2/2))  
   - k4=h*f(x(n)+h,y(n)+k3)  

3. **Update Solution Estimate:** Use these slopes to estimate the next value:  
   - y(n+1) = y(n) + (1/6)*(k1 + 2*k2 + 2*k3 + k4)  

4. **Repeat for Each Interval:** Move to the next interval by setting x*n+1 = x*n + h and repeat the process.  

**Applications:**  
- Widely used in engineering, physics, and biological modeling for solving ODEs.  
- Useful in simulations, where high accuracy over many intervals is essential.  

**Advantages:**  
- **High Accuracy:** The fourth-order version is highly accurate for many ODEs.  
- **Balanced Efficiency:** Good trade-off between computational cost and accuracy.  

**Disadvantages:**  
- **Fixed Step Size:** Standard RK4 can be inefficient for problems that need variable step sizes for better accuracy.  
- **Computationally Intensive:** More calculations per step than simpler methods like Euler's method.  

### <span style="font-size: 18px; font-weight: bold;">4. Matrix Inversion</span>  
1. **LU Decomposition:** Decompose matrix A into A = LU, where:  
   - L is a lower triangular matrix with ones on its diagonal.  
   - U is an upper triangular matrix.  
   - This decomposition can be done using techniques like Gaussian elimination.  

2. **Set up Identity Matrix:** To find the inverse A^-1,let A.A^-1=I be the identity matrix.  

3. **Solve for Inverse Columns:**  
   - Solve L.Y = I to find matrix Y (intermediate matrix).  
   - Solve U.A^-1 = Y to get A^-1.  
   - These are solved by considering each column of I and finding the corresponding column in X.  

4. **Combine Solutions:** The columns from step 3 form the inverse matrix A^-1.  

**Applications:**  
- System of Linear Equations: Used for efficiently solving linear systems where multiple right-hand sides need to be solved.  
- Engineering and Science: Common in numerical simulations and modeling where inversion of large matrices is required.  

**Advantages:**  
- Efficient for Repeated Solutions: LU decomposition can be reused for multiple right-hand sides, saving computation.  
- Simplifies Calculations: Once L and U are computed, inverting becomes a set of straightforward linear solves.  

**Disadvantages:**  
- Not Always Stable: LU decomposition may be numerically unstable for certain matrices (e.g., singular or near-singular matrices).  
- Inapplicable for Singular Matrices: LU decomposition and inversion cannot be used if A is singular (determinant is zero).  

## <span style="font-size: 18px; font-weight: bold;">Conclusion</span>  
This Numerical Methods Solver provides essential tools for solving linear and non-linear equations, differential equations, and matrix inversion, facilitating various applications. Each algorithm has strengths and weaknesses, making it crucial to choose the appropriate method based on the specific problem requirements.  
