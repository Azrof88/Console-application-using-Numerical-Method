# Nonlinear Equations

Nonlinear equations are equations where the relationship between variables is not linear. These equations include higher powers of variables (e.g., \(X^2\) or \(X^3\)), as well as trigonometric (e.g., sin(x), cos(x)), exponential, and logarithmic functions. In nonlinear equations, the dependent variable does not change in direct proportion to the independent variable.

For example, in the equation \(y = x^2 + 1\), as \(x\) increases, \(y\) changes in a non-linear manner. The highest degree (power) of any term in the equation is called the degree of the equation, and the total number of solutions corresponds to this degree.

## General Form of non linear equation

A common representation of a nonlinear equation of degree 3 is:


ax^3 + bx^2 + cx + d = 0


where \(x\) and \(y\) are variables and \(a\), \(b\), \(c\), and \(d\) are constants.

## Methods for Finding Roots of Nonlinear Equations

There are several numerical methods to find the roots of nonlinear equations. Some of the most common methods include:

### 1. Bisection Method

The Bisection Method, also known as the binary chopping or half-interval method, works under the condition that if a function \(f(x)\) is real and continuous within an interval \(a < x < b\) and \(f(a)\) and \(f(b)\) are opposite sign, that is \(f(a) \* f(b) < 0\) then there is at least one real root in the interval between \(a\) and \(b\).

1. Start with an interval \(temp\_a, temp\_b\) where \(temp\_a = -Xmax\) & \(temp\_b = -Xmax + 1\) if \(f(temp\_a)\) and \(f(temp\_b)\) have opposite signs. Then these values are taken as \(a = temp\_a\) and \(b = temp\_b\) and the first root is calculated by following steps 2 to 5.
2. Compute the midpoint \(x_0\) as:

   X_0 = \(a + b\)/2
 
3. Check if \(f(X_0) = 0\). If yes, then \(X_0\) is the root.
4. If \(f(X_0) \* f(a) < 0\), replace \(b\) with \(X_0\) (narrowing the interval to \([a, x_0]\)).
5. If \(f(X_0) \* f(b) < 0\), replace \(a\) with \(x_0\) (narrowing the interval to \([X_0, b]\)).
6. Repeat this process until all roots are found. After finding the first root using step 1, calculate a new interval to find the remaining roots.

### 2. False Position Method

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

### 4. Secant Method

The Secant Method is similar to the Newton-Raphson method but does not require the computation of the derivative. Instead, it approximates the derivative by using two previous points. The formula is:


X_n+1 = X_n - \( f(X_n) \* (X_n - X_n-1) \) \/ f(X_n) - f(X_n-1)

Steps:

1. Start with an initial guess \(X_n\) and \(X_n-1\) where \(X_n = a\) and \(X_n-1 = b\).
2. Compute the next approximation using the formula.
3. Repeat until the function value \(f(X_n+1)\) is sufficiently close to zero.

## Summary of Methods

| Method            | Efficiency | Approach                                |
|-------------------|------------|-----------------------------------------|
| Bisection Method  | Moderate   | Finds root by halving intervals         |
| False Position    | Faster     | Improves bisection with weighted midpoint |
| Newton-Raphson    | Fast       | Uses derivative for rapid convergence   |
| Secant Method     | Fast       | Similar to Newton-Raphson; no derivative |

Each method has its strengths and is suited for different types of equations and requirements for speed and accuracy.
