"""
Solving equations by approximating roots using Approximation Methods from Numerical Analysis

    - Bisection method
    - Methods of type Xn+1 = F(Xn)
    - Lagrange's method
    - Newton's method
    - Combined method
    - Secant method
    - Dekker-Brent method

"""


def bisection(f, a, b, eps=1e-5, n=20):

    """
     function - use lambda functions to evaluate your function
                all you need is math library and lambda functions
     a, b - segment in which root is isolated
     eps - max value for error 0.00001 by default
     n - max number of iterations, 20 by default

     returns x: approximation of root
    """

    # checking if segment is valid
    if f(a) * f(b) < 0:
        for iteration in range(n):
            x = (a+b) / 2  # half of the segment
            if f(x) <= eps:
                return x
            else:
                if f(x) < 0:
                    a = x
                else:
                    b = x