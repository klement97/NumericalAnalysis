import math


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
    This method finds a single root in an isolated segment using bisection method
    number of iterations required can be found by this formula
    ( (b-a) / 2**(n+1) ) < eps: where n is number of iterations required


    :param f: use lambda functions to evaluate your function
                all you need is math library and lambda functions
    :param a: bottom of the segment in which root is isolated
    :param b: top of the segment in which root is isolated
    :param eps: max value for error 0.00001 by default
    :param n: max number of iterations, 20 by default
    :return x: approximation of root
    """

    # checking if segment is valid
    if f(a) * f(b) < 0:
        for iteration in range(n):
            x = (a+b) / 2
            # checking if x is close enough to the root
            if f(x) <= eps:
                return x

            else:
                if f(x) < 0:
                    a = x
                else:
                    b = x
    else:
        print("There is no root inside the segment [{};{}]".format(a, b))


def fixed_point(f, x0, eps=1e-5, n=20):

    """
    Finding root with fixed point iteration method. User need to find the required function to use in.
    Function is found by the equation itself. Functions first derivative must be < 1
    for all x from the segment in order to converge to the root. Otherwise function will diverge.
    To prevent diverging there is a cap on number of iterations.
    Accelerating process using Aitken's schema.


    :param f: iterative function
    :param x0: first approximation,
    :param eps: max error,
    :param n: max number of iterations,
    :return x: approximation of the root
    """

    for i in range(n):
        x1 = f(x0)
        x2 = f(x1)
        a = (x0*x2 - x1**2) / (x2 - 2*x1 + x0)
        if f(a) <= eps:
            return a
        if math.fabs(x1 - x2) < eps * math.fabs(x1):
            return a
        else:
            x0 = a
    else:
        print("No convergence after {} iterations".format(n))
