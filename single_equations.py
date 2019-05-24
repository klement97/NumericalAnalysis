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

import math


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
            if math.fabs(f(x)) <= eps:
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
        if math.fabs(f(a)) <= eps:
            return a
        if math.fabs(x1 - x2) < eps * math.fabs(x1):
            return a
        else:
            x0 = a
    else:
        print("No convergence after {} iterations".format(n))


def lagrange(f, a, b, eps=1e-5, n=30):

    """
    Approximating root using Lagrange's Method.
    Root must be isolated in the segment [a,b]
    This method is faster than bisection and fixed-point most of the times.
    Method stand on selecting two points inside the segment [a,b] as a base approximation and
    another fixed point to use in formula and finding better approximations using Lagrange's formula.
    We will select:
                    x0 = a
                    c = b

    :param f: function
    :param a: bottom of the segment
    :param b: top of the segment
    :param eps: max error
    :param n: max number of iterations
    :return x: approximation of root
    """
    x0 = a
    c = b
    for i in range(n):
        x = x0 - (f(x0) / (f(x0) - f(c))) * (x0 - c)
        if math.fabs(f(x)) <= eps:
            return x
        else:
            x0 = x


def newton(f, df, x0, eps=1e-5, n=10):
    """
    Newton's method for finding roots of an equation.
    Approximation found by this method is doubled in precision after every iteration.
    This is the fastest method in most cases.
    Only disadvantage of this method is that it requires first derivative of the function.

    :param f: function ,
    :param df: first derivative of function,
    :param x0: top of the segment,
    :param eps: max error,
    :param n: max number of iterations,
    :return x: approximation found by iterative formula.
    """

    for i in range(n):  # first critter: max number of iterations
        x = x0 - f(x0)/df(x0)

        if math.fabs(f(x)) <= eps:  # second critter: function at this point is close enough to zero
            return x
        else:
            x0 = x


def combined(f, df, a, b, eps=1e-5, n=10):
    """
    This method combines three methods: newton, lagrange and bisection.
    Using lagrange to approximate root from left,
    newton to approximate root from right
    and bisection to take half of the two approximations.
    Usually newton goes faster than lagrange so taking their mean is faster
    than using only one of the methods.

    :param f: function
    :param df: first derivative of function which is used in newtons method
    :param a: bottom of the segment
    :param b: top of the segment
    :param eps: max error: 1e-5 by default
    :param n: max number of iterations: 10 by default
    :return x: approximation with specified error
    """

    if f(a) * f(b) < 0:
        for i in range(n):
            x_l = a - (f(a) / (f(a) - f(b))) * (a - b)
            x_n = b - f(b)/df(b)
            x = (x_l + x_n) / 2

            if math.fabs(f(x)) <= eps:
                return x
            else:
                a = x_l
                b = x_n
    else:
        print("There is no root inside this segment.")
