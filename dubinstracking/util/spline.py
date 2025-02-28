import numpy as np


def fit(p: np.ndarray, type: int = 0):
    """
    fit a 3rd order polynomial on the knots and points p

    Parameters
    ----------
    p: np.ndarray
        a Nx2 array of knots and points with the first column being the knots and the 2nd column being the points

    type: int
        0 starts and ends with zero acceleration, 1 started and ends with zero speed, 2 is periodic constraints

    Returns
    -------
        np.ndarray
        a one dimensional array of coefficients for the spline fit
    """
    k = 3 + 1
    N = p.shape[0] - 1
    n = (p.shape[0] - 1) * k
    A = np.zeros([n, n])
    # first system with zero velocity boundary condition
    A[:2, :4] = [
        [p[0, 0] ** 3, p[0, 0] ** 2, p[0, 0], 1],
        [p[1, 0] ** 3, p[1, 0] ** 2, p[1, 0], 1],
    ]
    # starts and ends with zero acceleration
    if type == 0:
        A[2, :4] = [6 * p[0, 0], 2, 0, 0]
        A[3, (n - k):n] = [6 * p[N, 0], 2, 0, 0]
    # starts and ends with zero speed
    elif type == 1:
        A[2, :k] = [3 * p[0, 0] ** 2, 2 * p[0, 0], 1, 0]
        A[3, (n - k):n] = [3 * p[N, 0] ** 2, 2 * p[N, 0], 1, 0]
    # periodic
    elif type == 2:
        A[2, :k] = [3 * p[0, 0] ** 2, 2 * p[0, 0], 1, 0]
        A[2, (n - k):n] = [-3 * p[N, 0] ** 2, -2 * p[N, 0], -1, 0]
        A[3, :k] = [6 * p[0, 0], 2, 0, 0]
        A[3, (n - k):n] = [-6 * p[N, 0], -2, 0, 0]

    B = np.zeros(n)
    B[:k] = [
        p[0, 1],
        p[1, 1],
        0,
        0
    ]
    for i in range(1, N):
        # system matrix
        A[(i * k):((i + 1) * k), (i * k):((i + 1) * k)] = [
            [p[i, 0] ** 3, p[i, 0] ** 2, p[i, 0], 1],
            [p[i + 1, 0] ** 3, p[i + 1, 0] ** 2, p[i + 1, 0], 1],
            [-3 * p[i, 0] ** 2, -2 * p[i, 0], -1, 0],
            [-6 * p[i, 0], -2, 0, 0]
        ]
        A[(i * k + 2):((i + 1) * k), ((i - 1) * k):(i * k)] = [
            [3 * p[i, 0] ** 2, 2 * p[i, 0], 1, 0],
            [6 * p[i, 0], 2, 0, 0]
        ]
        B[(i * k):((i + 1) * k)] = [
            p[i, 1],
            p[i + 1, 1],
            0,
            0
        ]
    x = np.linalg.solve(A, B)
    return x


def fit2func(p: np.ndarray, c: np.ndarray):
    """
    create a function from a fit c of a 3rd order polynomial on the knots and points p

    Parameters
    ----------
    p: np.ndarray
        a Nx2 array of knots and points with the first column being the knots and the 2nd column being the points

    c: np.ndarray
        a 4n array of fit coefficients
    
    Returns
    -------
    (float) -> float
    """
    k = 4

    def f(t):
        i = 0
        while t > p[i + 1, 0]:
            i += 1
        r = c[i * k] * t ** 3 + c[i * k + 1] * t ** 2 + c[i * k + 2] * t + c[i * k + 3]
        return r
    return f


def fit2deriv(p, c):
    """
    create a function of the derivative of a fit c of a 3rd order polynomial on the knots and points p

    Parameters
    ----------
    p: np.ndarray
        a Nx2 array of knots and points with the first column being the knots and the 2nd column being the points

    c: np.ndarray
        a 4n array of fit coefficients
    
    Returns
    -------
    (float) -> float
    """
    k = 4
    def f(t):
        i = 0
        while t > p[i + 1, 0]:
            i += 1
        r = 3 * c[i * k] * t ** 2 + 2 * c[i * k + 1] * t + c[i * k + 2]
        return r
    return f


def fit2deriv2(p, c):
    """
    create a function of the second derivative of a fit c of a 3rd order polynomial on the knots and points p

    Parameters
    ----------
    p: np.ndarray
        a Nx2 array of knots and points with the first column being the knots and the 2nd column being the points

    c: np.ndarray
        a 4n array of fit coefficients
    
    Returns
    -------
    (float) -> float
    """
    k = 4

    def f(t):
        i = 0
        while t > p[i + 1, 0]:
            i += 1
        r = 6 * c[i * k] * t + 2 * c[i * k + 1]
        return r
    return f
