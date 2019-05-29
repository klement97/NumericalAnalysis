import numpy as np

a = np.arange(15).reshape(3, 5)

# print(a)
# print("Shape: ", a.shape)
# print("Dimension: ", a.ndim)
# print("Size", a.size)
# print("Data type: ", a.dtype)
# print("Item size", a.itemsize)
# print("Data", a.data)
# print("Sum of columns:\n", a.sum(axis=0))
# print("Sum of rows:\n", a.sum(1))
# print("Min by columns: ", a.min(axis=0))
# print("Min by rows: ", a.min(1))
# print("Max by columns: ", a.max(0))
# print("Max by rows: ", a.max(1))
# print("Cumulative sum by columns:\n", a.cumsum(0))
# print("Cumulative sum by rows:\n", a.cumsum(1))


def gauss_nn(a, b):
    n = a.shape[0]  # number of rows
    a = np.array(a)
    b = np.array(b)
    g = np.eye(n, n)

    if a[0, 0] == 0:
        print("Enter a valid matrix. Element in the first row and first column is zero.")
        return -1

    for k in range(0, n-2):
        for i in range(k+1, n-1):
            g[i, k] = a[i, k] / a[k, k]
            for j in range(k+1, n):
                a[i, j] = a[i, j] - g[i, k] * a[k, j]

    return a


a = np.array([[1, 2, 3], [3, 7, 8], [4, 5, -1]])
b = gauss_nn(a, 1)
print(b)
