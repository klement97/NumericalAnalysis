import numpy as np


a = np.arange(15).reshape(3, 5)

print(a)
print("Shape: ", a.shape)
print("Dimension: ", a.ndim)
print("Size", a.size)
print("Data type: ", a.dtype)
print("Item size", a.itemsize)
print("Data", a.data)
print("Sum of columns:\n", a.sum(axis=0))
print("Sum of rows:\n", a.sum(1))
print("Min by columns: ", a.min(axis=0))
print("Min by rows: ", a.min(1))
print("Max by columns: ", a.max(0))
print("Max by rows: ", a.max(1))
print("Cumulative sum by columns:\n", a.cumsum(0))
print("Cumulative sum by rows:\n", a.cumsum(1))
