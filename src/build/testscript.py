import struker
import numpy as np
string = struker.DiracStringKernel()
a = ["hi", "hello", "bye", "hi"]
K = string.computeKernelMatrix(a, {})
"""for i in range(K.shape[0]):
    K[i, i] = -1*i
    if i < K.shape[0] - 1:
        K[i, i + 1] = 1.0
        K[i + 1, i] = 1.0
K[3, 0] = -5.0
K[0, 3] = -5.0"""
Knorm = string.normaliseHilbert(K)
print(K)
print("\n")
print(np.round(Knorm))