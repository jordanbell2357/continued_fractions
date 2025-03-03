# continued_fractions

```python
from fractions import Fraction
import numpy as np

a = [1, 1, 1, 1, 1, 1, 1, 1]
N = len(a)

A = np.empty(N * 4, dtype='uint').reshape(N, 2, 2)
for k in range(N):
    A[k] = np.array([[a[k], 1], [1, 0]])

P = np.empty(N * 4, dtype='uint').reshape(N, 2, 2)
P[0] = A[0]
for k in range(1, N):
    P[k] = P[k-1] @ A[k]

Fraction(*P[N - 1][:,0])
```
