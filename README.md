# continued_fractions

```python
import numpy as np

x = np.sqrt(3)
N = 10

value_array = np.empty(N)

t = x

for i in range(N):
    value_array[i] = t
    t = np.modf(np.reciprocal(t))[0]

a = np.floor(np.reciprocal(np.roll(value_array, -1)))
```

```python
import numpy as np

N = len(a) - 1 # subtract 1 because we used np.roll

A = np.empty(N * 4, dtype='uint').reshape(N, 2, 2)
for k in range(N):
    A[k] = np.array([[a[k], 1], [1, 0]])

P = np.empty(N * 4, dtype='uint').reshape(N, 2, 2)
P[0] = A[0]
for k in range(1, N):
    P[k] = P[k-1] @ A[k]
```

```python
from fractions import Fraction
from decimal import Decimal

convergent_fractions = [Fraction(*P[n][:,0]) for n in range(N)]
convergent_decimals = [Decimal(int(f.numerator)) / Decimal(int(f.denominator)) for f in convergent_fractions]
```

```python
import numpy as np

x = np.sqrt(3)
N = 20

value_array = np.empty(N, dtype='longdouble')

t = x

for i in range(N):
    value_array[i] = t
    t = np.modf(np.reciprocal(t))[0]

a = np.floor(np.reciprocal(np.roll(value_array, -1)))[:-1]

L = len(a)

A = np.empty(L * 4, dtype='uint').reshape(L, 2, 2)
for k in range(L):
    A[k] = np.array([[a[k], 1], [1, 0]])

P = np.empty(L * 4, dtype='uint').reshape(L, 2, 2)
P[0] = A[0]
for k in range(1, L):
    P[k] = P[k-1] @ A[k]

p = np.empty(L, dtype='uint')
q = np.empty(L, dtype='uint')

for k in range(L):
    p[k] = P[k][0,0]
    q[k] = P[k][1,0]

convergent_array = p / q

convergent_array - x
```
