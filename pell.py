import math
from fractions import Fraction
import numpy as np


def sqrt_continued_fraction(d: int) -> tuple[list[int], list[int]]:
    a0 = int(math.isqrt(d))  # floor(sqrt(d))
    if a0 ** 2 == d:
        # d is a square, so sqrt(d) is an integer and
        # the continued fraction has length 1
        return ([a0], [])

    # Initialize (M, D, a)
    M = 0
    D = 1
    a = a0

    expansion = [a0]
    
    # Keep track of triples (M, D, a) to detect repetition
    triples_dict = {}
    # Map: (M, D, a) -> index in the expansion

    i = 0
    while True:
        M = D * a - M
        D = (d - M**2) // D
        a = (a0 + M) // D
        
        i += 1

        # If we see a triple (M, D, a) twice, we've found the start of the period
        if (M, D, a) in triples_dict:
            # We now know the block from triples_dict[(M,D,a)] to i-1 is the repeating cycle
            start_period = triples_dict[(M, D, a)]
            non_periodic_part = expansion[:start_period]
            periodic_part = expansion[start_period:]
            return non_periodic_part, periodic_part
        else:
            triples_dict[(M, D, a)] = i
            expansion.append(a)


def convergents_from_partial_quotients(a):
    N = len(a)

    A = np.empty(N * 4, dtype='uint').reshape(N, 2, 2)
    for k in range(N):
        A[k] = np.array([[a[k], 1], [1, 0]])

    P = np.empty(N * 4, dtype='uint').reshape(N, 2, 2)
    P[0] = A[0]
    for k in range(1, N):
        P[k] = P[k-1] @ A[k]

    p = np.empty(N, dtype='uint')
    q = np.empty(N, dtype='uint')

    for k in range(N):
        p[k] = P[k][0,0]
        q[k] = P[k][1,0]

    return p, q

# x^2 - d y^2 = 1
# x https://oeis.org/A002350 and y https://oeis.org/A002349
# Example usage:
if __name__ == "__main__":
    n = 19
    non_periodic_part, periodic_part = sqrt_continued_fraction(n)
    print(non_periodic_part, periodic_part)
    r = len(periodic_part)
    r = r if r % 2 == 0 else 2 * r
    partial_quotients = non_periodic_part + 2 * periodic_part
    p, q = convergents_from_partial_quotients(partial_quotients)
    print(p[r-1], q[r-1])

    convergent_fraction_list = [Fraction(p[k], q[k]) for k in range(len(p))]
    print(convergent_fraction_list)
