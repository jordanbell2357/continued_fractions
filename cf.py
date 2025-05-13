import decimal
from fractions import Fraction
import numpy as np

def decimal_to_cf(x: decimal.Decimal, num_terms: int = 20) -> list[int]:
    """Compute the first `num_terms` partial quotients of the real number x."""

    cf = []
    for _ in range(num_terms):
        a = int(x)          # floor
        cf.append(a)
        frac_part = x - a   # x - floor(x)
        if frac_part == 0:
            # x is an integer, so the continued fraction terminates
            break
        x = 1 / frac_part   # 1 / (x - floor(x))

    return cf


def cf_fraction(x: Fraction) -> list[int]:
    """Compute the partial quotients of the fraction p / q."""

    cf = []
    while True:
        a = int(x)          # floor
        cf.append(a)
        frac_part = x - a   # x - floor(x)
        if frac_part == 0:
            # x is an integer, so the continued fraction terminates
            break
        x = 1 / frac_part   # 1 / (x - floor(x))

    return cf


def convergents_from_partial_quotients(a: list[int]) -> list[Fraction]:
    N = len(a)
    if N == 0:
        return []
    
    p_prev2, p_prev1 = 0, 1
    q_prev2, q_prev1 = 1, 0
    convergent_list = []
    
    for k in range(N):
        p_k = a[k] * p_prev1 + p_prev2
        q_k = a[k] * q_prev1 + q_prev2
        convergent_list.append(Fraction(p_k, q_k))

        p_prev2, p_prev1 = p_prev1, p_k
        q_prev2, q_prev1 = q_prev1, q_k
    
    return convergent_list


def convergents_from_partial_quotients_np(a: list[int]) -> list[Fraction]:
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

    fraction_list = [Fraction(int(p[k]), int(q[k])) for k in range(N)]

    return fraction_list


def determinants_from_partial_quotients_np(a: list[int]) -> list[int]:
    N = len(a)

    A = np.empty(N * 4, dtype='uint').reshape(N, 2, 2)
    for k in range(N):
        A[k] = np.array([[a[k], 1], [1, 0]])

    P = np.empty(N * 4, dtype='uint').reshape(N, 2, 2)
    P[0] = A[0]
    for k in range(1, N):
        P[k] = P[k-1] @ A[k]

    return A, [int(np.linalg.det(A[k])) for k in range(N)]


# Example usage:
if __name__ == "__main__":
    # Example 1: gcd
    p = 18
    q = 4
    partial_quotients = cf_fraction(Fraction(p, q))
    print(f"Partial quotients for {Fraction(p, q)}")
    print(partial_quotients)
    print(f"gcd({p}, {q}) =", partial_quotients[-1])

    print()

    # Example 2: Bezout coefficients
    p = 18
    q = 4
    partial_quotients = cf_fraction(Fraction(p, q))
    convergents = convergents_from_partial_quotients(partial_quotients)
    n = len(convergents) - 1
    p_second_last = convergents[n - 1].numerator
    q_second_last = convergents[n - 1].denominator
    s = (-1)**(n + 1) * q_second_last
    t = (-1)**n * p_second_last
    print(f"Bezout coefficients of s * {p} + t * {q} = gcd({p}, {q})")
    print(f"s * p + t * q = {s} * {p} + {t} * {q} = {s * p + t * q}")
    print(f"gcd({p}, {q}) = {s * p + t * q}")
    
    print()

    # Example 3: sqrt(n)
    d = 7
    num_terms = 20
    decimal.getcontext().prec = 50

    sqrtd = decimal.Decimal(d).sqrt()
    cf_sqrtd = decimal_to_cf(sqrtd, num_terms=num_terms)
    print(f"Partial quotients for sqrt({d})")
    print(f"[a_0, ..., a_{num_terms -1}]", '=', cf_sqrtd)

    print()

    # Example 4: Convergents
    d = 7
    num_terms = 20
    decimal.getcontext().prec = 50

    sqrtd = decimal.Decimal(d).sqrt()
    cf_sqrtd = decimal_to_cf(sqrtd, num_terms=num_terms)

    print(f"Convergents for sqrt({d})")
    convergents = convergents_from_partial_quotients(cf_sqrtd)
    print(convergents)

    print("Using NumPy matrix multiplication")
    convergents_np = convergents_from_partial_quotients_np(cf_sqrtd)
    print(convergents_np)

    print()

    # Example 5: Determinants
    d = 7
    num_terms = 20
    decimal.getcontext().prec = 50

    sqrtd = decimal.Decimal(d).sqrt()
    cf_sqrtd = decimal_to_cf(sqrtd, num_terms=num_terms)
    convergents_sqrtd = convergents_from_partial_quotients(cf_sqrtd)
    matrices, determinants = determinants_from_partial_quotients_np(convergents_sqrtd)
    print("Determinants of 2x2 matrices of convergents")
    print(determinants)

    print()

    # Example 6: Cotes continued fraction for e
    num_terms = 40
    decimal.getcontext().prec = 50
    e = decimal.Decimal(1).exp()
    cf_e = decimal_to_cf(e, num_terms=num_terms)
    print("Partial quotients of e")
    print(cf_e)
