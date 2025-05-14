import itertools as it
import functools as ft
import operator
from fractions import Fraction
import decimal
import math

def fraction_tuple_to_cf(fraction_tuple: tuple[int, int]) -> list[int]:
    x = Fraction(*fraction_tuple)
    cf = []
    while True:
        a = int(x) # floor
        cf.append(a)
        frac_part = x - a # x - floor(x)
        if frac_part == 0:
            # x is an integer, so the continued fraction terminates
            break
        x = 1 / frac_part # 1 / (x - floor(x))
    return cf

def decimal_to_cf(x: decimal.Decimal, num_terms: int = 20) -> list[int]:
    cf = []
    for _ in range(num_terms):
        a = int(x)
        cf.append(a)
        frac_part = x - a
        if frac_part == 0:
            break
        x = 1 / frac_part
    return cf

def cf_to_convergent_list(cf: list[int]) -> list[tuple[int, int]]:
    N = len(cf)
    if N == 0:
        return []
    p_prev2, p_prev1 = 0, 1
    q_prev2, q_prev1 = 1, 0
    convergent_list = []
    for k in range(N):
        p_k = cf[k] * p_prev1 + p_prev2
        q_k = cf[k] * q_prev1 + q_prev2
        convergent_list.append((p_k, q_k))
        p_prev2, p_prev1 = p_prev1, p_k
        q_prev2, q_prev1 = q_prev1, q_k
    return convergent_list

def mediant(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (x[0] + y[0], x[1] + y[1])

def mediant_left_inverse(x: tuple[int, int], m: tuple[int, int]):
    return (m[0] - x[0], m[1] - x[1])

def mediant_right_inverse(y: tuple[int, int], m: tuple[int, int]):
    return (m[0] - y[0], m[1] - y[1])


def stern_diatomic(n):
    @ft.lru_cache
    def stern_diatomic_recurse(n):
        if n == 0:
            return 0
        elif n == 1:
            return 1
        else:
            if n % 2 == 0:
                return stern_diatomic_recurse(n // 2)
            else:
                m = (n - 1) // 2
                return stern_diatomic_recurse(m) + stern_diatomic_recurse(m + 1)
    return stern_diatomic_recurse(n)


if __name__ == "__main__":
    print(mediant((0, 1), (1, 2)))
    print(mediant_left_inverse((0, 1), (1, 3)))
    print(mediant_right_inverse((1, 2), (1, 3)))


# # Example usage:
# if __name__ == "__main__":
#     # Example 1: gcd
#     p = 18
#     q = 4

#     partial_quotients = fraction_tuple_to_cf((p, q))
#     print(f"Partial quotients for {Fraction(p, q)}")
#     print(partial_quotients)
#     print(f"gcd({p}, {q}) =", partial_quotients[-1])

#     print()

#     # Example 2: Bezout coefficients
#     partial_quotients = fraction_tuple_to_cf((p, q))
#     convergents = cf_to_convergent_list(partial_quotients)
#     n = len(convergents) - 1
#     p_second_last, q_second_last = convergents[n - 1]
#     s = (-1)**(n + 1) * q_second_last
#     t = (-1)**n * p_second_last
#     print(f"Bezout coefficients of s * {p} + t * {q} = gcd({p}, {q})")
#     print(f"s * p + t * q = {s} * {p} + {t} * {q} = {s * p + t * q}")
#     print(f"gcd({p}, {q}) = {s * p + t * q}")
    
#     print()

#     # Example 3: sqrt(n)
#     d = 7
#     num_terms = 20
#     decimal.getcontext().prec = 50

#     sqrtd = decimal.Decimal(d).sqrt()
#     cf_sqrtd = decimal_to_cf(sqrtd, num_terms=num_terms)
#     print(f"Partial quotients for sqrt({d})")
#     print(f"[a_0, ..., a_{num_terms -1}]", '=', cf_sqrtd)

#     print()

#     # Example 4: Convergents
#     d = 7
#     num_terms = 20
#     decimal.getcontext().prec = 50
#     sqrtd = decimal.Decimal(d).sqrt()
#     cf_sqrtd = decimal_to_cf(sqrtd, num_terms=num_terms)
#     print(f"Convergents for sqrt({d})")
#     convergents = cf_to_convergent_list(cf_sqrtd)
#     print(convergents)

#     # Example 5: Cotes continued fraction for e
#     num_terms = 40
#     decimal.getcontext().prec = 50
#     e = decimal.Decimal(1).exp()
#     cf_e = decimal_to_cf(e, num_terms=num_terms)
#     print("Partial quotients of e")
#     print(cf_e)