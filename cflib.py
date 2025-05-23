import functools as ft
from fractions import Fraction
import decimal
from decimal import Decimal
import math
import typing


def stern_diatomic(n: int) -> int:
    @ft.lru_cache
    def stern_diatomic_recurse(n: int) -> int:
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


def euclid(x: int, y: int) -> tuple[list[int], list[int], list[int], list[int]]:
    r0 = 1 * x - 0 * y
    r  = 0 * x - (-1) * y

    a0, a = 1, 0
    b0, b = 0, 1

    q_list = []
    r_list = [r0]
    a_list = [a0]
    b_list = [b0]
    
    while r != 0:
        q = r0 // r
        q_list.append(q)
        
        r0, r = r, r0 - q * r
        a0, a = a, a0 - q * a
        b0, b = b, b0 - q * b
        
        r_list.append(r0)
        a_list.append(a0)
        b_list.append(b0)

    r_list.append(r)
    a_list.append(a)
    b_list.append(b)
    
    return q_list, r_list, a_list, b_list


def decimal_to_cf(x: Decimal, num_terms: int = 20) -> list[int]:
    cf = []
    for _ in range(num_terms):
        a = int(x)
        cf.append(a)
        frac_part = x - a
        if frac_part == 0:
            break
        x = 1 / frac_part
    return cf


def gauss_transformation(x: Decimal) -> Decimal:
    return (1 / x) % 1 if x != 0 else 0


def distance_to_nearest_integer(x: Decimal) -> Decimal:
    floor = math.floor(x)
    ceil = math.ceil(x)
    return min(x - floor, ceil - x)


def distance_to_nearest_integer_sum(x: Decimal, n: int) -> Decimal:
    return sum(distance_to_nearest_integer(k * x) for k in range(1, n + 1))


def distance_to_nearest_integer_reciprocal_sum(x: Decimal, n: int) -> Decimal:
    return sum(1 / distance_to_nearest_integer(k * x) for k in range(1, n + 1))


def decimal_to_cf_unit_interval(x: Decimal, num_terms: int = 20) -> list[int]:
    x_list = []
    for _ in range(num_terms - 1):
        x_list.append(x)
        x = gauss_transformation(x)
    a_list = [0] + [int(1 / x) for x in x_list]
    return a_list


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


def fraction_tuple_to_convergent_list(fraction_tuple: tuple[int, int]) -> list[tuple[int, int]]:
    cf = fraction_tuple_to_cf(fraction_tuple)
    convergent_list = cf_to_convergent_list(cf)
    return convergent_list


def cf_to_fraction_tuple(cf: list[int]) -> tuple[int, int]:
    convergent_list = cf_to_convergent_list(cf)
    fraction_tuple = convergent_list[-1]
    return fraction_tuple


def modular_inverse(n: int, x: int) -> int:
    if math.gcd(n, x) > 1:
        return None
    else:
        _, _, _, b_list = euclid(n, x)
        bezout_coefficient = b_list[-2]
        return bezout_coefficient % n


class EEA:
    def __init__(self: typing.Self, x: int, y: int) -> None:
        self.x = x
        self.y = y
        q_list, r_list, a_list, b_list = euclid(self.x, self.y)
        self.q_list = q_list
        self.r_list = r_list
        self.a_list = a_list
        self.b_list = b_list
        self.x_list = [(-1) ** (k + 1) * t for k, t in enumerate(self.b_list)]
        self.y_list = [(-1) ** k * s for k, s in enumerate(self.a_list)]
        self.gcd = self.r_list[-2]
        self.bezout_x = self.a_list[-2]
        self.bezout_y = self.b_list[-2]
        self.cf = self.q_list
        self.convergent_list = list(zip(self.x_list, self.y_list))[2:]

    def __repr__(self) -> str:
        return f"{type(self).__name__}(x={self.x}, y={self.y})"
    
    def __len__(self):
        return len(self.q_list)


# Example usage:
if __name__ == "__main__":
    # Example: Euclid
    x, y = 240, 18

    eea = EEA(x, y)
    q_list, r_list, a_list, b_list = euclid(x, y)
    assert eea.q_list == q_list
    assert eea.r_list == r_list
    assert eea.a_list == a_list
    assert eea.b_list == b_list

    eea = EEA(x, y)
    assert len(eea) == len(eea.q_list)

    eea = EEA(x, y)
    assert eea.gcd == math.gcd(x, y)

    eea = EEA(x, y)
    assert eea.bezout_x * x + eea.bezout_y * y == math.gcd(x, y)

    eea = EEA(x, y)
    assert eea.convergent_list[-1] == Fraction(x, y).as_integer_ratio()

    eea = EEA(x, y)
    assert eea.cf == fraction_tuple_to_cf((x, y))

    eea = EEA(x, y)
    assert eea.convergent_list == fraction_tuple_to_convergent_list((x, y))


    # Example: modular inverse
    n = 74
    x = 13
    assert math.gcd(n, x) > 1 or modular_inverse(n, x) * x % n == 1


    # Example: Gauss transformation
    x = Decimal(0.34534267)
    num_terms = 5
    n = 20
    precision = 20
    decimal.getcontext().prec = precision

    cf = decimal_to_cf(x, num_terms)
    a_list = decimal_to_cf_unit_interval(x, num_terms)
    assert x < 0 or x >= 1 or a_list == cf

    # Example: distance_to_nearest_integer_sum
    assert distance_to_nearest_integer_sum(x, n) >= n / 4 - math.log(n) ** 2
    assert distance_to_nearest_integer_sum(x, n) <= n / 4 + math.log(n) ** 2

    # Example: distance_to_nearest_integer_reciprocal_sum
    assert distance_to_nearest_integer_reciprocal_sum(x, n) / (Decimal(n) * Decimal(n).ln()) <= 3

    # Example: Decimal Convergents
    d = 7
    num_terms = 20
    precision = 20
    decimal.getcontext().prec = precision
    sqrtd = Decimal(d).sqrt()
    cf_sqrtd = decimal_to_cf(sqrtd, num_terms=num_terms)
    print(f"Convergents for sqrt({d})")
    convergents = cf_to_convergent_list(cf_sqrtd)
    print(convergents)

    print()

    # Example: Decimal Cotes continued fraction for e
    num_terms = 40
    precision = 20
    decimal.getcontext().prec = precision
    e = Decimal(1).exp()
    cf_e = decimal_to_cf(e, num_terms=num_terms)
    print("Partial quotients of e")
    print(cf_e)

    print()

    # Example: gcd
    p = 18
    q = 4

    partial_quotients = fraction_tuple_to_cf((p, q))
    print(f"Partial quotients for {Fraction(p, q)}")
    print(partial_quotients)
    print(f"gcd({p}, {q}) =", partial_quotients[-1])

    print()

    # Example: Bezout coefficients
    partial_quotients = fraction_tuple_to_cf((p, q))
    convergents = cf_to_convergent_list(partial_quotients)
    n = len(convergents) - 1
    p_second_last, q_second_last = convergents[n - 1]
    s = (-1)**(n + 1) * q_second_last
    t = (-1)**n * p_second_last
    print(f"Bezout coefficients of s * {p} + t * {q} = gcd({p}, {q})")
    print(f"s * p + t * q = {s} * {p} + {t} * {q} = {s * p + t * q}")
    print(f"gcd({p}, {q}) = {s * p + t * q}")
    
    print()

    # Example: sqrt(n)
    d = 7
    num_terms = 20
    precision = 20
    decimal.getcontext().prec = precision
    sqrtd = Decimal(d).sqrt()
    cf_sqrtd = decimal_to_cf(sqrtd, num_terms=num_terms)
    print(f"Partial quotients for sqrt({d})")
    print(f"[a_0, ..., a_{num_terms -1}]", '=', cf_sqrtd)
