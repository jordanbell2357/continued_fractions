from collections.abc import Generator
from fractions import Fraction
import decimal
from decimal import Decimal
import bisect
import itertools as it
import functools as ft
import operator
import math
import time


def reduced_fraction_mediant(x: Fraction, y: Fraction) -> Fraction:
    numerator = x.numerator + y.numerator
    denominator = x.denominator + y.denominator
    return Fraction(numerator, denominator)


def totient(n: int) -> int:
    return sum(math.gcd(n, k) == 1 for k in range(1, n + 1))


def farey_generator(n: int) -> Generator[Fraction]:
    a, b, c, d = 0, 1, 1, n

    yield Fraction(a, b)

    while 0 <= c <= n:
        k = (n + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
        yield Fraction(a, b)


def farey_list(n: int) -> list[Fraction]:
    return list(farey_generator(n))


def next_farey_list(farey_list_order_n: list[Fraction]) -> list[Fraction]:
    a_iter, b_iter = it.tee(farey_list_order_n)
    next(b_iter, None)
    farey_list_order_n_plus_one = ft.reduce(
        operator.iadd,
        [([a, reduced_fraction_mediant(a, b)] if b is not None and reduced_fraction_mediant(a, b).denominator <= n + 1 else [a])
         for a, b in it.zip_longest(a_iter, b_iter)],
        [])
    return farey_list_order_n_plus_one


def next_farey_generator(farey_generator_order_n: Generator[Fraction]) -> Generator[Fraction]:
    a_iter, b_iter = it.tee(farey_generator_order_n)
    next(b_iter, None)
    farey_generator_order_n_plus_one = it.chain(*[(iter([a, reduced_fraction_mediant(a, b)]) if b is not None and reduced_fraction_mediant(a, b).denominator <= n + 1 else iter([a]))
         for a, b in it.zip_longest(a_iter, b_iter)])
    return farey_generator_order_n_plus_one


def find_farey_generator_index(n: int, target_fraction: Fraction) -> int:
    farey_generator_order_n = farey_generator(n)
    for k, farey_fraction in enumerate(farey_generator_order_n):
        if target_fraction == farey_fraction:
            return k
    return None

def find_farey_list_index(n: int, target_fraction: Fraction) -> int:
    farey_list_order_n = farey_list(n)
    left_index = bisect.bisect_left(farey_list_order_n, target_fraction)
    right_index = bisect.bisect_right(farey_list_order_n, target_fraction)
    if left_index == right_index:
        return None
    else:
        return left_index


def farey_diff(n: int, k: int) -> Fraction:
    farey_list_order_n = farey_list(n)
    l = len(farey_list_order_n)
    farey_fraction = farey_list_order_n[k]
    return farey_fraction - Fraction(k, l)


def farey_abs_diff_sum_fraction(n: int) -> Fraction:
    farey_list_order_n = farey_list(n)
    l = len(farey_list_order_n)
    return sum(abs(f - Fraction(k, l)) for k, f in enumerate(farey_list_order_n))


def farey_abs_diff_sum_decimal(n: int, precision: int = 10) -> Decimal:
    decimal.getcontext().prec = precision
    farey_list_order_n = farey_list(n)
    l = len(farey_list_order_n)
    sum_decimal = sum(Decimal((abs(f - Fraction(k, l)).numerator)) / Decimal((abs(f - Fraction(k, l)).denominator)) \
                      for k, f in enumerate(farey_list_order_n))
    return sum_decimal


def farey_abs_diff_sum_float(n: int) -> float:
    farey_list_order_n = farey_list(n)
    l = len(farey_list_order_n)
    return sum(abs(f - k / l) for k, f in enumerate(farey_list_order_n))



if __name__ == "__main__":
    precision = 20
    n = 13
    k = 4
    target_fraction = Fraction(2, 17)
    decimal.getcontext().prec = precision

    start_time = time.perf_counter()
    fraction_sum = float(farey_abs_diff_sum_fraction(n))
    end_time = time.perf_counter()
    print("Exact value", fraction_sum, "Time", end_time - start_time)

    start_time = time.perf_counter()
    decimal_sum = float(farey_abs_diff_sum_decimal(n, precision))
    end_time = time.perf_counter()
    print("Decimal approximation", decimal_sum, "Time", end_time - start_time, "Approximation error", decimal_sum - fraction_sum)

    start_time = time.perf_counter()
    float_sum = farey_abs_diff_sum_float(n)
    end_time = time.perf_counter()
    print("Float approximation", float_sum, "Time", end_time - start_time, "Approximation error", float_sum - fraction_sum)

    farey_generator_order_n = farey_generator(n)
    farey_generator_order_n_plus_one = next_farey_generator(farey_generator_order_n)
    farey_list_order_n_plus_one = list(farey_generator_order_n_plus_one)
    assert farey_list_order_n_plus_one == farey_list(n + 1)

    farey_list_order_n = farey_list(n)
    assert next_farey_list(farey_list_order_n) == farey_list(n + 1)

    assert len(farey_list(n + 1)) == len(farey_list(n)) + totient(n + 1)

    farey_list_order_n = farey_list(n)
    assert k >= n or find_farey_list_index(n, farey_list_order_n[k]) == k

    farey_list_order_n = farey_list(n)
    farey_fraction = farey_list_order_n[k]
    assert find_farey_generator_index(n, farey_fraction) == find_farey_list_index(n, farey_fraction)

    assert find_farey_generator_index(n, target_fraction) is None and find_farey_list_index(n, target_fraction) is None \
        or find_farey_generator_index(n, target_fraction) == find_farey_list_index(n, target_fraction)




