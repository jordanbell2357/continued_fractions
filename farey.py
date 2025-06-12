from collections.abc import Generator
from fractions import Fraction
import decimal
from decimal import Decimal
import bisect
import itertools as it
import functools as ft
import operator
import math
import cmath
import time
import typing

import prime_numbers


def find_ordered_list_index(ordered_list: list[typing.Any], target_item: typing.Any) -> int:
    left_index = bisect.bisect_left(ordered_list, target_item)
    right_index = bisect.bisect_right(ordered_list, target_item)
    if left_index == right_index:
        return None
    else:
        return left_index

def check_ordered_list_membership(ordered_list: list[typing.Any], target_item: typing.Any) -> int:
    left_index = bisect.bisect_left(ordered_list, target_item)
    right_index = bisect.bisect_right(ordered_list, target_item)
    if left_index == right_index:
        return False
    else:
        return True


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


def find_farey_index(n: int, target_fraction: Fraction) -> int:
    farey_generator_order_n = farey_generator(n)
    for k, farey_fraction in enumerate(farey_generator_order_n):
        if target_fraction == farey_fraction:
            return k
    return None


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
    return math.fsum(abs(f - k / l) for k, f in enumerate(farey_list_order_n))


def mertens_function(n: int) -> int:
    farey_list_order_n = farey_list(n)
    exponential_sum = -1 + sum(cmath.exp(complex(0, 2 * math.pi * a)) for a in farey_list_order_n)
    return round(exponential_sum.real)    


class Farey:
    def __init__(self, n: int) -> None:
        self.n = n
        self.farey_list = farey_list(self.n)
        self.num_terms = len(self.farey_list)
        exponential_sum = -1 + sum(cmath.exp(complex(0, 2 * math.pi * a)) for a in self.farey_list)
        self.mertens_function = round(exponential_sum.real)
    
    def __repr__(self) -> str:
        return f"{type(self).__name__}(n={self.n})"
    
    def __eq__(self, other: typing.Self) -> bool:
        return self.n == other.n
    
    def __lt__(self, other: typing.Self) -> bool:
        return self.n < other.n
    
    # __iter__ implements __contains__ but we choose to implement __contains__ independently
    def __contains__(self, item: Fraction) -> bool:
        return check_ordered_list_membership(self.farey_list, item)
    
    def __iter__(self):
        return iter(self.farey_list)


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
    target_fraction = farey_list_order_n[k]
    assert find_ordered_list_index(farey_list_order_n, target_fraction) == find_farey_index(n, target_fraction)

    farey_list_order_n = farey_list(n)
    assert k >= n or find_farey_index(n, farey_list_order_n[k]) == k

    assert find_farey_index(n, target_fraction) is None and find_farey_index(n, target_fraction) is None \
        or find_farey_index(n, target_fraction) == find_farey_index(n, target_fraction)
    
    f = Farey(n)
    assert f.n == n

    f = Farey(n)
    assert f.farey_list == farey_list(n)

    f = Farey(n)
    assert find_farey_index(n, target_fraction) is None and target_fraction not in f \
        or find_farey_index(n, target_fraction) is not None and target_fraction in f

    f = Farey(n)
    assert f.farey_list == [farey_fraction for farey_fraction in f]

    f = Farey(n)
    assert f.mertens_function == mertens_function(n)

    n = 23
    f = Farey(n)
    assert f.mertens_function == prime_numbers.mertens(n)




