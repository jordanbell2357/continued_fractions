import math
import cmath
from fractions import Fraction
import itertools as it
import functools as ft
import operator
from collections import Counter
import decimal
from decimal import Decimal

import farey


def sieve_eratosthenes(x: float) -> list[int]:
    n = int(x)

    if n < 2:
        return []

    # one byte per number: 1 = prime candidate, 0 = composite
    flags = bytearray(b"\x01") * (n + 1)
    flags[0:2] = b"\x00\x00"                     # 0 and 1 are not prime

    for p in range(2, math.isqrt(n) + 1):        # stop at ⌊√n⌋
        if flags[p]:                             # p is still marked prime
            flags[p * p::p] = b"\x00" * ((n - p*p)//p + 1)
            # slice-assignment wipes out multiples

    # itertools.compress keeps numbers whose flag byte is still 1
    return list(it.compress(range(n + 1), flags))


def isprime(n: int) -> bool:
    n = abs(n)
    if n % 2 == 0:
        if n == 2:
            return True
        else:
            return False
    else:
        for k in range(3, math.isqrt(n) + 1, 2):
            if n % k == 0:
                return False
        return True


def generate_primes():
    yield 2
    for k in it.count(start=3, step=2):
        if isprime(k):
            yield k


def make_prime_list(n: int) -> list[int]:
    prime_gen = generate_primes()
    prime_list = [p for p, _ in zip(prime_gen, range(n))]
    return prime_list


def prime_counting_function(x: float) -> int:
    prime_gen = generate_primes()
    prime_count = 0
    for p in prime_gen:
        if p > x:
            break
        else:
            prime_count += 1
    return prime_count


def prime_power_counting_function(x: float, precision: int = 20) -> Decimal:
    decimal_context = decimal.getcontext()
    decimal_context.prec = precision
    endpoint = math.floor(math.log2(x)) + 1
    return sum(Decimal(prime_counting_function(x ** Fraction(1, n))) / Decimal(n) for n in range(1, endpoint))


def make_prime_factor_counter(n: int) -> Counter:
    prime_gen = generate_primes()
    prime_factor_counter = Counter()
    for p in prime_gen:
        if n in [0, 1]:
            break
        while n % p == 0:
            prime_factor_counter.update([p])
            n = n // p
    return prime_factor_counter


def valuation(p: int, n: int) -> int:
    v = 0
    while n % p == 0:
        v += 1
        n = n // p
    return v


def von_mangoldt(n: int, precision: int = 20) -> Decimal:
    decimal.getcontext().prec = precision
    prime_factor_counter = make_prime_factor_counter(n)
    if len(prime_factor_counter) != 1:
        return Decimal(0)
    else:
        p, _ = prime_factor_counter.popitem()
        decimal_log_p = Decimal(p).ln()
        return decimal_log_p


def chebyshev_theta(x: float, precision: int = 20) -> Decimal:
    decimal.getcontext().prec = precision
    log_prime_sum = Decimal(0)
    prime_gen = generate_primes()
    for p in prime_gen:
        if p <= x:
            log_prime_sum += Decimal(p).ln()
        else:
            break
    return log_prime_sum


def chebyshev_psi(x: float, precision: int = 20) -> Decimal:
    decimal.getcontext().prec = precision
    von_mangoldt_sum = sum([von_mangoldt(n, precision) for n in range(math.floor(x) + 1)])
    return von_mangoldt_sum


def mobius_complex(n: int) -> complex:
    return sum(cmath.exp(complex(0, 2 * math.pi * Fraction(k, n))) for k in range(n) if math.gcd(n, k) == 1)


def mobius(n: int) -> int:
    return round(mobius_complex(n).real)


def mobius_prime(n: int) -> int:
    prime_factor_counter = make_prime_factor_counter(n)
    if any(v > 1 for v in prime_factor_counter.values()):
        return 0
    else:
        return (-1) ** len(prime_factor_counter)


def mertens(n: int) -> int:
    return sum(mobius(k) for k in range(1, n + 1))


def prime_little_omega(n: int) -> int:
    prime_factor_counter = make_prime_factor_counter(n)
    return len(prime_factor_counter)


def prime_big_omega(n: int) -> int:
    prime_factor_counter = make_prime_factor_counter(n)
    return sum(prime_factor_counter.values())


def liouville(n: int) -> int:
    return (-1) ** prime_big_omega(n)


def make_factor_list(n: int) -> int:
    prime_factor_counter = make_prime_factor_counter(n)
    exponent_range_list = [range(prime_factor_counter[p] + 1) for p in prime_factor_counter]
    factor_list = [ft.reduce(operator.mul, [p ** e for p, e in zip(prime_factor_counter.keys(), exponent_list)], 1)
                   for exponent_list in it.product(*exponent_range_list)]
    return sorted(factor_list)


def number_of_divisors(n: int) -> int:
    prime_factor_counter = make_prime_factor_counter(n)
    divisor_count = math.prod(a + 1 for _, a in prime_factor_counter.items())
    return divisor_count


def sum_of_divisors(n: int) -> int:
    prime_factor_counter = make_prime_factor_counter(n)
    divisor_sum = math.prod(Fraction(p ** (a + 1) - 1, p - 1) for p, a in prime_factor_counter.items())
    return int(divisor_sum)


if __name__ == "__main__":
    n = 14
    x = 180.45
    precision = 40
    display_precision = 20

    assert all(isprime(p) for p in sieve_eratosthenes(n))

    assert len(sieve_eratosthenes(x)) == prime_counting_function(x)

    assert len(make_prime_list(n)) == n

    prime_factor_counter = make_prime_factor_counter(n)
    assert all(valuation(p, n) == v for p, v in prime_factor_counter.items())

    prime_factor_counter = make_prime_factor_counter(n)
    assert math.prod(p ** a for p, a in prime_factor_counter.items()) == n

    assert len(make_factor_list(n)) == number_of_divisors(n)

    assert sum(make_factor_list(n)) == sum_of_divisors(n)

    assert mobius(n) == mobius_prime(n)

    assert prime_big_omega(n) == prime_little_omega(n) or mobius(n) == 0

    assert liouville(n) == (-1) ** prime_big_omega(n)

    assert mobius(n) == liouville(n) if prime_big_omega(n) == prime_little_omega(n) else 0

    assert decimal.Context(prec=display_precision).create_decimal(chebyshev_psi(x, precision)) == \
        decimal.Context(prec=display_precision).create_decimal(sum(chebyshev_theta(x ** Fraction(1, n), precision) for n in range(1, math.floor(math.log2(x)) + 1)))

    assert decimal.Context(prec=display_precision).create_decimal(prime_power_counting_function(x, precision)) == \
        decimal.Context(prec=display_precision).create_decimal(sum(von_mangoldt(n, precision) / Decimal(n).ln() for n in range(2, math.floor(x) + 1)))

    f = farey.Farey(n)
    assert f.mertens_function == mertens(n)

    assert farey.mertens_function(n) == mertens(n)

