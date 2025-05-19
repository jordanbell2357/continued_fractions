import math
import cmath
from fractions import Fraction
import itertools as it
import functools as ft
import operator
import decimal
from decimal import Decimal
from collections import Counter
from collections import abc
import typing


def sieve_eratosthenes(n: int) -> list[int]:
    """Indices of bytearray are 0 to n. We set the bytes at each index initially to 1. For each index
    whose byte is still set to 1, we set the bytes at all multiples of the index to 0.
    The indices whose bytes are still equal to 1 at termination, are the prime numbers between
    0 and n.
    """

    if n < 2:
        return []
    sieve_bytearray = bytearray(b"\x01" * (n + 1))
    sieve_bytearray[0:2] = b"\x00\x00" # 0 and 1 are not prime
    for q in range(2, math.isqrt(n) + 1):
        if sieve_bytearray[q] == 1:
            sieve_bytearray[q * q::q] = b"\x00" * ((n - q * q) // q + 1) # set bytes at multiples of index to 0
    return list(it.compress(range(n + 1), sieve_bytearray)) # return those indices whose byte is 1


def sieve_eratosthenes_list(n: int) -> list[int]:
    if n < 2:
        return []
    sieve_list = [True] * (n + 1)
    sieve_list[0] = False
    sieve_list[1] = False
    for q in range(2, math.isqrt(n) + 1):
        if sieve_list[q] == 1:
            sieve_list[q * q::q] = [False] * ((n - q * q) // q + 1)
    return [k for k in range(n + 1) if sieve_list[k] == True]


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


def make_factor_list(n: int) -> list[int]:
    prime_factor_counter = make_prime_factor_counter(n)
    exponent_range_list = [range(prime_factor_counter[p] + 1) for p in prime_factor_counter]
    factor_list = [ft.reduce(operator.mul, [p ** e for p, e in zip(prime_factor_counter.keys(), exponent_list)], 1)
                   for exponent_list in it.product(*exponent_range_list)]
    return sorted(factor_list)


def make_prime_factor_list(n: int) -> list[int]:
    prime_factor_counter = make_prime_factor_counter(n)
    return prime_factor_counter.keys()


def number_of_divisors(n: int) -> int:
    prime_factor_counter = make_prime_factor_counter(n)
    divisor_count = math.prod(a + 1 for _, a in prime_factor_counter.items())
    return divisor_count


def sum_of_divisors(n: int) -> int:
    prime_factor_counter = make_prime_factor_counter(n)
    divisor_sum = math.prod(Fraction(p ** (a + 1) - 1, p - 1) for p, a in prime_factor_counter.items())
    return int(divisor_sum)


def bernoulli(n: int) -> Fraction:
    return sum(Fraction(1, k + 1) * sum(math.comb(k, j) * (-1)**j * j**n for j in range(k + 1)) for k in range(n + 1))


def faulhaber(m: int, n: int) -> int:
    return Fraction(1, m + 1) * sum((-1)**k * math.comb(m + 1, k) * bernoulli(k) * Fraction(pow(n, m + 1), pow(n, k)) for k in range(m + 1))


def ramanujan_sum_c(q: int, n: int) -> int:
    return round(sum(cmath.exp(complex(0, 2 * math.pi * Fraction(a, q) * n)) for a in range(1, q + 1) if math.gcd(q, a) == 1).real)


def dedekind_psi(n: int) -> int:
    prime_factor_list = make_prime_factor_list(n)
    return n * math.prod(Fraction(1, 1) + Fraction(1, p) for p in prime_factor_list)


def dirichlet_convolution(f: abc.Callable, g: abc.Callable, n: int) -> complex:
    factor_list = make_factor_list(n)
    return sum(f(d) * g(n // d) for d in factor_list)


def dirichlet_convolution_identity(n: int) -> int:
    return 1 if n == 1 else 0


if __name__ == "__main__":
    n = 4
    x = 180.45
    precision = 40
    display_precision = 20
    m = 3
    q = 14

    assert all(dirichlet_convolution(lambda _: 1, mobius, k) == dirichlet_convolution_identity(k) for k in range(1, n + 1))

    assert all(dirichlet_convolution(number_of_divisors, mobius, k) == 1 for k in range(1, n + 1))

    assert all(dirichlet_convolution(sum_of_divisors, mobius, k) == k for k in range(1, n + 1))

    assert all(dirichlet_convolution(lambda _: 1, lambda _: 1, k) == number_of_divisors(k) for k in range(1, n + 1))

    assert all(dirichlet_convolution(lambda x: x, lambda _: 1, k) == sum_of_divisors(k) for k in range(1, n + 1))

    assert sieve_eratosthenes_list(n) == sieve_eratosthenes(n)

    assert all(isprime(p) for p in sieve_eratosthenes(n))

    assert len(sieve_eratosthenes(int(x))) == prime_counting_function(x)

    assert ramanujan_sum_c(q, 1) == mobius(q)

    factor_list = make_factor_list(math.gcd(n, q))
    assert ramanujan_sum_c(q, n) == sum(mobius(q // d) * d for d in factor_list)

    assert faulhaber(m, n) == sum(k**m for k in range(n + 1))

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

    assert prime_big_omega(n) == prime_little_omega(n) and mobius(n) == liouville(n) or \
        prime_big_omega(n) > prime_little_omega(n) and mobius(n) == 0
    
    assert dirichlet_convolution(von_mangoldt, lambda _: 1, Decimal(n)) == Decimal(n).ln()

    assert decimal.Context(prec=display_precision).create_decimal(chebyshev_psi(x, precision)) == \
        decimal.Context(prec=display_precision).create_decimal(sum(chebyshev_theta(x ** Fraction(1, n), precision)
                                                                   for n in range(1, math.floor(math.log2(x)) + 1)))

    assert decimal.Context(prec=display_precision).create_decimal(prime_power_counting_function(x, precision)) == \
        decimal.Context(prec=display_precision).create_decimal(sum(von_mangoldt(n, precision) / Decimal(n).ln()
                                                                   for n in range(2, math.floor(x) + 1)))
    

