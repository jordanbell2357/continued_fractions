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


EULER_CONSTANT_30_DIGITS = Decimal("0.577215664901532860606512090082")
EULER_CONSTANT_FLOAT = float(EULER_CONSTANT_30_DIGITS)


def sieve_eratosthenes_list(n: int) -> list[int]:
    if n < 2:
        return []
    sieve_list = [True] * (n + 1)
    sieve_list[0] = False
    sieve_list[1] = False
    for q in range(2, math.isqrt(n) + 1):
        if sieve_list[q] == 1:
            sieve_list[q * q::q] = [False] * ((n - q * q) // q + 1)
    return list(it.compress(range(n + 1), sieve_list))
    # # same as
    # return [k for k in range(n + 1) if sieve_list[k] == True]


def sieve_eratosthenes(n: int) -> list[int]:
    """Indices of bytearray are 0 to n. We set the bytes at each index initially to 1.
    For each index whose byte is still set to 1 when "encountered", we set the bytes at all
    multiples of the index to 0. Termination happens when there are no more indices to "encounter".
    At termination, the indices whose bytes are still equal to 1, are the prime numbers between 0 and n.
    """

    if n < 2:
        return []
    sieve_bytearray = bytearray(b"\x01" * (n + 1))
    sieve_bytearray[0:2] = b"\x00\x00" # 0 and 1 are not prime
    for q in range(2, math.isqrt(n) + 1):
        if sieve_bytearray[q] == 1:
            sieve_bytearray[q * q::q] = b"\x00" * ((n - q * q) // q + 1) # set bytes at multiples of index to 0
    return list(it.compress(range(n + 1), sieve_bytearray)) # return those indices whose byte is 1


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
    

def characteristic_function_primes(n: int) -> int:
    return int(isprime(n))


def generate_primes() -> abc.Generator[int]:
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


def mobius(n: int) -> int:
    complex_sum = sum(cmath.exp(complex(0, 2 * math.pi * Fraction(k, n))) for k in range(n) if math.gcd(n, k) == 1)
    return round(complex_sum.real)


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
    if n == 0:
        return 0
    prime_factor_counter = make_prime_factor_counter(n)
    divisor_count = math.prod(a + 1 for _, a in prime_factor_counter.items())
    return divisor_count


def sum_of_divisors(n: int) -> int:
    if n == 0:
        return 0
    prime_factor_counter = make_prime_factor_counter(n)
    divisor_sum = math.prod(Fraction(p ** (a + 1) - 1, p - 1) for p, a in prime_factor_counter.items())
    return int(divisor_sum)


def divisor_summatory_function(x: float) -> int:
    return sum(number_of_divisors(n) for n in range(1, int(x) + 1))


def divisor_summatory_function_hyperbola(x: float) -> int:
    return sum(math.floor(x / n) for n in range(1, int(x) + 1))


def bernoulli(n: int) -> Fraction:
    return sum(Fraction(1, k + 1) * sum(math.comb(k, j) * (-1)**j * j**n for j in range(k + 1)) for k in range(n + 1))


def faulhaber(m: int, n: int) -> int:
    return Fraction(1, m + 1) * sum((-1)**k * math.comb(m + 1, k) * bernoulli(k) * Fraction(pow(n, m + 1), pow(n, k)) for k in range(m + 1))


def euler_totient(n: int) -> int:
    prime_factor_list = make_prime_factor_list(n)
    return n * math.prod(Fraction(1, 1) - Fraction(1, p) for p in prime_factor_list)


def dirichlet_convolution(f: abc.Callable, g: abc.Callable, n: int) -> int:
    factor_list = make_factor_list(n)
    return sum(f(d) * g(n // d) for d in factor_list)


def harmonic_sum_fraction(n: int) -> Fraction:
    return sum(Fraction(1, k) for k in range(1, n + 1))


def harmonic_sum_decimal(n: int, precision: int) -> Decimal:
    decimal.getcontext().prec = precision
    return sum(Decimal(1) / Decimal(k) for k in range(1, n + 1))


def euler_constant(n: int, m: int, precision: int, method: str = "decimal") -> Decimal:
    """Uses Euler-Maclaurin summation formula."""
    decimal.getcontext().prec = precision
    if method == "decimal":
        harmonic_sum_decimal_n = harmonic_sum_decimal(n, precision)
    elif method == "fraction":
        harmonic_sum_fraction_n = harmonic_sum_fraction(n)
        harmonic_sum_decimal_n = Decimal(harmonic_sum_fraction_n.numerator) / Decimal(harmonic_sum_fraction_n.denominator)
    euler_maclaurin_sum_fraction_m = sum(bernoulli(2 * k) / (2 * k * n**(2 * k)) for k in range(1, m + 1))
    euler_maclaurin_sum_decimal_m = Decimal(euler_maclaurin_sum_fraction_m.numerator) / Decimal(euler_maclaurin_sum_fraction_m.denominator)
    return harmonic_sum_decimal_n - Decimal(n).ln() - 1 / Decimal(2 * n) + euler_maclaurin_sum_decimal_m


if __name__ == "__main__":
    # Example: prime numbers
    n = 200

    assert sieve_eratosthenes_list(n) == sieve_eratosthenes(n)

    assert all(isprime(p) for p in sieve_eratosthenes(n))

    assert len(sieve_eratosthenes(n)) == prime_counting_function(n)

    assert len(make_prime_list(n)) == n

    prime_factor_counter = make_prime_factor_counter(n)
    assert all(valuation(p, n) == v for p, v in prime_factor_counter.items())

    prime_factor_counter = make_prime_factor_counter(n)
    assert math.prod(p ** a for p, a in prime_factor_counter.items()) == n

    assert len(make_factor_list(n)) == number_of_divisors(n)

    assert sum(make_factor_list(n)) == sum_of_divisors(n)

    # sum of empty list is 0 (list of integers d satisfying 0 < d <= 0)
    assert number_of_divisors(0) == 0


    # Example: Legendre's formula for v_p(n!)
    n = 150
    p = 3
    assert valuation(p, math.factorial(n)) == sum(int(n / p ** i) for i in range(1, int(math.log(n, p) + 1)))



    # Example: arithmetic functions
    m = 15
    n = 107
    assert math.gcd(m, n) > 1 or prime_little_omega(m * n) == prime_little_omega(m) + prime_little_omega(n)

    assert prime_big_omega(m * n) == prime_big_omega(m) + prime_big_omega(n)

    assert prime_big_omega(n) == prime_little_omega(n) or mobius(n) == 0

    assert mobius(n) == mobius_prime(n)

    assert math.gcd(m, n) > 1 or mobius(m * n) == mobius(m) * mobius(n)

    assert prime_big_omega(n) == prime_little_omega(n) and mobius(n) == liouville(n) or \
        prime_big_omega(n) > prime_little_omega(n) and mobius(n) == 0
    
    assert liouville(n) == (-1) ** prime_big_omega(n)

    assert dirichlet_convolution(prime_little_omega, mobius, n) == characteristic_function_primes(n)


    # Example: Bernoulli numbers
    n = 20
    m = 4
    assert faulhaber(m, n) == sum(k**m for k in range(n + 1))


    # Example: Euler's constant
    n = 100
    m = 8
    precision = 20
    assert math.isclose(float(euler_constant(n, m, precision)), EULER_CONSTANT_FLOAT)


    # Mertens' third theorem
    N = 200_000
    assert math.isclose(round(1 / math.log(N) * math.prod(p / (p - 1) for p in sieve_eratosthenes(N)), 3), round(math.exp(EULER_CONSTANT_FLOAT), 3))


    # Example: divisor summatory function
    n = 150
    C = 2

    assert divisor_summatory_function_hyperbola(n) == divisor_summatory_function(n)

    assert divisor_summatory_function(n) <= n * math.log(n) + n * (2 * EULER_CONSTANT_FLOAT - 1) + C * math.isqrt(n)


    # Example: Dirichlet convolution
    n = 620

    assert all(dirichlet_convolution(euler_totient, lambda _: 1, k) == k for k in range(1, n + 1))

    assert all(dirichlet_convolution(number_of_divisors, mobius, k) == 1 for k in range(1, n + 1))

    assert all(dirichlet_convolution(sum_of_divisors, mobius, k) == k for k in range(1, n + 1))

    assert all(dirichlet_convolution(lambda _: 1, lambda _: 1, k) == number_of_divisors(k) for k in range(1, n + 1))

    assert all(dirichlet_convolution(lambda x: x, lambda _: 1, k) == sum_of_divisors(k) for k in range(1, n + 1))