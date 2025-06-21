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


def sieve_eratosthenes(n: int) -> list[int]:
    """
    Indices of bytearray are 0 to n. We set the bytes at each index initially to 1.
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


def squarefull_and_squarefree_parts(d: int) -> tuple[int, int]:
        prime_factor_counter_d = make_prime_factor_counter(d)
        squarefull_part_d = 1
        squarefree_part_d = 1
        for p, k in prime_factor_counter_d.items():
            squarefull_part_d *= p ** (k // 2 * 2)
            squarefree_part_d *= p ** (k % 2)
        return squarefull_part_d, squarefree_part_d


def is_squarefree(d: int) -> bool:
    _, squarefree_part_d = squarefull_and_squarefree_parts(d)
    return squarefree_part_d == d


def legendre_symbol(a: int, p: int) -> int:
    """
    We extend the usual Legendre symbol to include p=2,
    following the extension by the Kronecker symbol.
    """
    if not isprime(p):
        raise ValueError(f"p must be prime: {p=}")
    if p == 2:
        r = a % 8
        if r in [0, 2, 4, 6]:
            return 0
        elif r in [1, 7]:
            return 1
        elif r in [3, 5]:
            return -1
    if a % p == 0:
        return 0
    if any((a - x ** 2) % p == 0 for x in range(1, p)):
        return 1
    else:
        return -1


def kronecker_symbol(a: int, n: int) -> int:
    """
    Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
    Definition 1.4.8, p. 28, 
    """
    if n == 0:
        return 1 if a in [-1, 1] else 0
    elif n == -1:
        return -1 if a < 0 else 1
    elif n == 1:
        return 1

    u = n // abs(n)
    n = n // u
    prime_factor_counter_n  = make_prime_factor_counter(n)
    legendre_symbol_product_n = math.prod(legendre_symbol(a, p) ** v for p, v in prime_factor_counter_n.items())
    return kronecker_symbol(a, u) * legendre_symbol_product_n


def solve_quadratic_congruence(a: int, n: int) -> int:
    for b in range(n):
        if (b ** 2 - a) % n == 0:
            return b
    return None


def square_root_mod_p(d: int, p: int) -> int:
    """Return t with t² ≡ d (mod p) when it exists (Tonelli–Shanks)."""
    if p == 2:
        return d % 2                                 # trivial root
    # simple Tonelli–Shanks; assumes (d/p)=+1
    leg = pow(d, (p - 1) // 2, p)
    if leg != 1:
        raise ValueError(f"{d} is not a quadratic residue mod {p}.")
    # find q·2^s with q odd
    s, q = 0, p - 1
    while q % 2 == 0:
        s += 1
        q //= 2
    # find z a non-square
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(d, q, p), pow(d, (q + 1) // 2, p)
    while t != 1:
        i, tmp = 1, pow(t, 2, p)
        while tmp != 1:
            tmp = pow(tmp, 2, p)
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, pow(b, 2, p), t * pow(b, 2, p) % p, r * b % p
    return r



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


    # Example: numer of divisors
    n = 150
    C = 2

    assert sum(number_of_divisors(k) for k in range(1, n + 1)) == sum(n // k for k in range(1, n + 1))

    assert sum(number_of_divisors(k) for k in range(1, n + 1)) <= n * math.log(n) + n * (2 * EULER_CONSTANT_FLOAT - 1) + C * math.isqrt(n)


    # Example: Dirichlet convolution
    n = 620

    assert all(dirichlet_convolution(euler_totient, lambda _: 1, k) == k for k in range(1, n + 1))

    assert all(dirichlet_convolution(number_of_divisors, mobius, k) == 1 for k in range(1, n + 1))

    assert all(dirichlet_convolution(sum_of_divisors, mobius, k) == k for k in range(1, n + 1))

    assert all(dirichlet_convolution(lambda _: 1, lambda _: 1, k) == number_of_divisors(k) for k in range(1, n + 1))

    assert all(dirichlet_convolution(lambda x: x, lambda _: 1, k) == sum_of_divisors(k) for k in range(1, n + 1))