import math
import cmath
from fractions import Fraction
import itertools as it
from collections import Counter


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

def make_prime_list(n: int) -> list[int]:
    prime_list = []
    if n > len(prime_list):
        prime_list.append(2)
    for k in it.count(start=3, step=2):
        if len(prime_list) >= n:
            break
        if isprime(k):
            prime_list.append(k)
    return prime_list

def primes():
    yield 2
    for k in it.count(start=3, step=2):
        if isprime(k):
            yield k

def factor(n: int) -> Counter:
    prime_counter = Counter()
    prime_gen = primes()
    for p in prime_gen:
        if n == 1:
            break
        while n % p == 0:
            prime_counter.update([p])
            n = n // p
    return prime_counter

def mu_complex(n: int) -> complex:
    return sum(cmath.exp(2j * math.pi * Fraction(k, n)) for k in range(n) if math.gcd(n, k) == 1)

def mu(n: int) -> int:
    return int(mu_complex(n).real)


