import math
from fractions import Fraction
import decimal
from decimal import Decimal
import itertools as it
import functools as ft
from collections import abc

import cflib


def sqrt_periodic_cf(d: int) -> tuple[list[int], list[int]]:
    a0 = math.isqrt(d)  # floor(sqrt(d))
    if a0 ** 2 == d:
        # d is a square, so sqrt(d) is an integer and the continued fraction has length 1
        return ([a0], [])
    # Initialize (M, D, a)
    M = 0
    D = 1
    a = a0
    expansion = []
    expansion.append(a0)
    # Keep track of triples (M, D, a) to detect repetition
    triple_dict = {}
    # Map: (M, D, a) -> index in the expansion
    i = 0
    while True:
        M = D * a - M
        D = (d - M**2) // D
        a = (a0 + M) // D
        i += 1
        # If we see a triple (M, D, a) twice, we've found the start of the period
        if (M, D, a) in triple_dict:
            # The terms from index triple_dict[(M,D,a)] to index i-1 are the repeating cycle
            start_period = triple_dict[(M, D, a)]
            initial_part = expansion[:start_period]
            periodic_part = expansion[start_period:]
            return initial_part, periodic_part
        else:
            triple_dict[(M, D, a)] = i
            expansion.append(a)


def integral_quadratic_polynomial_to_periodic_cf(
        a: int,
        b: int,
        c: int,
        epsilon: int = 1) -> tuple[list[int], list[int]]:
    """
    Return (initial_part, periodic_part) for the quadratic irrational

            (-b ± √D) / (2a),   where  D = b² - 4ac  ≥ 0  and  a ≠ 0.

    epsilon:
        1  →  use  (-b + √D)/(2a)   (like sqrt_periodic_cf)
       -1  →  use  (-b - √D)/(2a)

    The algorithm is the same three-term recurrence coded for sqrt_periodic_cf:
    only the *starting triple* (P, Q, a₀) differs.

    If D is a perfect square the root is rational so the continued fraction
    terminates after its first term, we return  ([q0], []) just like the
    sqrt(d) routine.

    Implements construction in proof of Theorem 177, Hardy and Wright, Chapter X, pp. 185-187, sixth edition:
    The continued fraction which represents a quadratic surd is periodic.
    """
    if a == 0:
        raise ValueError("Not a quadratic polynomial: a must be non-zero")

    D = b * b - 4 * a * c
    if D < 0:
        raise ValueError("Discriminant is negative: roots are complex")

    isqrt_D = math.isqrt(D)
    if isqrt_D * isqrt_D == D: # Rational root implies finite continued fraction
        numerator = -b + epsilon * isqrt_D
        denominator = 2 * a
        if denominator == 0:
            raise ZeroDivisionError("Denominator 2a vanished (a=0)")
        if denominator < 0: # Make denominator positive
            numerator, denominator = -numerator, -denominator
        q0, r = divmod(numerator, denominator)
        if r != 0: # safety check; should never occur
            q0 = math.floor(numerator / denominator)
        return ([q0], [])

    # We want the reduced form (P + √D) / Q with Q > 0 and Q | D - P²
    P = -b
    Q = 2 * a
    P = epsilon * P
    if Q < 0:
        P, Q = -P, -Q

    gcd = math.gcd(math.gcd(abs(P), Q), D)
    P //= gcd
    Q //= gcd
    D //= gcd
    isqrt_D = math.isqrt(D)

    a0 = (P + isqrt_D) // Q
    expansion = []
    expansion.append(a0)

    #  We only need to record (P, Q) to recognise when the state repeats
    state_dict = {(P, Q): 0}

    while True:
        P = a0 * Q - P
        Q = (D - P * P) // Q
        a0 = (P + isqrt_D) // Q

        state = (P, Q)
        if state in state_dict: # cycle found
            start = state_dict[state]
            return expansion[:start], expansion[start:]

        state_dict[state] = len(expansion)
        expansion.append(a0) 


def periodic_cf_generator(d: int) -> abc.Generator[int]:
    initial_part, periodic_part = sqrt_periodic_cf(d)
    cf = it.chain(initial_part, it.cycle(periodic_part))
    yield from cf


def solve_pell_equation(d: int) -> tuple[int, int]:
    initial_part, periodic_part = sqrt_periodic_cf(d)
    r = len(periodic_part)
    if r % 2 == 0:
        partial_quotients = initial_part + periodic_part
        fraction_list = cflib.cf_to_convergent_list(partial_quotients)
        x, y = fraction_list[r - 1]
    else:
        partial_quotients = initial_part + 2 * periodic_part
        fraction_list = cflib.cf_to_convergent_list(partial_quotients)
        x, y = fraction_list[2 * r - 1]
    return x, y


def periodic_cf_to_integral_quadratic_polynomial(
        initial_part: list[int],
        periodic_part: list[int]) -> tuple[int, int, int]:
    """Implements construction in proof of Theorem 176 of Hardy and Wright, Chapter X, pp. 184-185, sixth edition:
    Theorem 176. A periodic continued fraction is a quadratic surd,
    i.e. an irrational root of a quadratic equation with integral coefficients."""
    initial_convergent_list = cflib.cf_to_extended_convergent_list(initial_part)
    initial_p1, initial_q1 = initial_convergent_list[-2]
    initial_p2, initial_q2 = initial_convergent_list[-1]

    periodic_convergent_list = cflib.cf_to_extended_convergent_list(periodic_part)
    periodic_p1, periodic_q1 = periodic_convergent_list[-2]
    periodic_p2, periodic_q2 = periodic_convergent_list[-1]

    a = periodic_q2 * initial_q1 ** 2 - (periodic_q1 - periodic_p2) * initial_q1 * initial_q2 - periodic_p1 * initial_q2 ** 2
    b = -2 * periodic_q2 * initial_q1 * initial_p1 + (periodic_q1 - periodic_p2) * (initial_p1 * initial_q2 + initial_q1 * initial_p2) + \
        2 * periodic_p1 * initial_p2 * initial_q2
    c = periodic_q2 * initial_p1 ** 2 - (periodic_q1 - periodic_p2) * initial_p1 * initial_p2 - periodic_p1 * initial_p2 ** 2

    gcd = ft.reduce(math.gcd, [b, c], a)
    return a // gcd, b // gcd, c // gcd


def is_reduced_surd(a, b, c):
    D = b ** 2 - 4 * a * c
    P = -b
    Q = 2 * a
    r = (P + math.sqrt(D)) / Q
    r_conjugate = (P - math.sqrt(D)) / Q
    if r > 1 and -1 < r_conjugate < 1:
        return True
    return False


# Example usage:
if __name__ == "__main__":
    d = 5
    num_terms = 30
    precision = 30

    print(integral_quadratic_polynomial_to_periodic_cf(1, -1, -1))

    decimal.getcontext().prec = precision

    initial_part, periodic_part = sqrt_periodic_cf(d)
    print(f"Continued fraction for sqrt({d})")
    print(f"Initial part:\t{initial_part}")
    print(f"Periodic part:\t{periodic_part}")

    print()
    
    print(f"Pell equation x² - {d}y² = 1 has following fundamental solution:")
    print("(x, y)=", solve_pell_equation(d), sep="")

    print()

    print(f"Using {num_terms=} of continued fraction for sqrt({d}), the decimal expansion with {precision=} is:")
    cf_generator = periodic_cf_generator(d)
    cf = [a for a, _ in zip(cf_generator, range(num_terms))]
    convergent_list = cflib.cf_to_convergent_list(cf)
    fraction_tuple = convergent_list[-1]
    print(Decimal(fraction_tuple[0]) / Decimal(fraction_tuple[1]))

    print()

    initial_part, periodic_part = sqrt_periodic_cf(d)
    a, b, c = periodic_cf_to_integral_quadratic_polynomial(initial_part, periodic_part)
    print(f"""sqrt({d}) is a root of the following quadratic polynomial:\n{a}x² \t{b:+}x \t{c:+}""")

    print()

    initial_part, periodic_part = [], [1]
    a, b, c = periodic_cf_to_integral_quadratic_polynomial(initial_part, periodic_part)
    print(f"The value of the continued fraction is a root of the following quadratic polynomial:\n{a}x² \t{b:+}x \t{c:+}")

    print()

    # Galois' theorem on quadratic surds: a quadratic surd is reduced if and only if it has a purely periodic cf
    # A reduced surd (P + sqrt(D)) / Q is one that is greater than 1 and whose conjugate (P - sqrt(D)) / Q is
    # strictly between -1 and 0.
    print(is_reduced_surd(1, -1, -1))

