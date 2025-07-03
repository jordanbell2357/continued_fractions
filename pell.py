import math
import decimal
from decimal import Decimal
from fractions import Fraction
import itertools as it
from collections import abc
from dataclasses import dataclass

import cflib


@dataclass
class PeriodicCF:
    initial_part: list[int]
    periodic_part: list[int]


def sqrt_periodic_cf(d: int) -> PeriodicCF:
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
            # The terms from index triple_dict[(M, D, a)] to index i - 1 are the repeating cycle
            start_period = triple_dict[(M, D, a)]
            initial_part = expansion[:start_period]
            periodic_part = expansion[start_period:]
            return PeriodicCF(initial_part, periodic_part)
        else:
            triple_dict[(M, D, a)] = i
            expansion.append(a)


def integral_quadratic_polynomial_to_periodic_cf(
        a: int,
        b: int,
        c: int,
        epsilon: int = 1) -> PeriodicCF:
    """
    Return PeriodicCF(initial_part, periodic_part) for the quadratic irrational

            (-b ± √D) / (2a),   where  D = b² - 4ac  ≥ 0  and  a ≠ 0.

    epsilon:
        +1 → use  (-b + √D)/(2a)
        -1 → use  (-b - √D)/(2a)

    Implements construction in Theorem 177 of Hardy and Wright, Chapter X, pp. 185-187, sixth edition:
    "The continued fraction which represents a quadratic surd is periodic."

    Chakravala method used by Bhāskara II. Theorem of Lagrange.

    Manfred Einsiedler and Thomas Ward, Ergodic Theory with a view towards Number Theory,
    Graduate Texts in Mathematics, Volume 259, Springer, 2011.
    Theorem 3.13 (Lagrange). Let u be an irrational positive real number.
    Then the continued fraction expansion of u is eventually periodic if and only
    if u is a quadratic irrational.
    """
    if a == 0:
        raise ValueError("Not a quadratic polynomial: a must be non-zero")

    D0 = b * b - 4 * a * c
    if D0 < 0:
        raise ValueError("Discriminant is negative: roots are complex")

    # rational case
    isqrt_D0 = math.isqrt(D0)
    if isqrt_D0 * isqrt_D0 == D0:
        numerator   = -b + epsilon * isqrt_D0
        denominator = 2 * a
        if denominator < 0:
            numerator, denominator = -numerator, -denominator
        q0 = numerator // denominator
        return PeriodicCF([q0], [])

    # --- reduction so that Q | (D - P^2) remains true ---
    P0 = epsilon * (-b)
    Q0 = 2 * a
    if Q0 < 0:
        P0, Q0 = -P0, -Q0

    g = math.gcd(abs(P0), Q0)
    P = P0 // g
    Q = Q0 // g
    D = D0 // (g * g)
    isqrt_D = math.isqrt(D)
    # --- end reduction ---

    i = 0 # first partial quotient
    partial_quotient = (P + isqrt_D) // Q
    expansion = [partial_quotient]
    state_dict = {}
    state_dict[(P, Q, partial_quotient)] = 0

    while True:
        # recurrence step
        P   = partial_quotient * Q - P
        Q   = (D - P * P) // Q
        partial_quotient  = (P + isqrt_D) // Q
        i  += 1

        state = (P, Q, partial_quotient)
        if state in state_dict:
            start = state_dict[state]
            initial_part = expansion[:start]
            periodic_part = expansion[start:i]
            return PeriodicCF(initial_part, periodic_part)

        state_dict[state] = i
        expansion.append(partial_quotient)


def periodic_cf_generator(periodic_cf: PeriodicCF) -> abc.Generator[int]:
    initial_part, periodic_part = periodic_cf.initial_part, periodic_cf.periodic_part
    cf = it.chain(initial_part, it.cycle(periodic_part))
    yield from cf


def solve_pell_equation(d: int) -> tuple[int, int]:
    """
    x^2 - dy^2 = 1
    """
    periodic_cf = sqrt_periodic_cf(d)
    initial_part, periodic_part = periodic_cf.initial_part, periodic_cf.periodic_part
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


def periodic_cf_to_integral_quadratic_polynomial(periodic_cf: PeriodicCF) -> tuple[int, int, int]:
    """
    Implements construction in proof of Theorem 176 of Hardy and Wright, Chapter X, pp. 184-185, sixth edition:
    "A periodic continued fraction is a quadratic surd, i.e. an irrational root of a quadratic equation with integral coefficients."
    Theorem of Euler.
    """

    initial_part, periodic_part = periodic_cf.initial_part, periodic_cf.periodic_part
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

    gcd = math.gcd(a, b, c)
    return a // gcd, b // gcd, c // gcd


def is_real_reduced_surd(P: int, Q: int, D: int) -> bool:
    """
    Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
    Corollary 5.6.7, p. 266:
    A quadratic surd (P+√D)/Q is called "reduced" when it is greater than 1,
    and its conjugate (P-√D)/Q is strictly between -1 and 0.
    
    Galois' theorem for quadratic surds: a quadratic surd is reduced if and only if its
    continued fraction is purely periodic.

    Peter M. Neumann, The mathematical writings of Évariste Galois, European Mathematical Society, 2011.
    II.1 A theorem on continued fractions, pp. 35-47
    Proof of a theorem on periodic continued fractions
    "If one of the roots of an equation of arbitrary degree is an immediately periodic continued fraction,
    the equation will necessarily have another root that is likewise periodic."
    """
    r = (P + math.sqrt(D)) / Q
    r_conjugate = (P - math.sqrt(D)) / Q
    if D > 0 and r > 1 and -1 < r_conjugate < 0:
        return True
    return False


# Example usage:
if __name__ == "__main__":
    d = 19
    num_terms = 15
    precision = 30
    decimal.getcontext().prec = precision

    # Example
    periodic_cf = sqrt_periodic_cf(d)
    initial_part, periodic_part = periodic_cf.initial_part, periodic_cf.periodic_part
    print(f"Continued fraction for √{d}:")
    print(f"Initial part:\t{initial_part}")
    print(f"Periodic part:\t{periodic_part}")
    print()

    # Example
    print(f"Using {num_terms=} of the continued fraction for √{d}, the decimal expansion with {precision=} is:")
    cf_generator = periodic_cf_generator(sqrt_periodic_cf(d))
    cf = [a for a, _ in zip(cf_generator, range(num_terms))]
    convergent_list = cflib.cf_to_convergent_list(cf)
    fraction_tuple = convergent_list[-1]
    print(Decimal(fraction_tuple[0]) / Decimal(fraction_tuple[1]))
    print()
    
    # Example
    print(f"Pell equation x² - {d}y² = 1 has the following fundamental solution:")
    print("(x, y) = ", solve_pell_equation(d), sep="")
    print()

    # Example
    periodic_cf = sqrt_periodic_cf(d)
    a, b, c = periodic_cf_to_integral_quadratic_polynomial(periodic_cf)
    print(f"√{d} is a root of the following quadratic polynomial:\n{a}x² \t{b:+}x \t{c:+}")
    print()

    # Example
    initial_part, periodic_part = [], [1, 2]
    periodic_cf = PeriodicCF(initial_part, periodic_part)
    a, b, c = periodic_cf_to_integral_quadratic_polynomial(periodic_cf)
    D = b ** 2 - 4 * a * c
    P = -b
    Q = 2 * a
    print(f"The value of {periodic_cf} is the quadratic surd ({P}+√{D})/{Q} which has minimal polynomial:\n{a}x² \t{b:+}x \t{c:+}")
    print()

    # Example
    coefficients = (7, -8, -3)
    assert periodic_cf_to_integral_quadratic_polynomial(integral_quadratic_polynomial_to_periodic_cf(*coefficients)) == \
        coefficients

    # Example
    coefficients = (7, -8, -4)
    a, b, c = coefficients
    D = b ** 2 - 4 * a * c
    P = -b
    Q = 2 * a
    assert is_real_reduced_surd(P, Q, D) == True
    assert periodic_cf_to_integral_quadratic_polynomial(integral_quadratic_polynomial_to_periodic_cf(*coefficients)) == coefficients

    # Example
    d = 5
    print(f"Pell equation x² - {d}y² = 4 has the following fundamental solution:")
    x, y = solve_pell_equation(d)
    x *= 2
    y *= 2
    assert x ** 2 - d * y ** 2 == 4


