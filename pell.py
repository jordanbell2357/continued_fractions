import math
from fractions import Fraction
import decimal
from decimal import Decimal
import itertools as it
from collections import abc

import cflib


def sqrt_periodic_continued_fraction(d: int) -> tuple[list[int], list[int]]:
    a0 = math.isqrt(d)  # floor(sqrt(d))
    if a0 ** 2 == d:
        # d is a square, so sqrt(d) is an integer and the continued fraction has length 1
        return ([a0], [])
    # Initialize (M, D, a)
    M = 0
    D = 1
    a = a0
    expansion = [a0]
    # Keep track of triples (M, D, a) to detect repetition
    triples_dict = {}
    # Map: (M, D, a) -> index in the expansion
    i = 0
    while True:
        M = D * a - M
        D = (d - M**2) // D
        a = (a0 + M) // D
        i += 1
        # If we see a triple (M, D, a) twice, we've found the start of the period
        if (M, D, a) in triples_dict:
            # The terms from index triples_dict[(M,D,a)] to index i-1 are the repeating cycle
            start_period = triples_dict[(M, D, a)]
            initial_part = expansion[:start_period]
            periodic_part = expansion[start_period:]
            return initial_part, periodic_part
        else:
            triples_dict[(M, D, a)] = i
            expansion.append(a)


def sqrt_periodic_continued_fraction_generator(d: int) -> abc.Generator[int]:
    initial_part, periodic_part = sqrt_periodic_continued_fraction(d)
    cf = it.chain(initial_part, it.cycle(periodic_part))
    yield from cf


def solve_pell_equation(d: int) -> tuple[int, int]:
    initial_part, periodic_part = sqrt_periodic_continued_fraction(d)
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


def periodic_cf_to_quadratic_coefficients(initial_part: list[int],
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
    return a, b, c


# Example usage:
if __name__ == "__main__":
    d = 5
    num_terms = 30
    precision = 30

    decimal.getcontext().prec = precision

    initial_part, periodic_part = sqrt_periodic_continued_fraction(d)
    print(f"Continued fraction for sqrt({d})")
    print(f"Initial part:\t{initial_part}")
    print(f"Periodic part:\t{periodic_part}")

    print()
    
    print(f"Pell equation x² - {d}y² = 1 has following fundamental solution:")
    print("(x, y)=", solve_pell_equation(d), sep="")

    print()

    print(f"Using {num_terms=} of continued fraction for sqrt({d}), the decimal expansion with {precision=} is:")
    cf_generator = sqrt_periodic_continued_fraction_generator(d)
    cf = [a for a, k in zip(cf_generator, range(num_terms))]
    convergent_list = cflib.cf_to_convergent_list(cf)
    fraction_tuple = convergent_list[-1]
    print(Decimal(fraction_tuple[0]) / Decimal(fraction_tuple[1]))

    print()

    initial_part, periodic_part = sqrt_periodic_continued_fraction(d)
    a, b, c = periodic_cf_to_quadratic_coefficients(initial_part, periodic_part)
    print(f"""sqrt({d}) is a root of the following quadratic polynomial:\n{a}x² \t{b:+}x \t{c:+}""")

    print()

    initial_part, periodic_part = [], [1]
    a, b, c = periodic_cf_to_quadratic_coefficients(initial_part, periodic_part)
    print(f"The value of the continued fraction is a root of the following quadratic polynomial:\n{a}x² \t{b:+}x \t{c:+}")

    print()
