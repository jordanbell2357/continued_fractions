import math
from fractions import Fraction
import decimal
import itertools as it

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
            non_periodic_part = expansion[:start_period]
            periodic_part = expansion[start_period:]
            return non_periodic_part, periodic_part
        else:
            triples_dict[(M, D, a)] = i
            expansion.append(a)


def pell_decimal_approximation(d, non_periodic_part, periodic_part, num_terms, num_decimals):
    decimal.getcontext().prec = num_decimals
    l = num_terms - len(non_periodic_part)
    partial_quotients = non_periodic_part + math.ceil(l / len(periodic_part)) * periodic_part
    convergents = cflib.cf_to_convergent_list(partial_quotients)
    convergent = Fraction(*convergents[num_terms])
    convergent_decimal = decimal.Decimal(convergent.numerator) / decimal.Decimal(convergent.denominator)
    sqrt_d_decimal = decimal.Decimal(d).sqrt()
    print(f"sqrt({d}) \t\t=", sqrt_d_decimal)
    print(f"p_{num_terms} / q_{num_terms} \t\t=", convergent_decimal)
    print(f"sqrt({d}) - p_{num_terms} / q_{num_terms} \t=", convergent_decimal - sqrt_d_decimal)


def pell_equation(d, non_periodic_part, periodic_part):
    r = len(periodic_part)
    if r % 2 == 0:
        partial_quotients = non_periodic_part + periodic_part
        fraction_list = cflib.cf_to_convergent_list(partial_quotients)
        x, y = fraction_list[r - 1]
    else:
        partial_quotients = non_periodic_part + 2 * periodic_part
        fraction_list = cflib.cf_to_convergent_list(partial_quotients)
        x, y = fraction_list[2 * r - 1]

    assert x**2 - d * y**2 == 1

    print(f"{d=}")
    print(f"Fundamental solution {x=}, {y=}")
    print(f"x^2 - {d} * y^2 = {x**2} - {d} * {y**2} =  {x**2 - d * y**2}")


def main(d: int, num_terms: int, num_decimals) -> None:
    non_periodic_part, periodic_part = sqrt_periodic_continued_fraction(d)

    print(f"Periodic continued fraction for sqrt({d})")
    print(f"non-periodic part \t= {non_periodic_part}")
    print(f"periodic part \t\t= {periodic_part}")

    print()

    print(f"Decimal approximation of sqrt({d}) using first {num_terms} partial quotients.")
    pell_decimal_approximation(d, non_periodic_part, periodic_part, num_terms, num_decimals)

    print()
    
    print(f"Pell equation x^2 - {d} * y^2 = 1")
    pell_equation(d, non_periodic_part, periodic_part)

    print()

    print(f"Partial quotients a_0,...,a_{num_terms-1} of sqrt({d})")
    cf = it.chain(non_periodic_part, it.cycle(periodic_part))
    print([a for a, _ in zip(cf, range(num_terms))])


# Example usage:
if __name__ == "__main__":
    d = 3
    num_terms = 30
    num_decimals = 30

    main(d, num_terms, num_decimals)