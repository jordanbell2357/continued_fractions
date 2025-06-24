import math
import functools as ft

import cflib


def modular_inverse(n: int, x: int) -> int:
    if math.gcd(n, x) > 1:
        return None
    else:
        eea = cflib.EEA(n, x)
        bezout_coefficient = eea.bezout_y
        return bezout_coefficient % n


def solve_linear_congruence(a: int, b: int, m: int) -> int:
    """
    Hardy and Wright, fifth edition.
    Theorem 57, Chapter V, p. 62.
    """
    return [x for x in range(m) if (a * x - b) % m == 0]


def crt_two_moduli(congruence_tuple_1: tuple[int, int], congruence_tuple_2: tuple[int, int]) -> int:
    # assumes moduli m1 and m2 are relatively prime, no assumptions on a1 and a2
    a1, n1 = congruence_tuple_1
    a2, n2 = congruence_tuple_2
    gcd = math.gcd(n1, n2)
    if gcd > 1:
        raise ValueError(f"{n1=} and {n2=} have {gcd=} > 0")
    # m1, m2 satisfy m1 * n1 + m2 * n2 == 1
    eea = cflib.EEA(n1, n2)
    m1, m2 = eea.bezout_x, eea.bezout_y
    a = a1 * m2 * n2 + a2 * m1 * n1
    lcm = n1 * n2
    return a % lcm, lcm


def crt(congruence_tuple_list: list[tuple[int, int]]) -> int:
    """
    Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
    Algorithm 1.3.12 (Inductive Chinese), p. 21.
    """
    return ft.reduce(crt_two_moduli, congruence_tuple_list, (1, 1))



if __name__ == "__main__":
    # Example: modular inverse
    n = 74
    x = 13
    assert math.gcd(n, x) == 1

    assert modular_inverse(n, x) * x % n == 1


    # Example: Chinese Remainder Theorem for two moduli
    congruence_tuple_1 = (4, 30)
    congruence_tuple_2 = (1, 1)
    assert math.gcd(congruence_tuple_1[1], congruence_tuple_2[1]) == 1

    assert crt_two_moduli(congruence_tuple_1, congruence_tuple_2) == congruence_tuple_1


    # Example: Chinese Remainder Theorem for two moduli
    congruence_tuple_1 = (2, 15)
    congruence_tuple_2 = (18, 23)
    assert math.gcd(congruence_tuple_1[1], congruence_tuple_2[1]) == 1

    a, n = crt_two_moduli(congruence_tuple_1, congruence_tuple_2)
    assert n == congruence_tuple_1[1] * congruence_tuple_2[1]
    assert a % congruence_tuple_1[1] == congruence_tuple_1[0]
    assert a % congruence_tuple_2[1] == congruence_tuple_2[0]


    # Example: Chinese remainder theorem for four moduli
    congruence_tuple_list = [(1, 2), (4, 9), (4, 11), (17, 35)]
    assert math.gcd(congruence_tuple_list[0][1], congruence_tuple_list[1][1], congruence_tuple_list[2][1], congruence_tuple_list[3][1]) == 1

    a, n = crt(congruence_tuple_list)
    assert n == math.prod(moduli for _, moduli in congruence_tuple_list)
    assert all(a % congruence_tuple[1] == congruence_tuple[0] for congruence_tuple in congruence_tuple_list)
