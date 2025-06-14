import math
import itertools as it
import functools as ft
import operator
from collections import abc
from fractions import Fraction
from numbers import Rational
import re
import typing

import cflib


def detM2(alpha: Rational, beta: Rational, gamma: Rational, delta: Rational) -> Rational:
    return alpha * delta - beta * gamma


class GL2Z(abc.Hashable):
    __slots__ = ("alpha", "beta", "gamma", "delta")

    def __init__(self: typing.Self, alpha: int, beta: int, gamma: int, delta: int) -> None:
        if not all(isinstance(x, int) for x in [alpha, beta, gamma, delta]):
            raise TypeError("alpha, beta, gamma, delta must all be integers.")
        det = detM2(alpha, beta, gamma, delta)
        if det not in [-1, 1]:
            raise ValueError(f"Determinant must -1 or 1: {det=}.")
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta

    def __repr__(self: typing.Self) -> str:
        return f"{type(self).__name__}({self.alpha}, {self.beta}, {self.gamma}, {self.delta})"
    
    def __eq__(self: typing.Self, other: typing.Self) -> bool:
        return self.alpha == other.alpha and self.beta == other.beta and self.gamma == other.gamma and self.delta == other.delta
    
    def __hash__(self: typing.Self) -> int:
        return hash((self.alpha, self.beta, self.gamma, self.delta))
    
    def __neg__(self: typing.Self) -> typing.Self:
        return type(self)(-self.alpha, -self.beta, -self.gamma, -self.delta)

    def __mul__(self: typing.Self, other: typing.Self) -> typing.Self:
        alpha = self.alpha * other.alpha + self.beta * other.gamma
        beta = self.alpha * other.beta + self.beta * other.delta
        gamma = self.gamma * other.alpha + self.delta * other.gamma
        delta = self.gamma * other.beta + self.delta * other.delta
        return type(self)(alpha, beta, gamma, delta)
    
    def inverse(self) -> typing.Self:
        det = self.det # either 1 or -1
        alpha = self.delta // det
        beta = -self.beta // det
        gamma = -self.gamma // det
        delta = self.alpha // det
        return type(self)(alpha, beta, gamma, delta)
    
    def __truediv__(self: typing.Self, other: typing.Self) -> typing.Self:
        return self * other.inverse()
    
    def __pow__(self: typing.Self, exponent: int) -> typing.Self:
        I = type(self)(1, 0, 0, 1)
        if exponent == 0:
            return I
        elif exponent > 0:
            return ft.reduce(operator.mul, (self for _ in range(exponent)), I)
        elif exponent < 0:
            self_inv = self.inverse()
            return ft.reduce(operator.mul, (self_inv for _ in range(abs(exponent))), I)
    
    def __iter__(self: typing.Self) -> abc.Iterator[int]:
        return iter((self.alpha, self.beta, self.gamma, self.delta))
    
    def __contains__(self: typing.Self, other: int) -> bool:
        return any(other == entry for entry in self)
    
    def __getitem__(self: typing.Self, index: int) -> int:
        return (*self, )[index]
    
    def __len__(self: typing.Self) -> int:
        return len(type(self).__slots__)
    
    def __abs__(self: typing.Self) -> int: # maximum absolute value of entries (ℓ∞ norm)
        return max(abs(entry) for entry in self)
    
    def transpose(self: typing.Self) -> typing.Self:
        return type(self)(self.alpha, self.gamma, self.beta, self.delta)
    
    @property
    def det(self: typing.Self) -> int:
        return detM2(*self)
    
    @property
    def trace(self: typing.Self) -> int:
        return self.alpha + self.delta
    
    @classmethod
    def N(cls, n: int) -> typing.Self:
        return cls(1, n, 0, 1) # T ** n


I = GL2Z(1, 0, 0, 1)
P = GL2Z(0, 1, 1, 0) # det == -1. P⁻¹ = P
R = GL2Z(0, 1, -1, 0) # S⁻¹
S = GL2Z(0, -1, 1, 0)
T = GL2Z(1, 1, 0, 1)
U = GL2Z(1, -1, 0, 1) # T⁻¹
V = GL2Z(0, -1, 1, 1) # S * T
W = GL2Z(1, 1, -1, 0) # V⁻¹


ALPHABET_DICT = {
    "I": I,
    "P": P,
    "R": R,
    "S": S,
    "T": T,
    "U": U,
    "V": V,
    "W": W
}

ALPHABET = ALPHABET_DICT.keys()

REVERSE_ALPHABET_DICT = {v: k for k, v in ALPHABET_DICT.items()}

MATRIX_ALPHABET = REVERSE_ALPHABET_DICT.keys()

"""
Wilhelm Magnus, Abraham Karrass, Donald Solitar,
Combinatorial group theory: Presentations of Groups in Terms of Generators and Relations,
second revised edition, Dover Publications, 1976.
Section 1.4, Problem 24, pp. 46-47.
Section 3.2, Theorem 3.2, p. 131.
Section 3.5, Corollary N4, p. 169.
"""

RELATIONS_GL2Z = {
    # length 1 word
    "I": "",
    # length 2 words, namely inverses
    "PP": "",
    "RS": "",
    "SR": "",
    "TU": "",
    "UT": "",
    "VW": "",
    "WV": "",
    # length 3 words
    ## non-conjugations
    "RRR": "S",
    "SSS": "R",
    ## conjugations
    "PRP": "S",
    "PSP": "R",
    "PTP": "TST",
    "PVP": "W",
    "PWP": "V",
    # higher length words
    "RRRR": "",
    "SSSS": "", # S⁴ = I
    "STSTST": "SS", # (ST)³ = S²
}

RELATIONS_GL2Z_REGEX = re.compile('|'.join(sorted(RELATIONS_GL2Z, key=len, reverse=True)))


def rewrite_word(word: str) -> str:
    """Apply one relation from left, if possible."""
    return RELATIONS_GL2Z_REGEX.sub(lambda m: RELATIONS_GL2Z[m.group()], word, count=1)


def reduce_word(word: str) -> str:
    """Repeatedly apply rewrite_once until fixed point reached."""
    while True:
        new_word = rewrite_word(word)
        if new_word == word:
            return word
        word = new_word


def minimum_matrix_list_from_matrix_alphabet_it(target_matrix: GL2Z, matrix_alphabet: list[GL2Z]) -> GL2Z:
    for word_length in it.count():
        matrix_word_iterator = it.product(matrix_alphabet, repeat=word_length)
        for matrix_word in matrix_word_iterator:
            if math.prod(matrix_word, start=I) == target_matrix:
                return list(matrix_word)


def minimum_word_from_alphabet_dp(target_matrix: GL2Z, alphabet: list[str]) -> str:
    product_dict = {(): I}
    for word_length in it.count():
        for letter_tuple in it.product(alphabet, repeat=word_length):
            word = "".join(letter_tuple)
            if word_length == 0:
                matrix_product = I
            else:
                initial_word, terminal_letter = word[:-1], word[-1]
                matrix_product = product_dict[initial_word] * ALPHABET_DICT[terminal_letter]
            product_dict[word] = matrix_product
            if matrix_product == target_matrix:
                return word
            

def minimum_word_from_alphabet_dp_max_len(target_matrix: GL2Z, alphabet: list[str], max_len: int) -> str:
    product_dict = {(): I}
    for word_length in range(max_len + 1):
        for letter_tuple in it.product(alphabet, repeat=word_length):
            word = "".join(letter_tuple)
            if word_length == 0:
                matrix_product = I
            else:
                initial_word, terminal_letter = word[:-1], word[-1]
                matrix_product = product_dict[initial_word] * ALPHABET_DICT[terminal_letter]
            product_dict[word] = matrix_product
            if matrix_product == target_matrix:
                return word
            

def word_to_matrix_list(word: str) -> list[GL2Z]:
    if not all(letter in ALPHABET for letter in word):
        raise ValueError(f"{word=} must belong to {ALPHABET=}.")
    matrix_list = []
    for letter in word:
        matrix_list.append(ALPHABET_DICT[letter])
    return matrix_list

def word_to_matrix(word: str) -> GL2Z:
    matrix_list = word_to_matrix_list(word)
    matrix_product = math.prod(matrix_list, start=I)
    return matrix_product

def matrix_list_to_word(matrix_list: list[GL2Z]) -> str:
    if not all(matrix in MATRIX_ALPHABET for matrix in matrix_list):
        raise ValueError(f"{matrix_list=} must belong to {MATRIX_ALPHABET=}.")
    word = ""
    for matrix in matrix_list:
        word += REVERSE_ALPHABET_DICT[matrix]
    return word


def max_word_len_linf(m: GL2Z) -> int:
    """
    Search radius for words in {P,S,T} equal to m. 
    Uses the ℓ∞ norm bound L = 3·⌈log_ϕ ‖core‖∞⌉ + 2, plus 1 if a single leading P is required.
    """

    EXTREME_AND_MEAN_RATIO = (1 + math.sqrt(5)) / 2

    prepend_P_flag = m.det == -1

    matrix_product = P * m if prepend_P_flag else m

    if matrix_product == I:
        return 1 if prepend_P_flag else 0

    matrix_product_norm = abs(matrix_product) # ℓ∞‐norm

    L = 3 * math.ceil(math.log(matrix_product_norm, EXTREME_AND_MEAN_RATIO)) + 2

    # For matrix_product_norm = 1 the formula gives L = 2, but S³ and R need 3.
    if matrix_product_norm == 1 and L < 3:
        L = 3

    return L + (1 if prepend_P_flag else 0)


def ball_SL2Z(radius: int) -> abc.Generator[GL2Z]:
    for alpha in range(-radius, radius + 1):
        for beta in range(-radius, radius + 1):
            if math.gcd(alpha, beta) == 1:
                eea = cflib.EEA(alpha, beta)
                gamma = -eea.bezout_y
                delta = eea.bezout_x
                matrix = GL2Z(alpha, beta, gamma, delta)
                yield matrix


def ball_GL2Z(radius: int) -> abc.Generator[GL2Z]:
    for alpha in range(-radius, radius + 1):
        for beta in range(-radius, radius + 1):
            if math.gcd(alpha, beta) == 1:
                eea = cflib.EEA(alpha, beta)
                gamma = -eea.bezout_y
                delta = eea.bezout_x
                matrix = GL2Z(alpha, beta, gamma, delta)
                yield matrix
                yield P * matrix


def upper_half_plane_action(matrix: GL2Z, tau: complex) -> complex:
    if matrix.det != 1:
        raise ValueError(f"{matrix=} must belong to SL₂(𝐙).")
    if tau.imag <= 0:
        raise ValueError(f"{tau=} must belong to upper half plane 𝓗.")
    return (matrix.alpha * tau + matrix.beta) / (matrix.gamma * tau + matrix.delta)


def transformation_to_fundamental_domain(tau: complex) -> tuple[list[GL2Z], GL2Z, list[int]]:
    """
    Henri Cohen, A Course in Computational Algebraic Number Theory.
    Algorithm 7.4.2, p. 395:
    "Given 𝜏 ∈ 𝓗, this algorithm outputs the unique 𝜏' equivalent to 𝜏 under the action of SL₂(𝐙) and which belongs to
    the standard fundamental domain 𝓕, as well as the matrix A ∈ SL₂(𝐙) such that 𝜏' = 𝜏."
    """
    A = I
    matrix_list = []
    exponent_list = []
    n = math.floor(tau.real + 1 / 2)
    tau = tau - n
    N = GL2Z(1, -n, 0, 1) # T ** (-n) == U ** n
    A = N * A
    matrix_list.append(N)
    exponent_list.append(n)
    m = abs(tau ** 2)
    while m < 1:
        tau = -tau.conjugate() / m
        A = S * A
        matrix_list.append(S)
        n = round(tau.real)
        tau = tau - n
        N = GL2Z(1, -n, 0, 1) # T ** (-n) == U ** n
        A = N * A
        matrix_list.append(N)
        exponent_list.append(n)
        m = abs(tau ** 2)
    return matrix_list, A, exponent_list


class GL2Q(abc.Hashable):
    def __init__(self, alpha: Rational, beta: Rational, gamma: Rational, delta: Rational) -> None:
        det = detM2(alpha, beta, gamma, delta)
        if det == 0:
            raise ValueError(f"Determinant must be nonzero: {det=}")
        self.alpha = Fraction(alpha)
        self.beta = Fraction(beta)
        self.gamma = Fraction(gamma)
        self.delta = Fraction(delta)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.alpha}, {self.beta}, {self.gamma}, {self.delta})"

    def __eq__(self: typing.Self, other: typing.Self) -> bool:
        return self.alpha == other.alpha and self.beta == other.beta and self.gamma == other.gamma and self.delta == other.delta
    
    def __hash__(self: typing.Self) -> int:
        return hash((self.alpha, self.beta, self.gamma, self.delta))
    
    def __neg__(self: typing.Self) -> typing.Self:
        return type(self)(-self.alpha, -self.beta, -self.gamma, -self.delta)

    def __mul__(self: typing.Self, other: typing.Self) -> typing.Self:
        alpha = self.alpha * other.alpha + self.beta * other.gamma
        beta = self.alpha * other.beta + self.beta * other.delta
        gamma = self.gamma * other.alpha + self.delta * other.gamma
        delta = self.gamma * other.beta + self.delta * other.delta
        return type(self)(alpha, beta, gamma, delta)
    
    def inverse(self) -> typing.Self:
        det = self.det
        alpha = Fraction(self.delta, det)
        beta = Fraction(-self.beta, det)
        gamma = Fraction(-self.gamma, det)
        delta = Fraction(self.alpha, det)
        return type(self)(alpha, beta, gamma, delta)
    
    def __truediv__(self: typing.Self, other: typing.Self) -> typing.Self:
        return self * other.inverse()
    
    def __pow__(self: typing.Self, exponent: int) -> typing.Self:
        I = type(self)(1, 0, 0, 1)
        if exponent == 0:
            return I
        elif exponent > 0:
            return ft.reduce(operator.mul, (self for _ in range(exponent)), I)
        elif exponent < 0:
            self_inv = self.inverse()
            return ft.reduce(operator.mul, (self_inv for _ in range(abs(exponent))), I)
    
    def __iter__(self: typing.Self) -> abc.Iterator[Fraction]:
        return iter((self.alpha, self.beta, self.gamma, self.delta))
    
    def __contains__(self: typing.Self, other: int) -> bool:
        return any(other == entry for entry in self)
    
    def __getitem__(self: typing.Self, index: int) -> Fraction:
        return (*self, )[index]
    
    def __len__(self: typing.Self) -> int:
        return len(type(self).__slots__)
    
    def transpose(self: typing.Self) -> typing.Self:
        return type(self)(self.alpha, self.gamma, self.beta, self.delta)
    
    @property
    def det(self: typing.Self) -> Fraction:
        return detM2(*self)
    
    @property
    def trace(self: typing.Self) -> Fraction:
        return self.alpha + self.delta
    
    def transpose_action_GL2Z(self, matrix: GL2Z) -> typing.Self:
        g = type(self)(*matrix)
        h = self
        return g.transpose() * h * g
    
    def action_GL2Z_on_real_quadratic_field_image(self, matrix: GL2Z) -> typing.Self:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        § 4.2.3 "The Matrix (or Regular) Representation of an Algebraic Number", p. 160.
        
        Given left regular representation L of 𝐐(√d) in GL₂(𝐐)
        x + y * √d
        alpha = x
        beta = x * d
        gamma = y
        delta = x

        This function implements action of GL₂(𝐙) on GL₂(𝐐) that is compatible
        with the action of GL₂(𝐙) by linear fractional transformations on 𝐐(√d).
        In other words, the diagram

        𝐐(√d)  -> 𝐐(√d)
         |         |
         v         v
        GL₂(𝐐) -> GL₂(𝐐)

        s in 𝐐(√d)
        g = (a, b, c, d) in GL₂(𝐙)
        g·s = (alpha s + beta)/(gamma s + delta) = s.lft_GL2Z(g)

        s  ———→  g·s
        │         │
       L│         │L
        ↓         ↓
        L(s) ———→ L(g·s)

        commutes.

        Given: self = L(s) for s = x + y√d.
        Return L(g·s) where g = (alpha, beta, gamma, delta) ∈ GL₂(𝐙).
        """
        # 1) recover x, y, and the field parameter d:
        x, dy, y, _ = self.alpha, self.beta, self.gamma, self.delta
        if y == 0:
            raise ValueError("y must be nonzero to uniquely determine d from L(s).")
        d = Fraction(dy, y) # d == int(d)

        alpha, beta, gamma, delta = matrix.alpha, matrix.beta, matrix.gamma, matrix.delta

        # L(a·s + b)
        # a·s + b = (a*x + b) + (a*y)·√d
        x_numerator = alpha * x + beta
        y_numerator = alpha * y
        result_numerator = type(self)(x_numerator, d * y_numerator, y_numerator, x_numerator)

        # L(c·s + d):
        # c·s + d = (c*x + d0) + (c*y)·√d
        x_denominator = gamma * x + delta
        y_denominator = gamma * y
        result_denominator = type(self)(x_denominator, d * y_denominator, y_denominator, x_denominator)

        return result_numerator * result_denominator.inverse()



if __name__ == "__main__":
    assert P ** (-1) == P

    assert S ** 2 == -I

    assert S * T == V

    assert V ** 3 == -I

    assert S * R == I and R * S == I
    assert S ** (-1) == R
    assert S ** 3 == R
    assert R ** 3 == S

    assert T * U == I and U * T == I
    assert T ** (-1) == U

    assert P * V * P == W

    m = U # any matrix in GL₂(𝐙)
    assert U / U == I

    m = U # any matrix in GL₂(𝐙)
    assert m ** (-1) == m.inverse()

    tau = complex(13.5, 0.3) # in 𝓗
    matrix_list, m, exponent_list = transformation_to_fundamental_domain(tau)
    tau_F = upper_half_plane_action(m, tau)
    assert -1 / 2 <= tau_F.real <= 1 / 2 and abs(tau_F) >= 1 # in 𝓕

    tau = complex(13.5, 0.3) # in 𝓗
    matrix_list, m, exponent_list = transformation_to_fundamental_domain(tau)
    assert ft.reduce(operator.mul, reversed(matrix_list), I) == m
    assert math.prod(reversed(matrix_list), start=I) == ft.reduce(operator.mul, reversed(matrix_list), I)

    tau = complex(13.5, 0.3) # in 𝓗
    matrix_list0, _, exponent_list = transformation_to_fundamental_domain(tau)
    matrix_list = []
    for exponent in exponent_list[0:1]:
        matrix_list.append(U ** exponent)
    for exponent in exponent_list[1:]:
        matrix_list.append(S)
        matrix_list.append(U ** exponent)
    assert matrix_list == matrix_list0

    m = U
    matrix_alphabet = [S, T]
    alphabet = ["S", "T"]
    matrix_list = minimum_matrix_list_from_matrix_alphabet_it(m, matrix_alphabet)
    word = minimum_word_from_alphabet_dp(m, alphabet)
    assert math.prod(matrix_list, start=I) == m
    assert word_to_matrix_list(word) == matrix_list
    assert matrix_list_to_word(matrix_list) == word
    assert word == "SSSTSTS"

    word = "SSTSSSTSSTSTSTSSTTSSSSSS" # any string from alphabet {I, P, S, R, T, U, V, W}
    reduced_word = reduce_word(word)
    assert word_to_matrix(reduced_word) == word_to_matrix(word)

    m = R
    max_word_len = max_word_len_linf(m)
    word  = minimum_word_from_alphabet_dp_max_len(target_matrix=m, alphabet=['P', 'S', 'T'], max_len=max_word_len)
    assert len(word) <= max_word_len

    radius = 5
    assert len(list(ball_GL2Z(radius))) == 2 * len(list(ball_SL2Z(radius)))
