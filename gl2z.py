import math
import itertools as it
import functools as ft
import operator
from collections import abc
from fractions import Fraction
from numbers import Rational
import re
import typing
import random

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
        if not isinstance(other, GL2Z):
            return NotImplemented
        return self.alpha == other.alpha and self.beta == other.beta and self.gamma == other.gamma and self.delta == other.delta
    
    def __hash__(self: typing.Self) -> int:
        return hash((self.alpha, self.beta, self.gamma, self.delta))
    
    def __neg__(self: typing.Self) -> typing.Self:
        return type(self)(-self.alpha, -self.beta, -self.gamma, -self.delta)

    def __mul__(self: typing.Self, other: typing.Self) -> typing.Self:
        if not isinstance(other, GL2Z):
            return NotImplemented
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
        if not isinstance(other, GL2Z):
            return NotImplemented
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
    
    @classmethod
    def elementary_matrix_interchange_rows(cls) -> typing.Self:
        """
        Charles C. Sims, Computation with finitely presented groups, Encyclopedia of Mathematics and Its Applications, volume 48,
        Cambridge University Press, 1994.
        Chapter 8: Abelian groups, pp. 319-382.
        
        Section 8.1, "Free abelian groups", p. 321:
        There are three types of integer row operations:
        (1) Interchange two rows.
        (2) Multiply a row by -1.
        (3) Add an integral multiple of one row to another row.

        Section 8.2, "Elementary matrices", p. 331:
        Proposition 2.2. Let B be an m-by-m integer matrix. The following are
        equivalent:
        (a) det£ = ±l.
        (b) B is in GL(m,Z).
        (c) B is row equivalent to the m-by-m identity matrix.
        (d) S(B) = Zm.
        (e) B is a product of elementary matrices.
        """
        return cls(0, 1, 1, 0)
        
    @classmethod
    def elementary_matrix_multiply_row(cls, i: int, sgn: int) -> typing.Self:
        if sgn not in [-1, 1]:
            raise ValueError(f"{sgn=} must be 1 or -1.")
        if i == 1:
            return cls(sgn, 0, 0, 1)
        elif i == 2:
            return cls(1, 0, 0, sgn)
        else:
            raise ValueError(f"Row number {i=} must be 1 or 2.")
        

    @classmethod
    def elementary_matrix_add_multiple_of_row(cls, target_row: int, multiplicand_row: int, c: int) -> typing.Self:
        if target_row == 1 and multiplicand_row == 2:
            return cls(1, c, 0, 1)
        elif target_row == 2 and multiplicand_row == 1:
            return cls(1, 0, c, 1)
        else:
            raise ValueError(f"Row numbers {target_row=} and {multiplicand_row=} must be distinct and 1 or 2.")
    
    @classmethod
    def elementary_matrix_interchange_columns(cls) -> typing.Self:
        return GL2Z(0, 1, 1, 0)

    @classmethod
    def elementary_matrix_multiply_column(cls, j: int, sgn: int) -> typing.Self:
        if sgn not in [-1, 1]:
            raise ValueError(f"{sgn=} must be ±1.")
        if j == 1:
            return GL2Z(sgn, 0, 0, 1)
        elif j == 2:
            return GL2Z(1, 0, 0, sgn)
        else:
            raise ValueError(f"Column number {j} must be 1 or 2.")

    @classmethod
    def elementary_matrix_add_multiple_of_column(cls, target_col: int, multiplicand_column: int, c: int) -> typing.Self:
        if target_col == 1 and multiplicand_column == 2:
            return GL2Z(1, 0, c, 1)
        elif target_col == 2 and multiplicand_column == 1:
            return GL2Z(1, c, 0, 1)
        else:
            raise ValueError("Column numbers must be distinct and either 1 or 2.")


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
    __slots__ = ("alpha", "beta", "gamma", "delta")

    def __init__(self, alpha: Rational, beta: Rational, gamma: Rational, delta: Rational) -> None:
        if not all(isinstance(x, Rational) for x in [alpha, beta, gamma, delta]):
            raise TypeError(f"{alpha=}, {beta=}, {gamma=}, {delta=} must all be rational.")
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
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, GL2Q):
            return NotImplemented
        return self.alpha == other.alpha and self.beta == other.beta and self.gamma == other.gamma and self.delta == other.delta
    
    def __hash__(self: typing.Self) -> int:
        return hash((self.alpha, self.beta, self.gamma, self.delta))
    
    def __neg__(self: typing.Self) -> typing.Self:
        return type(self)(-self.alpha, -self.beta, -self.gamma, -self.delta)

    def __mul__(self: typing.Self, other: typing.Self | GL2Z) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, GL2Q):
            return NotImplemented
        alpha = self.alpha * other.alpha + self.beta * other.gamma
        beta = self.alpha * other.beta + self.beta * other.delta
        gamma = self.gamma * other.alpha + self.delta * other.gamma
        delta = self.gamma * other.beta + self.delta * other.delta
        return type(self)(alpha, beta, gamma, delta)
    
    def __rmul__(self: typing.Self, other: typing.Self | GL2Z) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, GL2Q):
            return NotImplemented
        alpha = other.alpha * self.alpha + other.beta * self.gamma
        beta = other.alpha * self.beta + other.beta * self.delta
        gamma = other.gamma * self.alpha + other.delta * self.gamma
        delta = other.gamma * self.beta + other.delta * self.delta
        return type(self)(alpha, beta, gamma, delta)
    
    def inverse(self) -> typing.Self:
        det = self.det
        alpha = Fraction(self.delta, det)
        beta = Fraction(-self.beta, det)
        gamma = Fraction(-self.gamma, det)
        delta = Fraction(self.alpha, det)
        return type(self)(alpha, beta, gamma, delta)
    
    def __truediv__(self: typing.Self, other: typing.Self) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, GL2Q):
            return NotImplemented
        return self * other.inverse()
    
    def __rtruediv__(self: typing.Self, other: typing.Self) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, GL2Q):
            return NotImplemented
        return other * self.inverse()
    
    def __pow__(self: typing.Self, exponent: int) -> typing.Self:
        if exponent == 0:
            return type(self)(1, 0, 0, 1)
        elif exponent > 0:
            return ft.reduce(operator.mul, (self for _ in range(exponent)), type(self)(1, 0, 0, 1))
        elif exponent < 0:
            self_inv = self.inverse()
            return ft.reduce(operator.mul, (self_inv for _ in range(abs(exponent))), type(self)(1, 0, 0, 1))
    
    def __iter__(self: typing.Self) -> abc.Iterator[Fraction]:
        return iter((self.alpha, self.beta, self.gamma, self.delta))
    
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


class M2Z(abc.Hashable):
    """
    𝐙-algebra M_2(𝐙).
    """

    __slots__ = ("a11", "a12", "a21", "a22")

    def __init__(self, a11: int, a12: int, a21: int, a22: int) -> None:
        if not all(a == int(a) for a in [a11, a12, a21, a22]):
            raise TypeError(f"{a11=}, {a12=}, {a21=}, {a22=} must all be integers.")
        self.a11 = int(a11)
        self.a12 = int(a12)
        self.a21 = int(a21)
        self.a22 = int(a22)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.a11}, {self.a12}, {self.a21}, {self.a22})"
    
    def __eq__(self, other: typing.Self | GL2Z) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, M2Z):
            return NotImplemented
        return all(getattr(self, a) == getattr(other, a) for a in type(self).__slots__)
    
    def __hash__(self) -> int:
        return hash((self.a11, self.a12, self.a21, self.a22))

    def __neg__(self: typing.Self) -> typing.Self:
        return type(self)(-self.a11, -self.a12, -self.a21, -self.a22)
    
    def __add__(self: typing.Self, other: typing.Self | GL2Z) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, M2Z):
            return NotImplemented
        a11 = self.a11 + other.a11
        a12 = self.a12 + other.a12
        a21 = self.a21 + other.a21
        a22 = self.a22 + other.a22
        return type(self)(a11, a12, a21, a22)
    
    def __radd__(self: typing.Self, other: typing.Self | GL2Z) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, M2Z):
            return NotImplemented
        a11 = other.a11 + self.a11
        a12 = other.a12 + self.a12
        a21 = other.a21 + self.a21
        a22 = other.a22 + self.a22
        return type(self)(a11, a12, a21, a22)
    
    def __sub__(self: typing.Self, other: typing.Self | GL2Z) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, M2Z):
            return NotImplemented
        a11 = self.a11 - other.a11
        a12 = self.a12 - other.a12
        a21 = self.a21 - other.a21
        a22 = self.a22 - other.a22
        return type(self)(a11, a12, a21, a22)
    
    def __rsub__(self: typing.Self, other: typing.Self | GL2Z) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if not isinstance(other, M2Z):
            return NotImplemented
        a11 = other.a11 - self.a11
        a12 = other.a12 - self.a12
        a21 = other.a21 - self.a21
        a22 = other.a22 - self.a22
        return type(self)(a11, a12, a21, a22)

    def __mul__(self: typing.Self, other: typing.Self | GL2Z | int) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if isinstance(other, int):
            other = type(self)(other, 0, 0, other)
        if not isinstance(other, M2Z):
            return NotImplemented
        a11 = self.a11 * other.a11 + self.a12 * other.a21
        a12 = self.a11 * other.a12 + self.a12 * other.a22
        a21 = self.a21 * other.a11 + self.a22 * other.a21
        a22 = self.a21 * other.a12 + self.a22 * other.a22
        return type(self)(a11, a12, a21, a22)
    
    def __rmul__(self: typing.Self, other: typing.Self | GL2Z | int) -> typing.Self:
        if isinstance(other, GL2Z):
            other = type(self)(other.alpha, other.beta, other.gamma, other.delta)
        if isinstance(other, int):
            other = type(self)(other, 0, 0, other)
        if not isinstance(other, M2Z):
            return NotImplemented
        a11 = other.a11 * self.a11 + other.a12 * self.a21
        a12 = other.a11 * self.a12 + other.a12 * self.a22
        a21 = other.a21 * self.a11 + other.a22 * self.a21
        a22 = other.a21 * self.a12 + other.a22 * self.a22
        return type(self)(a11, a12, a21, a22)
    
    @property
    def det(self) -> int:
        return detM2(*self)
    
    def __iter__(self: typing.Self) -> abc.Iterator[int]:
        return iter((self.a11, self.a12, self.a21, self.a22))
    
    def __contains__(self: typing.Self, other: int) -> bool:
        return any(other == entry for entry in self)
    
    def __getitem__(self: typing.Self, index: int) -> int:
        return (*self, )[index]
    
    def __len__(self: typing.Self) -> int:
        return len(type(self).__slots__)
    
    def transpose(self: typing.Self) -> typing.Self:
        return type(self)(self.a11, self.a21, self.a12, self.a22)

    def __pow__(self: typing.Self, exponent: int) -> typing.Self:
        if exponent == 0 and self.det != 0:
            return type(self)(1, 0, 0, 1)
        elif exponent == 0 and self.det == 0:
            raise ValueError(f"Exponent must be positive for matrix with determinant 0: {self=}")
        elif exponent > 0:
            return ft.reduce(operator.mul, (self for _ in range(exponent)), type(self)(1, 0, 0, 1))
        elif exponent < 0:
            raise ValueError(f"Exponent must be nonnegative.")
    
    @property
    def trace(self: typing.Self) -> int:
        return self.a11 + self.a22
    
    def __abs__(self: typing.Self) -> int: # maximum absolute value of entries (ℓ∞ norm)
        return max(abs(entry) for entry in self)
    
    @property
    def gcd(self: typing.Self) -> int:
        return math.gcd(*self)
    
    def entry(self, row_index: int, column_index: int) -> int:
        if row_index == 1 and column_index == 1:
            return self.a11
        elif row_index == 1 and column_index == 2:
            return self.a12
        elif row_index == 2 and column_index == 1:
            return self.a21
        elif row_index == 2 and column_index == 2:
            return self.a22
        else:
            raise ValueError(f"Both {row_index=} and {column_index=} must each be 1 or 2.")
        
    def try_to_GL2Q(self) -> GL2Q:
        if self.det != 0:
            return GL2Q(self.a11, self.a12, self.a21, self.a22)
        else:
            return None
    
    @property
    def hadamard_det_bound(self) -> float:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Proposition 2.2.4 (Hadamard's inequality), p. 51.
        """
        return math.prod(math.sqrt(sum(self.entry(i, j) ** 2 for j in range(1, 3))) for i in range(1, 3))

    @classmethod
    def GL2Z_to_M2Z(cls, m: GL2Z) -> typing.Self:
        return cls(m.alpha, m.beta, m.gamma, m.delta)


def hnf_2x2(A: M2Z) -> tuple[GL2Z, M2Z]:
    """
    Charles C. Sims, Computation with finitely presented groups, Encyclopedia of Mathematics and Its Applications, volume 48,
    Cambridge University Press, 1994.
    Chapter 8: Abelian groups, pp. 319-382.
    p. 323, Procedure ROW_REDUCE for integer row reduction.

    Given 2x2 matrix A in M2Z, find U in GL2Z and H in M2Z such that H=UA
    where H is in row Hermite normal form.

    Returns U, H.
    """

    if isinstance(A, GL2Z):
        A = M2Z(A.alpha, A.beta, A.gamma, A.delta)
    if not isinstance(A, M2Z):
        raise TypeError(f"{A=} must be an instance of M2Z.")
    U = M2Z(1, 0, 0, 1) # I in M2Z
    H = A

    i, j = 1, 1
    while i <= 2 and j <= 2:

        if all(H.entry(k, j) == 0 for k in range(i, 3)):
            j += 1
            continue

        if i == 1:
            while (0 < abs(H.entry(1, j)) <= abs(H.entry(2, j))
                   or 0 < abs(H.entry(2, j)) <= abs(H.entry(1, j))):
                if abs(H.entry(1, j)) <= abs(H.entry(2, j)):
                    q = H.entry(2, j) // H.entry(1, j)
                    T = GL2Z.elementary_matrix_add_multiple_of_row(2, 1, -q)
                else:
                    q = H.entry(1, j) // H.entry(2, j)
                    T = GL2Z.elementary_matrix_add_multiple_of_row(1, 2, -q)
                U = T * U
                H = T * H

        if H.entry(i, j) == 0: # unique non-zero entry in the column is in the other row
            T = GL2Z.elementary_matrix_interchange_rows()
            U = T * U
            H = T * H

        # make the pivot positive
        if H.entry(i, j) < 0:
            T = GL2Z.elementary_matrix_multiply_row(i, -1)
            U = T * U
            H = T * H

        # clear entries above the pivot (only possible when i == 2)
        if i == 2 and H.entry(1, j) != 0:
            q = H.entry(1, j) // H.entry(i, j)
            T = GL2Z.elementary_matrix_add_multiple_of_row(1, 2, -q)
            U = T * U
            H = T * H

        i += 1
        j += 1

    return GL2Z(U.a11, U.a12, U.a21, U.a22), H


def snf_2x2(m: M2Z) -> tuple[GL2Z, M2Z, GL2Z]:
    """
    Charles C. Sims, Computation with finitely presented groups, Encyclopedia of Mathematics and Its Applications, volume 48,
    Cambridge University Press, 1994.
    Chapter 8: Abelian groups, pp. 319-382.
    
    Section 8.3, "Finitely generated abelian groups", p. 333:
    Two m-by-n integer matrices A and B are equivalent over Z if one can be
    obtained from the other by a sequence of integer row and column operations.
    This is the same as saying that there are matrices P and Q in GL(m, Z)
    and GL(n,Z), respectively, such that A = PBQ. Equivalence of matrices
    over Z is an equivalence relation.

    Just as the matrices in row Hermite normal form are representatives for
    the equivalence classes of integer matrices under integer row equivalence,
    there is an easily identified set of representatives under equivalence. An
    integer matrix A is in Smith normal form if for some r >= 0 the entries
    di = Aii, 1 <= i <= r, are positive, A has no other nonzero entries, and di
    divides di+1, 1 <= i < r. The matrix A of Example 3.1 is in Smith normal
    form.

    Example 3.1, p. 332:
    2 0 0 0 0
    0 4 0 0 0
    0 0 12 0 0

    Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
    Algorithm 2.4.14 (Smith Normal Form), p. 77.
    """

    # ------------------------------------------------------------------
    # 0Input sanity check
    # ------------------------------------------------------------------
    if not isinstance(m, M2Z):
        raise TypeError(f"{m=} must be an instance of M2Z.")

    # Working copies (mutable through rebinding)
    H: M2Z = m
    U: M2Z = M2Z(1, 0, 0, 1)   # left accumulator (rows)
    V: M2Z = M2Z(1, 0, 0, 1)   # right accumulator (columns)

    while True:
        # --------------------------------------------------------------
        # 1Ensure the pivot H[0,0] is non‑zero by swapping rows/columns
        # --------------------------------------------------------------
        if H.a11 == 0:
            if H.a12 != 0:  # swap the two columns
                T = GL2Z.elementary_matrix_interchange_columns()
                H = H * T
                V = V * T
            elif H.a21 != 0:  # swap the two rows
                T = GL2Z.elementary_matrix_interchange_rows()
                H = T * H
                U = T * U
            elif H.a22 != 0:  # bottom‑right the only non‑zero
                # swap rows then columns to bring it to (1,1)
                T = GL2Z.elementary_matrix_interchange_rows()
                H = T * H
                U = T * U
                T = GL2Z.elementary_matrix_interchange_columns()
                H = H * T
                V = V * T
            else:  # zero matrix
                break

        # --------------------------------------------------------------
        # 2Clear the sub‑diagonal entry H[1,0] via extended gcd (row ops)
        # --------------------------------------------------------------
        if H.a21 != 0:
            q = H.a11 // H.a21
            T = GL2Z.elementary_matrix_add_multiple_of_row(1, 2, -q)
            H = T * H
            U = T * U
            # swap rows so that |H.a11| ≤ |H.a21| for next pass
            T = GL2Z.elementary_matrix_interchange_rows()
            H = T * H
            U = T * U
            continue  # repeat

        # --------------------------------------------------------------
        # 3Clear the super‑diagonal entry H[0,1] (column ops)
        # --------------------------------------------------------------
        if H.a12 != 0:
            q = H.a11 // H.a12
            T = GL2Z.elementary_matrix_add_multiple_of_column(1, 2, -q)  # col1 ← col1 − q·col2
            H = H * T
            V = V * T
            # swap columns so that |H.a11| ≤ |H.a12|
            T = GL2Z.elementary_matrix_interchange_columns()
            H = H * T
            V = V * T
            continue  # repeat

        # --------------------------------------------------------------
        # 4We now have H = [[d, 0], [0, e]].  Adjust signs & divisibility.
        # --------------------------------------------------------------
        d, e = H.a11, H.a22

        # make d positive
        if d < 0:
            T = GL2Z.elementary_matrix_multiply_row(1, -1)
            H = T * H
            U = T * U
            d = -d

        # make e positive
        if e < 0:
            T = GL2Z.elementary_matrix_multiply_column(2, -1)
            H = H * T
            V = V * T
            e = -e

        # enforce d | e
        if d != 0 and e % d != 0:
            q = e // d
            T = GL2Z.elementary_matrix_add_multiple_of_column(2, 1, -q)
            H = H * T
            V = V * T
            # loop back – this re‑introduces H.a12 ≠ 0 or revises e
            continue

        # done – H is in Smith form
        break

    # ------------------------------------------------------------------
    # 5Return GL2Z versions of the accumulators
    # ------------------------------------------------------------------
    U_gl = GL2Z(U.a11, U.a12, U.a21, U.a22)
    V_gl = GL2Z(V.a11, V.a12, V.a21, V.a22)
    return U_gl, H, V_gl


            

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

    # ------------------------------------------------------------------
    # ❶  Edge case – a relation that kills the word completely
    #     PP  →  ε   (identity matrix)
    # ------------------------------------------------------------------
    word = "PP"
    assert reduce_word(word) == ""                  # rewrites to empty string
    assert word_to_matrix(word) == I               # PP = I in GL₂ℤ

    # ------------------------------------------------------------------
    # ❷  Classic conjugation identity – PRP  →  S
    # ------------------------------------------------------------------
    word = "PRP"
    assert reduce_word(word) == "S"                 # single rewrite
    assert word_to_matrix(word) == S                # same group element

    # ------------------------------------------------------------------
    # ❸  Random word: idempotence of reduce_word and matrix invariance
    # ------------------------------------------------------------------
    random.seed(503)                           # deterministic reproducibility
    alphabet = list(ALPHABET)

    rand_word1 = "".join(random.choices(alphabet, k=random.randint(5, 25)))
    reduced1   = reduce_word(rand_word1)

    assert reduce_word(reduced1) == reduced1        # fixed-point reached
    assert word_to_matrix(rand_word1) == word_to_matrix(reduced1)

    # ------------------------------------------------------------------
    # ❹  Second random word: reduction never lengthens the string
    # ------------------------------------------------------------------
    rand_word2 = "".join(random.choices(alphabet, k=random.randint(5, 25)))
    reduced2   = reduce_word(rand_word2)

    assert len(reduced2) <= len(rand_word2)         # relations are length-non-increasing
    assert word_to_matrix(rand_word2) == word_to_matrix(reduced2)


    m = R
    max_word_len = max_word_len_linf(m)
    word  = minimum_word_from_alphabet_dp_max_len(target_matrix=m, alphabet=['P', 'S', 'T'], max_len=max_word_len)
    assert len(word) <= max_word_len

    radius = 5
    assert len(list(ball_GL2Z(radius))) == 2 * len(list(ball_SL2Z(radius)))

    m = M2Z(1, 0, 12, -5)
    assert I * m == m and m * I == m

    max_size = 20
    a11, a12, a21, a22 = random.randint(1, max_size), random.randint(1, max_size), random.randint(1, max_size), random.randint(1, max_size)
    m_A = M2Z(a11, a22, a21, a22)
    assert abs(m_A.det) <= m_A.hadamard_det_bound

    max_size = 20
    a11, a12, a21, a22 = random.randint(1, max_size), random.randint(1, max_size), random.randint(1, max_size), random.randint(1, max_size)
    m_A = M2Z(a11, a22, a21, a22)
    m_U, m_H = hnf_2x2(m_A)
    assert m_H == m_U * m_A
    assert m_U.det in [-1, 1]
    assert m_H.det == m_A.det or m_H.det == -m_A.det


    
    # 1. Identity matrix – already in SNF
    m_U, m_D, m_V = snf_2x2(M2Z(1, 0, 0, 1))
    assert m_D == M2Z(1, 0, 0, 1) and m_D == m_U * M2Z(1,0,0,1) * m_V
    assert m_U.det in (-1, 1) and m_V.det in (-1, 1)

    # 2. Full-rank example (|det| = 2) – divisibility d₁ | d₂
    m = M2Z(3, 5, 7, 11)          # det = -2
    m_U, m_D, m_V = snf_2x2(m)
    assert m_D == M2Z(1, 0, 0, 2)   # expected SNF
    assert m_D == m_U * m * m_V

    # 3. Rank-one example – zero appears on the second diagonal
    m_A = M2Z(2, 4, 4, 8)           # det = 0, gcd = 2
    m_U, m_D, m_V = snf_2x2(m_A)
    assert m_D == M2Z(2, 0, 0, 0)   # SNF diag(2, 0)
    assert m_D == m_U * m_A * m_V



    # ---------------------------------------------------------------------------
    # 3.  Smith form of a unimodular 2×2 matrix is the identity.
    # ---------------------------------------------------------------------------
    g = P
    _, g_SNF_, _ = snf_2x2(M2Z.GL2Z_to_M2Z(g))
    assert g_SNF_ == M2Z(1, 0, 0, 1)             # SNF must be diag(1,1)

    # ---------------------------------------------------------------------------
    # 4.  First diagonal entry of the SNF is  gcd of *all* entries.
    # ---------------------------------------------------------------------------
    B = M2Z( 6, 10,
            9, 15)                             # gcd = 3
    U, D, V = snf_2x2(B)
    assert D.a11 == B.gcd      # d₁  =  gcd(entries)

