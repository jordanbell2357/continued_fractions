import math
import itertools as it
import functools as ft
import operator
from collections import abc
import re
import typing


def detM2(alpha: int, beta: int, gamma: int, delta: int) -> int:
    return alpha * delta - beta * gamma


class M2Z(abc.Hashable):
    __slots__ = ("alpha", "beta", "gamma", "delta")

    def __init__(self: typing.Self, alpha: int, beta: int, gamma: int, delta: int) -> None:
        if not all(isinstance(x, int) for x in [alpha, beta, gamma, delta]):
            raise TypeError("alpha, beta, gamma, delta must all be integers.")
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
    
    def __iter__(self: typing.Self) -> abc.Iterator:
        return iter((self.alpha, self.beta, self.gamma, self.delta))
    
    def __contains__(self: typing.Self, other: int) -> bool:
        return any(other == entry for entry in self)
    
    def __getitem__(self: typing.Self, index: int) -> int:
        return (*self, )[index]
    
    def __len__(self: typing.Self) -> int:
        return len(type(self).__slots__)
    
    # maximum absolute value of entries (ℓ∞ norm)
    def __abs__(self: typing.Self) -> int:
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
    def I(cls) -> typing.Self:
        return cls(1, 0, 0, 1)
    
    # det(P) == -1
    # P * P == I
    @classmethod
    def P(cls) -> typing.Self:
        return cls(0, 1, 1, 0)

    # R == S⁻¹
    @classmethod
    def R(cls) -> typing.Self:
        return cls(0, 1, -1, 0)

    @classmethod
    def S(cls) -> typing.Self:
        return cls(0, -1, 1, 0)
    
    @classmethod
    def T(cls) -> typing.Self:
        return cls(1, 1, 0, 1)
    
    # U == T⁻¹
    @classmethod
    def U(cls) -> typing.Self:
        return cls(1, -1, 0, 1)
    
    # V == S * T
    @classmethod
    def V(cls) -> typing.Self:
        return cls(0, -1, 1, 1)
    
    # W == V⁻¹
    @classmethod
    def W(cls) -> typing.Self:
        return cls(1, 1, -1, 0)
    
    def __str__(self: typing.Self) -> str:
        pass


class GL2Z(M2Z):
    """
    Wilhelm Magnus, Abraham Karrass, Donald Solitar,
    Combinatorial group theory: Presentations of Groups in Terms of Generators and Relations,
    second revised edition, Dover Publications, 1976.
    Section 1.4, Problem 24, pp. 46-47.
    Section 3.2, Theorem 3.2, p. 131.
    Section 3.5, Corollary N4, p. 169.
    """
    def __init__(self, alpha: int, beta: int, gamma: int, delta: int):
        det = detM2(alpha, beta, gamma, delta)
        if det not in [-1, 1]:
            raise ValueError(f"Determinant must -1 or 1: {det=}.")
        super().__init__(alpha, beta, gamma, delta)

    def inv(self) -> typing.Self:
        det = self.det # either 1 or -1
        alpha = self.delta // det
        beta = -self.beta // det
        gamma = -self.gamma // det
        delta = self.alpha // det
        return type(self)(alpha, beta, gamma, delta)
    
    def __truediv__(self: typing.Self, other: typing.Self) -> typing.Self:
        return self * other.inv()
    
    def __pow__(self: typing.Self, exponent: int) -> typing.Self:
        if exponent == 0:
            return type(self).I()
        elif exponent > 0:
            return ft.reduce(operator.mul, (self for _ in range(exponent)), type(self).I())
        elif exponent < 0:
            self_inv = self.inv()
            return ft.reduce(operator.mul, (self_inv for _ in range(abs(exponent))), type(self).I())


class SL2Z(GL2Z):
    def __init__(self: typing.Self, alpha: int, beta: int, gamma: int, delta: int) -> None:
        det = detM2(alpha, beta, gamma, delta)
        if det != 1:
            raise ValueError("Determinant must be 1.")
        super().__init__(alpha, beta, gamma, delta)


I = GL2Z.I()
P = GL2Z.P() # det == -1. P⁻¹ = P
R = SL2Z.R() # S⁻¹
S = SL2Z.S()
T = SL2Z.T()
U = SL2Z.U() # T⁻¹
V = SL2Z.V() # S * T
W = SL2Z.W() # V⁻¹


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


def upper_half_plane_action(matrix: SL2Z, tau: complex) -> complex:
    if tau.imag <= 0:
        raise ValueError(f"{tau=} must belong to upper half plane 𝓗.")
    return (matrix.alpha * tau + matrix.beta) / (matrix.gamma * tau + matrix.delta)


def transformation_to_fundamental_domain(tau: complex) -> tuple[list[SL2Z], SL2Z, list[int]]:
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
    N = SL2Z(1, -n, 0, 1) # T ** (-n) == U ** n
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
        N = SL2Z(1, -n, 0, 1) # T ** (-n) == U ** n
        A = N * A
        matrix_list.append(N)
        exponent_list.append(n)
        m = abs(tau ** 2)
    return matrix_list, A, exponent_list


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

def matrix_list_to_word(matrix_list: list[SL2Z]) -> str:
    if not all(matrix in MATRIX_ALPHABET for matrix in matrix_list):
        raise ValueError(f"{matrix_list=} must belong to {MATRIX_ALPHABET=}.")
    word = ""
    for matrix in matrix_list:
        word += REVERSE_ALPHABET_DICT[matrix]
    return word


PHI = (1 + math.sqrt(5)) / 2

def max_word_len_linf(m1: GL2Z, m2: GL2Z) -> int:
    """
    Guaranteed search radius for words in {P,S,T} solving  w·m1 = m2.
    Uses the ℓ∞ norm bound   L = 3·⌈log_ϕ ‖core‖∞⌉ + 2,  plus 1 if a
    single leading P is required.
    """
    # ── 1.  Do we have to prefix one P?  (determinant mismatch)
    prepend_P: bool = (m1.det != m2.det)

    # ── 2.  Comparison matrix that lies in SL₂(ℤ)
    core = (P * m2 if prepend_P else m2) * m1.inv()

    # If the core itself is I, no S,T letters are needed
    if core == I:
        return 1 if prepend_P else 0

    N = abs(core)                  # ℓ∞‐norm via __abs__

    # General logarithmic bound (Cohen-style Euclidean decomposition)
    L = 3 * math.ceil(math.log(N, PHI)) + 2

    # For N = 1 the formula gives L = 2, but S³ (and hence R) needs 3.
    if N == 1 and L < 3:
        L = 3

    return L + (1 if prepend_P else 0)



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

    m = U
    assert U / U == I

    m = U
    assert m ** (-1) == m.inv()

    tau = complex(13.5, 0.3) # in 𝓗
    matrix_list, m, exponent_list = transformation_to_fundamental_domain(tau)
    tau_F = upper_half_plane_action(m, tau)
    assert -1 / 2 <= tau_F.real <= 1 / 2 and abs(tau_F) >= 1 # in 𝓕

    tau = complex(13.5, 0.3) # in 𝓗
    matrix_list, m, exponent_list = transformation_to_fundamental_domain(tau)
    assert ft.reduce(operator.mul, reversed(matrix_list), SL2Z.I()) == m
    assert math.prod(reversed(matrix_list), start=SL2Z.I()) == ft.reduce(operator.mul, reversed(matrix_list), SL2Z.I())

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

    word = "SSTSSSTSSTSTSTSSTTSSSSSS"
    reduced_word = reduce_word(word)
    assert word_to_matrix(reduced_word) == word_to_matrix(word)
    

    m1 = I
    m2 = R
    max_word_len = max_word_len_linf(m1, m2)
    word  = minimum_word_from_alphabet_dp_max_len(
                target_matrix=m2 * m1.inv(),
                alphabet=['P', 'S', 'T'],
                max_len=max_word_len)
    print(max_word_len)
    print(word)

