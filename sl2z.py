import math
import functools as ft
import operator
import typing

def detM2(alpha: int, beta: int, gamma: int, delta: int) -> int:
    return alpha * delta - beta * gamma


class M2Z:
    def __init__(self, alpha: int, beta: int, gamma: int, delta: int):
        if not all(isinstance(x, int) for x in [alpha, beta, gamma, delta]):
            raise TypeError("alpha, beta, gamma, delta must all be integers")
        det = detM2(alpha, beta, gamma, delta)
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.det = det

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.alpha}, {self.beta}, {self.gamma}, {self.delta})"
    
    def __eq__(self, other: typing.Self) -> bool:
        return self.alpha == other.alpha and self.beta == other.beta and self.gamma == other.gamma and self.delta == other.delta
    
    def __neg__(self) -> typing.Self:
        return type(self)(-self.alpha, -self.beta, -self.gamma, -self.delta)

    def __int__(self) -> int:
        return self.det

    def __mul__(self, other: typing.Self) -> typing.Self:
        alpha = self.alpha * other.alpha + self.beta * other.gamma
        beta = self.alpha * other.beta + self.beta * other.delta
        gamma = self.gamma * other.alpha + self.delta * other.gamma
        delta = self.gamma * other.beta + self.delta * other.delta
        return type(self)(alpha, beta, gamma, delta)
    
    @classmethod
    def I(cls) -> typing.Self:
        return cls(1, 0, 0, 1)

    def __pow__(self, exponent: int) -> typing.Self:
        if exponent <= 0:
            raise ValueError(f"{exponent=} must be be positive integer.")
        matrix_product = type(self).I()
        for _ in range(exponent):
            matrix_product *= self
        return matrix_product
    

class GL2Z(M2Z):
    def __init__(self, alpha: int, beta: int, gamma: int, delta: int):
        det = detM2(alpha, beta, gamma, delta)
        if det not in [-1, 1]:
            raise ValueError(f"determinant must -1 or 1: {det=}")
        super().__init__(alpha, beta, gamma, delta)

    def inv(self) -> typing.Self:
        det = self.det # either 1 or -1
        alpha = self.delta // det
        beta = -self.beta // det
        gamma = -self.gamma // det
        delta = self.alpha // det
        return type(self)(alpha, beta, gamma, delta)


class SL2Z(GL2Z):
    def __init__(self, alpha, beta, gamma, delta):
        det = detM2(alpha, beta, gamma, delta)
        if det != 1:
            raise ValueError("determinant must be 1")
        super().__init__(alpha, beta, gamma, delta)

    def upper_half_plane_action(self, z: complex) -> complex:
        if z.imag <= 0:
            raise ValueError(f"{z=} must belong to upper half plane 𝓗")
        return (self.alpha * z + self.beta) / (self.gamma * z + self.delta)
    
    @classmethod
    def J(cls) -> typing.Self:
        return cls(0, -1, 1, 0)
    
    @classmethod
    def S(cls) -> typing.Self:
        return cls(0, -1, 1, 0)
    
    @classmethod
    def T(cls) -> typing.Self:
        return cls(1, 1, 0, 1)
    
    @classmethod
    def word_to_matrix(cls, word: str) -> typing.Self:
        I = cls.I()
        S = cls.S()
        T = cls.T()
        product_matrix = I
        for letter in word:
            if letter == "S":
                product_matrix *= S
            elif letter == "T":
                product_matrix *= T
        return product_matrix
    
    @classmethod
    def word_list_to_matrix(cls, word_list: list[str]) -> typing.Self:
        word = "".join(word_list)
        return cls.word_to_matrix(word)

    @classmethod
    def transformation_to_fundamental_domain(cls, tau: complex) -> tuple[list[typing.Self], typing.Self]:
        """
        Henri Cohen, A Course in Computational Algebraic Number Theory.
        Algorithm 7.4.2, p. 395:
        "Given 𝜏 ∈ 𝓗, this algorithm outputs the unique 𝜏' equivalent to 𝜏 under the action of SL₂(𝐙) and which belongs to
        the standard fundamental domain 𝓕, as well as the matrix A ∈ SL₂(𝐙) such that 𝜏' = 𝜏."
        """
        A = cls.I() # initialization
        matrix_list = []
        matrix_list.append(A)

        n = round(tau.real) # first pass
        tau = tau - n
        N = cls(1, -n, 0, 1)
        A = N * A
        matrix_list.append(N)
        m = abs(tau ** 2) # == tau * tau.conjugate()
        while m < 1: # iteration
            tau = -tau.conjugate() / m
            A = cls.J() * A
            matrix_list.append(cls.J())
            n = round(tau.real)
            tau = tau - n
            N = cls(1, -n, 0, 1)
            A = N * A
            matrix_list.append(N)
            m = abs(tau ** 2)
        return matrix_list, A


if __name__ == "__main__":
    assert SL2Z.S() * SL2Z.S() == -SL2Z.I()

    assert (SL2Z.S() * SL2Z.T()) ** 3 == -SL2Z.I()

    tau = complex(13.5, 0.3) # in 𝓗
    matrix_list, m = SL2Z.transformation_to_fundamental_domain(tau)

    assert ft.reduce(operator.mul, reversed(matrix_list), SL2Z.I()) == m
    assert math.prod(reversed(matrix_list), start=SL2Z.I()) == m

    tau_F = m.upper_half_plane_action(tau)
    assert abs(tau_F) >= 1 and -1 / 2 <= tau_F.real <= 1 / 2

    word = "STTSSTT"
    print(word)
    m = SL2Z.word_to_matrix("STTSSTT")
    print(m)

