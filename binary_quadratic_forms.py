import math
from fractions import Fraction
import itertools as it
import typing

import prime_numbers
import pell


class RealQuadraticNumber:
    def __init__(self, d: int, x: Fraction, y: Fraction) -> None:
        if not isinstance(d, int):
            raise TypeError("d must be integer")
        if not d > 0:
            raise ValueError("d must be > 0")
        prime_factor_counter_d = prime_numbers.make_prime_factor_counter(d)
        squarefull_part_d = 1
        squarefree_part_d = 1
        for p, k in prime_factor_counter_d.items():
            squarefull_part_d *= p ** (k // 2 * 2)
            squarefree_part_d *= p ** (k % 2)
        self.d = squarefree_part_d
        self.x = x
        self.y = y * math.isqrt(squarefull_part_d)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.d}, {self.x}, {self.y})"
    
    def __eq__(self, other: typing.Self) -> bool:
        return self.d == other.d and self.x == other.x and self.y == other.y
    
    def __neg__(self) -> typing.Self:
        return type(self)(self.d, -self.x, -self.y)
    
    def __add__(self, other: typing.Self) -> typing.Self:
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d, self.x + other.x, self.y + other.y)
    
    def __sub__(self, other: typing.Self) -> typing.Self:
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d, self.x - other.x, self.y - other.y)
    
    def __mul__(self, other: typing.Self) -> typing.Self:
        if isinstance(other, int):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d,
                          self.x * other.x + self.y * other.y * self.d,
                          self.x * other.y + self.y * other.x)
    
    def __rmul__(self, other: typing.Self) -> typing.Self:
        if isinstance(other, int):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d,
                          self.x * other.x + self.y * other.y * self.d,
                          self.x * other.y + self.y * other.x)
    
    def __truediv__(self, other: typing.Self) -> typing.Self:
        dividend = self
        divisor = other
        if isinstance(divisor, int):
            divisor = type(dividend)(dividend.d, divisor, 0)
        if dividend.d != divisor.d:
            raise ValueError("dividend and divisor must be members of same field")
        return type(dividend)(dividend.d,
                          dividend.x * divisor.x + dividend.y * divisor.y * dividend.d,
                          dividend.x * divisor.y + dividend.y * divisor.x)
    
    def __rtruediv__(self, other: typing.Self) -> typing.Self:
        dividend = other
        divisor = self
        if isinstance(dividend, int):
            dividend = type(divisor)(divisor.d, dividend, 0)
        if dividend.d != divisor.d:
            raise ValueError("dividend and divisor must be members of same field")
        return type(dividend)(dividend.d,
                          dividend.x * divisor.x + dividend.y * divisor.y * dividend.d,
                          dividend.x * divisor.y + dividend.y * divisor.x)
    
    def __pow__(self, exponent: int) -> typing.Self:
        if not isinstance(exponent, int):
            raise TypeError(f"{exponent=} must be integer")
        if exponent <= 0 and self.norm == 0:
            raise ValueError(f"{exponent=} must be positive for non-invertible field elements")
        if exponent == 0:
            return type(self)(self.d, 1, 0)
        product_number = exponent // abs(exponent)
        for k in range(abs(exponent)):
            product_number = product_number * self
        return product_number

    def __float__(self) -> float:
        return self.x + self.y * math.sqrt(self.d)
    
    def __abs__(self) -> float:
        return abs(float(self))
    
    def __str__(self) -> str:
        return f"{self.x}\t{self.y:+} * √{self.d}"

    def conjugate(self) -> typing.Self:
        return type(self)(self.d, self.x, -self.y)

    @property
    def norm(self) -> Fraction:
        return self.x ** 2 - self.d * self.y ** 2
    
    @property
    def trace(self) -> Fraction:
        return 2 * self.x
    
    @property
    def is_integral(self) -> bool:
        return self.norm == int(self.norm) and self.trace == int(self.trace)


class RealQuadraticField:
    def __init__(self, d: int) -> None:
        if not isinstance(d, int):
            return TypeError("d must be integer")
        if not d > 0:
            return ValueError("d must be > 0")
        self.d = d

    def integral_basis(self) -> tuple[RealQuadraticNumber, RealQuadraticNumber]:
        if self.d % 4 == 1:
            return RealQuadraticNumber(self.d, Fraction(1, 1), Fraction(0, 1)), RealQuadraticNumber(self.d, Fraction(1, 2), Fraction(1, 2))
        elif self.d % 4 in [2, 3]:
            return RealQuadraticNumber(self.d, Fraction(1, 1), Fraction(0, 1)), RealQuadraticNumber(self.d, 1, 1)
        
    @property
    def discriminant(self) -> int:
        b1, b2 = self.integral_basis()
        a = b1
        b = b2
        c = b1.conjugate()
        d = b2.conjugate()
        determinant = a * d - b * c
        return int(float(determinant * determinant))
    
    def fundamental_unit(self) -> RealQuadraticNumber:
        x, y = pell.solve_pell_equation(self.d)
        return RealQuadraticNumber(d, x, y)
    
    def regulator(self) -> float:
        return math.log(float(self.fundamental_unit()))


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

    def __int__(self) -> int:
        return self.det

    def __matmul__(self, other: typing.Self) -> typing.Self:
        alpha = self.alpha * other.alpha + self.beta * other.gamma
        beta = self.alpha * other.beta + self.beta * other.delta
        gamma = self.gamma * other.alpha + self.delta * other.gamma
        delta = self.gamma * other.beta + self.delta * other.delta
        return type(self)(alpha, beta, gamma, delta)
    

class GL2Z(M2Z):
    def __init__(self, alpha: int, beta: int, gamma: int, delta: int):
        det = detM2(alpha, beta, gamma, delta)
        if det == 0:
            raise ValueError("determinant must not be 0")
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

    @classmethod
    def I(cls) -> typing.Self:
        return cls(1, 0, 0, 1)
    
    @classmethod
    def S(cls) -> typing.Self:
        return cls(0, 1, -1, 0)
    
    @classmethod
    def T(cls) -> typing.Self:
        return cls(1, 1, 0, 1)
    
    @classmethod
    def word_list_to_matrix(cls, word_list: list[str]) -> typing.Self:
        I = cls.I()
        S = cls.S()
        T = cls.T()
        product_matrix = I
        for word in word_list:
            for letter in word:
                if letter == "S":
                    product_matrix @= S
                elif letter == "T":
                    product_matrix @= T
        return product_matrix


class IndefiniteBQF:
    """
    Henri Cohen, A Course in Computation Algebraic Number Theory, Chapter 5, Algorithms for Quadratic Fields:
    Definition 5.2.3, p. 225, for binary quadratic forms.
    Definition 5.6.2, p. 262, for reduced indefinite binary quadratic forms.
    p. 263, for (a, b, c) being reduced if and only if 0 < (-b+√D)/(2|a|) < 1 and (b+√D)/(2|a|) > 1.
    Definition 5.6.4, p. 263 for reduction operator.
    Algorithm 5.6.5, p. 263 for reduction algorithm for indefinite quadratic forms.
    """

    def __init__(self, a: int, b: int, c: int) -> None:
        discriminant = b ** 2 - 4 * a * c
        if not all(isinstance(x, int) for x in [a, b, c]):
            raise TypeError("a, b, c must all be integers")
        if not discriminant > 0:
            raise ValueError("discriminant must be positive")
        self.a = a
        self.b = b
        self.c = c
        self.D = discriminant

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.a}, {self.b}, {self.c})"

    def __eq__(self, other: typing.Self) -> bool:
        return self.a == other.a and self.b == other.b and self.c == other.c
    
    def __str__(self) -> str:
        return f"{self.a}x²\t{self.b:+}xy\t{self.c:+}y².\tD={self.D}"
    
    def transform(self, matrix: SL2Z) -> typing.Self:
        a = self.a * matrix.alpha ** 2 + self.b * matrix.alpha * matrix.gamma + self.c * matrix.gamma ** 2
        b = 2 * self.a * matrix.alpha * matrix.beta + \
            self.b * matrix.alpha * matrix.delta + \
            self.b * matrix.beta * matrix.gamma + \
            2 * self.c * matrix.gamma * matrix.delta
        c = self.a * matrix.beta ** 2 + self.b * matrix.beta * matrix.delta + self.c * matrix.delta ** 2
        return type(self)(a, b, c)
    
    @property
    def gcd(self) -> int:
        return math.gcd(self.a, self.b, self.c)
    
    @property
    def is_primitive(self) -> bool:
        return self.gcd == 1

    @property
    def is_reduced(self) -> bool:
        """
        Cohen 5.6.2: (a,b,c) is reduced iff
           0 < (-b + √D)/(2|a|) < 1    and    (b + √D)/(2|a|) > 1.
        Equivalently, if alpha = (-b + √D)/(2a), beta = (-b - √D)/(2a), then
           0 < alpha < 1  and  beta < -1.
        """
        sqrtD = math.sqrt(self.D)
        alpha = (-self.b + sqrtD) / (2 * self.a)
        beta = (-self.b - sqrtD) / (2 * self.a)
        return (0 < alpha < 1) and (beta < -1)
    

    def reduce_with_m(self) -> tuple[typing.Self, int]:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Chapter 5, Algorithms for Quadratic Fields:
        Cohen Algorithm 5.6.5.
        1) Compute t0 = (-b mod 2|c| in [0, 2|c|-1].
        2) If |c| > sqrt(D), choose r so that -|c| < r ≤ |c|.
           If |c| < sqrt(D), choose r so that sqrt(D) - 2|c| < r < sqrt(D).
           In either case, r ≡ -b (mod 2c).
        3) Record m.
        4) Return BQF (c, r, (r^2 - D)/(4c)) and m.
        """

        _, b, c, D = self.a, self.b, self.c, self.D # self.a not used
        sqrtD = math.sqrt(D)
        abs_c = abs(c)

        # Step 1
        t0 = ((-b) % (2 * abs_c))

        if abs_c > sqrtD:
            # We need -|c| < r ≤ |c|  AND  r ≡ t0 (mod 2|c|)
            # The interval I = (-abs_c, abs_c] has length 2|c|.
            # Exactly one integer in I is ≡ t0 (mod 2|c|).
            # r = t0 + 2|c|k.
            # We want the unique k such that -abs_c < t0 + 2|c|k ≤ abs_c.
            # Equivalent to k = floor((abs_c - t0) / (2|c|)).
            k = (abs_c - t0) // (2 * abs_c)
            r = t0 + 2 * abs_c * k
        else:
            # We need  (sqrtD - 2|c|) < r < sqrtD,  AND  r ≡ t0 (mod 2|c|)
            # Let lower = sqrtD - 2|c|,  upper = sqrtD.
            lower = sqrtD - 2 * abs_c
            upper = sqrtD
            # We look for integer r = t0 + 2|c|·k lying in (lower, upper).
            # Since upper - lower = 2|c|, there is exactly one such integer.
            # Solve k ≥ (lower - t0)/(2|c|).  Let k0 = ceil((lower - t0)/(2|c|)).
            k0 = math.ceil((lower - t0) / (2 * abs_c))
            r = t0 + 2 * abs_c * k0

        # Step 3
        m = (b + r) // (2 * c)   # must be an integer exactly

        # Step 4
        a1 = c
        b1 = r
        c1 = (r * r - D) // (4 * c)
        return type(self)(a1, b1, c1), m

    
    def reduced_with_m_list(self) -> tuple[typing.Self, list[int]]:
        bqf = self
        m_list = []
        while not bqf.is_reduced:
            bqf, mi = bqf.reduce_with_m()
            m_list.append(mi)
        return (bqf, m_list)

    def real_quadratic_number(self) -> RealQuadraticNumber:
        return RealQuadraticNumber(self.D, Fraction(-self.b, 2 * abs(self.a)), Fraction(1, 2 * abs(self.a)))

    def evaluate(self, x: int, y: int) -> int:
        return self.a * x ** 2 + self.b * x * y + self.c * y **2
    
    @staticmethod
    def m_list_to_word_list(m_list: list[int]) -> list[str]:
        word_list = []
        for m in m_list:
            word = "SSS" + "T" * m
            word_list.append(word)
        return word_list



if __name__ == "__main__":
    d = 17

    assert RealQuadraticField(d).fundamental_unit().norm == 1

    assert float(RealQuadraticField(d).fundamental_unit()) > 1

    m = 13

    d = 4 * m + 1
    assert RealQuadraticField(d).discriminant == d

    d = 4 * m + 2
    assert RealQuadraticField(d).discriminant == 4 * d

    d = 4 * m + 3
    assert RealQuadraticField(d).discriminant == 4 * d

    d = 4 * m + 1
    assert RealQuadraticNumber(d, Fraction(1, 2), Fraction(1, 2)).is_integral == True

    d = 4 * m + 2
    assert RealQuadraticNumber(d, Fraction(1, 2), Fraction(1, 2)).is_integral == False

    d = 4 * m + 3
    assert RealQuadraticNumber(d, Fraction(1, 2), Fraction(1, 2)).is_integral == False

    bqf = IndefiniteBQF(2, 0, -6)

    assert bqf.is_reduced and (0 < float(bqf.real_quadratic_number()) < 1 and -float(bqf.real_quadratic_number().conjugate()) > 1) or \
        not bqf.is_reduced and not (0 < float(bqf.real_quadratic_number()) < 1 and -float(bqf.real_quadratic_number().conjugate()) > 1)


    bqf = IndefiniteBQF(3, 11, 2)

    reduced_bqf, m_list = bqf.reduced_with_m_list()
    word_list = IndefiniteBQF.m_list_to_word_list(m_list)
    product_matrix = SL2Z.word_list_to_matrix(word_list)
    assert bqf.transform(product_matrix) == reduced_bqf

    d = 19
    fundamental_unit_d = RealQuadraticField(d).fundamental_unit()
    for k in range(-5, 5):
        print(f"ε^{k}:", fundamental_unit_d ** k, sep="\t")