import math
from fractions import Fraction
import itertools as it
from collections import abc
import typing

import sl2z
import prime_numbers
import pell


class RealQuadraticNumber:
    def __init__(self: typing.Self, d: int, x: Fraction, y: Fraction) -> None:
        if not isinstance(d, int):
            raise TypeError("d must be integer")
        if not d > 0:
            raise ValueError("d must be > 0")
        squarefull_part_d, squarefree_part_d = prime_numbers.squarefull_and_squarefree_parts(d)
        self.d = squarefree_part_d
        self.x = x
        self.y = y * math.isqrt(squarefull_part_d)

    def __repr__(self: typing.Self) -> str:
        return f"{type(self).__name__}({self.d}, {self.x}, {self.y})"
    
    def __eq__(self: typing.Self, other: typing.Self) -> bool:
        return self.d == other.d and self.x == other.x and self.y == other.y
    
    def __neg__(self: typing.Self) -> typing.Self:
        return type(self)(self.d, -self.x, -self.y)
    
    def __add__(self: typing.Self, other: typing.Self) -> typing.Self:
        if isinstance(other, int) or isinstance(other, Fraction):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d, self.x + other.x, self.y + other.y)
    
    def __radd__(self: typing.Self, other: typing.Self) -> typing.Self:
        if isinstance(other, int) or isinstance(other, Fraction):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d, other.x + self.x, other.y + self.y)
    
    def __sub__(self: typing.Self, other: typing.Self) -> typing.Self:
        if isinstance(other, int) or isinstance(other, Fraction):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d, self.x - other.x, self.y - other.y)
    
    def __rsub__(self: typing.Self, other: typing.Self) -> typing.Self:
        if isinstance(other, int) or isinstance(other, Fraction):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d, other.x - self.x, other.y - self.y)
    
    def __mul__(self: typing.Self, other: typing.Self) -> typing.Self:
        if isinstance(other, int) or isinstance(other, Fraction):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d,
                          self.x * other.x + self.y * other.y * self.d,
                          self.x * other.y + self.y * other.x)
    
    def __rmul__(self: typing.Self, other: typing.Self) -> typing.Self:
        if isinstance(other, int) or isinstance(other, Fraction):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("algebraic numbers must be members of same field")
        return type(self)(self.d,
                          self.x * other.x + self.y * other.y * self.d,
                          self.x * other.y + self.y * other.x)
    
    def __truediv__(self: typing.Self, other: typing.Self) -> typing.Self:
        dividend = self
        divisor = other
        if isinstance(divisor, int) or isinstance(divisor, Fraction):
            divisor = type(dividend)(dividend.d, divisor, 0)
        if dividend.d != divisor.d:
            raise ValueError("dividend and divisor must be members of same field")
        divisor_norm = divisor.norm
        if divisor_norm == 0:
            raise ZeroDivisionError(f"divisor must have nonzero norm: {divisor_norm=}")
        return type(dividend)(dividend.d,
                          Fraction(dividend.x * divisor.x - dividend.y * divisor.y * dividend.d, divisor_norm),
                          Fraction(-dividend.x * divisor.y + dividend.y * divisor.x, divisor_norm))
    
    def __rtruediv__(self: typing.Self, other: typing.Self) -> typing.Self:
        dividend = other
        divisor = self
        if isinstance(dividend, int) or isinstance(divisor, Fraction):
            dividend = type(divisor)(divisor.d, dividend, 0)
        if dividend.d != divisor.d:
            raise ValueError("dividend and divisor must be members of same field")
        divisor_norm = divisor.norm
        if divisor_norm == 0:
            raise ZeroDivisionError(f"divisor must have nonzero norm: {divisor_norm=}")
        return type(dividend)(dividend.d,
                          Fraction(dividend.x * divisor.x - dividend.y * divisor.y * dividend.d, divisor_norm),
                          Fraction(-dividend.x * divisor.y + dividend.y * divisor.x, divisor_norm))
    
    def __pow__(self: typing.Self, exponent: int) -> typing.Self:
        if not isinstance(exponent, int):
            raise TypeError(f"{exponent=} must be integer")
        if exponent <= 0 and self.norm == 0:
            raise ZeroDivisionError(f"{exponent=} must be positive for non-invertible field elements")
        if exponent == 0:
            return type(self)(self.d, 1, 0)
        u = exponent // abs(exponent)
        multiplicand = self if u == 1 else 1 / self
        product_number = type(self)(self.d, 1, 0)
        for k in range(abs(exponent)):
            product_number = product_number * multiplicand
        return product_number

    def __float__(self: typing.Self) -> float:
        return self.x + self.y * math.sqrt(self.d)
    
    def __abs__(self: typing.Self) -> float:
        return abs(float(self))
    
    def __str__(self: typing.Self) -> str:
        return f"{self.x}\t{self.y:+} * √{self.d}"

    def conjugate(self: typing.Self) -> typing.Self:
        return type(self)(self.d, self.x, -self.y)

    @property
    def norm(self: typing.Self) -> Fraction:
        return self.x ** 2 - self.d * self.y ** 2
    
    @property
    def trace(self: typing.Self) -> Fraction:
        return 2 * self.x
    
    @property
    def is_integral(self: typing.Self) -> bool:
        return self.norm == int(self.norm) and self.trace == int(self.trace)


class RealQuadraticField:
    def __init__(self: typing.Self, d: int) -> None:
        if not isinstance(d, int):
            raise TypeError("d must be integer")
        if not d > 0:
            raise ValueError("d must be > 0")
        _, squarefree_part_d = prime_numbers.squarefull_and_squarefree_parts(d)
        self.d = squarefree_part_d
        # D = discriminant
        self.D = self.d if self.d % 4 == 1 else 4 * self.d
    
    def fundamental_unit(self: typing.Self) -> RealQuadraticNumber:
        x, y = pell.solve_pell_equation(self.d)
        return RealQuadraticNumber(self.d, x, y)

    @property
    def regulator(self: typing.Self) -> float:
        return math.log(float(self.fundamental_unit()))
    
    @property
    def primitive_real_dirichlet_character(self: typing.Self) -> abc.Callable:
        return lambda n: prime_numbers.kronecker_symbol(self.D, n)


class IndefiniteBQF:
    """
    Henri Cohen, A Course in Computation Algebraic Number Theory, Chapter 5, Algorithms for Quadratic Fields:
    Proposition 5.1.1, p. 223, for integral basis and discriminant of real quadratic fields.
    Definition 5.2.3, p. 225, for binary quadratic forms.
    Definition 5.6.2, p. 262, for reduced indefinite binary quadratic forms.
    p. 263, for (a, b, c) being reduced if and only if 0 < (-b + √D) / (2|a|) < 1 and (b + √D) / (2|a|) > 1.
    Definition 5.6.4, p. 263 for reduction operator.
    Algorithm 5.6.5, p. 263 for reduction algorithm for indefinite quadratic forms.
    """

    def __init__(self: typing.Self, a: int, b: int, c: int) -> None:
        discriminant = b ** 2 - 4 * a * c
        if not all(isinstance(x, int) for x in [a, b, c]):
            raise TypeError("a, b, c must all be integers")
        if not discriminant > 0:
            raise ValueError(f"{discriminant=} must be positive")
        if c == 0:
            raise ValueError(f"an indefinite binary quadratic form ax²+bx+c must have nonzero c: {c=}")
        self.a = a
        self.b = b
        self.c = c
        self.D = discriminant

    def __repr__(self: typing.Self) -> str:
        return f"{type(self).__name__}({self.a}, {self.b}, {self.c})"

    def __eq__(self: typing.Self, other: typing.Self) -> bool:
        return self.a == other.a and self.b == other.b and self.c == other.c
    
    def __str__(self) -> str:
        return f"{self.a}x²\t{self.b:+}xy\t{self.c:+}y².\tD={self.D}"
    
    def transform(self: typing.Self, matrix: sl2z.SL2Z) -> typing.Self:
        a = self.a * matrix.alpha ** 2 + self.b * matrix.alpha * matrix.gamma + self.c * matrix.gamma ** 2
        b = 2 * self.a * matrix.alpha * matrix.beta + \
            self.b * matrix.alpha * matrix.delta + \
            self.b * matrix.beta * matrix.gamma + \
            2 * self.c * matrix.gamma * matrix.delta
        c = self.a * matrix.beta ** 2 + self.b * matrix.beta * matrix.delta + self.c * matrix.delta ** 2
        return type(self)(a, b, c)
    
    @property
    def gcd(self: typing.Self) -> int:
        return math.gcd(self.a, self.b, self.c)
    
    @property
    def is_primitive(self: typing.Self) -> bool:
        return self.gcd == 1

    @property
    def is_reduced(self: typing.Self) -> bool:
        """
        Cohen, Definition 5.6.2, p. 262: an indefinite binary quadratic form (a,b,c)
        is reduced if and only if
        0 < (-b + √D) / (2|a|) < 1 and (b + √D) / (2|a|) > 1.
        """
        sqrtD = math.sqrt(self.D)
        alpha = (-self.b + sqrtD) / (2 * abs(self.a))
        beta = (-self.b - sqrtD) / (2 * abs(self.a))
        return (0 < alpha < 1) and (beta < -1)
    

    def reduce_with_m(self: typing.Self) -> tuple[typing.Self, int]:
        """
        Cohen, Definition 5.6.4 and Algorithm 5.6.5, p. 263:
        1) Compute t0 = -b mod 2|c| in [0, 2|c|-1]. (Used in Step 2.)
        2) If |c| > sqrt(D), choose r so that -|c| < r ≤ |c|.
           If |c| < sqrt(D), choose r so that sqrt(D) - 2|c| < r < sqrt(D).
           In each case, r ≡ -b (mod 2c).
        3) Record m.
        4) Return BQF (c, r, (r^2 - D) / (4c)) and m.
        """

        _, b, c, D = self.a, self.b, self.c, self.D # self.a not used
        sqrtD = math.sqrt(D)
        abs_c = abs(c)

        # Step 1
        t0 = ((-b) % (2 * abs_c))

        # Step 2
        if abs_c > sqrtD:
            # We need -|c| < r ≤ |c|  AND  r ≡ t0 (mod 2|c|)
            # The interval I = (-abs_c, abs_c] has length 2|c|.
            # Exactly one integer r in I is ≡ t0 (mod 2|c|).
            # r = t0 + 2|c|k for some integer k.
            # We want the unique k such that -abs_c < t0 + 2|c|k ≤ abs_c.
            # k = floor((abs_c - t0) / (2|c|)).
            k = (abs_c - t0) // (2 * abs_c)
            r = t0 + 2 * abs_c * k
        else:
            # We need  (sqrtD - 2|c|) < r < sqrtD,  AND  r ≡ t0 (mod 2|c|)
            # Let lower = sqrtD - 2|c|,  upper = sqrtD.
            lower = sqrtD - 2 * abs_c
            upper = sqrtD
            # We look for integer r = t0 + 2|c|·k lying in (lower, upper).
            # Since upper - lower = 2|c|, there is exactly one such integer.
            # Solve k ≥ (lower - t0)/(2|c|).
            # k = ceil((lower - t0)/(2|c|)).
            k = math.ceil((lower - t0) / (2 * abs_c))
            r = t0 + 2 * abs_c * k

        # Step 3
        m = (b + r) // (2 * c)

        # Step 4
        a1 = c
        b1 = r
        c1 = (r * r - D) // (4 * c)
        return type(self)(a1, b1, c1), m

    
    def reduced_with_m_list(self: typing.Self) -> tuple[typing.Self, list[int]]:
        bqf = self
        m_list = []
        while not bqf.is_reduced:
            bqf, mi = bqf.reduce_with_m()
            m_list.append(mi)
        return (bqf, m_list)

    def real_quadratic_number(self: typing.Self) -> RealQuadraticNumber:
        return RealQuadraticNumber(self.D, Fraction(-self.b, 2 * abs(self.a)), Fraction(1, 2 * abs(self.a)))

    def evaluate(self: typing.Self, x: int, y: int) -> int:
        return self.a * x ** 2 + self.b * x * y + self.c * y **2
    
    @staticmethod
    def m_list_to_word_list(m_list: list[int]) -> list[str]:
        word_list = []
        for m in m_list:
            word = "S" + "T" * m
            word_list.append(word)
        return word_list
    
    @staticmethod
    def word_list_to_word(word_list: list[str]) -> str:
        return "".join(word_list)
    
    @staticmethod
    def word_to_rle_tuple_list(word: str) -> list[tuple[str, int]]:
        return [(letter, len(list(g))) for letter, g in it.groupby(word)]


def class_number(d: int) -> int:
    """
    Return the (narrow) class number of Q(√d) by:
      1) computing the fundamental discriminant Δ,
      2) collecting *all* primitive, reduced forms (a,b,c) of disc Δ,
      3) traversing each Cohen‐reduction cycle exactly once, and
      4) returning the cycle count.

    This is exactly Algorithm 5.6.5 + Section 5.7 in Cohen’s book.
    """

    import math
    from prime_numbers import squarefull_and_squarefree_parts
    from binary_quadratic_forms import IndefiniteBQF

    # (a) Build fundamental discriminant Δ from d:
    _, d0 = squarefull_and_squarefree_parts(d)
    if d0 % 4 == 1:
        Δ = d0
    else:
        Δ = 4 * d0

    # (b) Gather every primitive, reduced form (a,b,c), allowing a<0 or a>0:
    A_max = int(math.isqrt(Δ))
    B     = int(math.isqrt(Δ))

    reduced_forms = set()
    for a in range(-A_max, A_max + 1):
        if a == 0:
            continue
        for b in range(-B, B + 1):
            num   = b*b - Δ
            denom = 4 * a
            if denom == 0 or (num % denom) != 0:
                continue
            c = num // denom

            if math.gcd(a, b, c) != 1:
                continue

            Q = IndefiniteBQF(a, b, c)
            if Q.is_reduced:
                reduced_forms.add((a, b, c))

    # (c) Walk each reduction‐with‐m cycle exactly once; each cycle = one class:
    visited = set()
    class_count = 0

    for (a0, b0, c0) in reduced_forms:
        if (a0, b0, c0) in visited:
            continue

        class_count += 1
        curr = IndefiniteBQF(a0, b0, c0)
        start = (curr.a, curr.b, curr.c)

        while True:
            key = (curr.a, curr.b, curr.c)
            if key in visited:
                break
            visited.add(key)

            next_bqf, _m = curr.reduce_with_m()
            curr = next_bqf

            if (curr.a, curr.b, curr.c) == start:
                # we’ve come back around to the same (a0,b0,c0):
                visited.add(start)
                break

    return class_count


if __name__ == "__main__":
    print(class_number(170))

    d = 48

    assert math.isclose(RealQuadraticField(d).regulator, math.log(float(RealQuadraticField(d).fundamental_unit())))

    assert RealQuadraticField(d).fundamental_unit().norm == 1

    assert float(RealQuadraticField(d).fundamental_unit()) > 1

    m = 13

    d = 4 * m + 1
    _, d0 = prime_numbers.squarefull_and_squarefree_parts(d)
    D = d0
    assert RealQuadraticField(d).D == D

    d = 4 * m + 2
    _, d0 = prime_numbers.squarefull_and_squarefree_parts(d)
    D = 4 * d0
    assert RealQuadraticField(d).D == D

    d = 4 * m + 3
    _, d0 = prime_numbers.squarefull_and_squarefree_parts(d)
    D = 4 * d0
    assert RealQuadraticField(d).D == D

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
    product_matrix = sl2z.SL2Z.word_list_to_matrix(word_list)
    assert bqf.transform(product_matrix) == reduced_bqf

    d = 19
    ε = RealQuadraticField(d).fundamental_unit()
    print(f"Powers of fundamental unit {ε=}")
    for k in range(-10, 10 + 1):
        print(f"ε^{k}:", ε ** k, sep="\t")
