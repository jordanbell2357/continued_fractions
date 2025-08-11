import math
import decimal
from decimal import Decimal
from fractions import Fraction
from numbers import Number
from numbers import Rational
from collections import abc
import itertools as it
import functools as ft
import operator
import typing

import gl2z
import prime_numbers
import pell


class RationalQuadraticPolynomial(typing.NamedTuple):
    a: Rational
    b: Rational
    c: Rational

    def __add__(self: typing.Self, other: typing.Self) -> typing.Self:
        return type(self)(self.a + other.a, self.b + other.b, self.c + other.c)

    def __str__(self: typing.Self) -> str:
        return f"{self.a}x² {self.b:+}x {self.c:+}"

    def evaluate(self: typing.Self, x: Number) -> Number:
        return self.a * x * x + self.b * x + self.c
    
    def is_integral(self: typing.Self) -> bool:
        return all(e == int(e) for e in [self.a, self.b, self.c])


class IntegralQuadraticPolynomial(RationalQuadraticPolynomial):
    def __init__(self: typing.Self, a: int, b: int, c: int) -> None:
        if not all(isinstance(e, int) for e in [a, b, c]):
            raise TypeError(f"{a}, {b}, {c} must be integers.")
        super().__init__(a, b, c)


class RealQuadraticNumber(Number):
    __slots__ = ["d", "x", "y"]

    def __init__(self: typing.Self, d: int, x: Rational, y: Rational) -> None:
        if not isinstance(d, int):
            raise TypeError("d must be integer.")
        if not d > 1:
            raise ValueError("d must be > 1.")
        squarefull_part_d, squarefree_part_d = prime_numbers.squarefull_and_squarefree_parts(d)
        if squarefree_part_d == 1:
            raise ValueError("d must not be a perfect square.")
        self.d = squarefree_part_d
        if isinstance(x, int):
            x = Fraction(x, 1)
        self.x = x
        if isinstance(y, int):
            y = Fraction(y, 1)
        self.y = y * math.isqrt(squarefull_part_d)

    def __repr__(self: typing.Self) -> str:
        return f"{type(self).__name__}({self.d}, {self.x}, {self.y})"
    
    def __eq__(self: typing.Self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        return self.d == other.d and self.x == other.x and self.y == other.y
    
    def __hash__(self: typing.Self) -> int:
        return hash((self.d, self.x, self.y))
    
    def __neg__(self: typing.Self) -> typing.Self:
        return type(self)(self.d, -self.x, -self.y)
    
    def __add__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("Algebraic numbers must be members of same field 𝐐(√d) to be added.")
        return type(self)(self.d, self.x + other.x, self.y + other.y)
    
    def __radd__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("Algebraic numbers must be members of same field 𝐐(√d) to be added.")
        return type(self)(self.d, other.x + self.x, other.y + self.y)
    
    def __sub__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("Algebraic numbers must be members of same field 𝐐(√d) to be subtracted.")
        return type(self)(self.d, self.x - other.x, self.y - other.y)
    
    def __rsub__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("Algebraic numbers must be members of same field 𝐐(√d) to be subtracted.")
        return type(self)(self.d, other.x - self.x, other.y - self.y)
    
    def __mul__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("Algebraic numbers must be members of same field 𝐐(√d) to be multiplied.")
        return type(self)(self.d,
                          self.x * other.x + self.y * other.y * self.d,
                          self.x * other.y + self.y * other.x)
    
    def __rmul__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        if self.d != other.d:
            raise ValueError("Algebraic numbers must be members of same field 𝐐(√d) to be multiplied.")
        return type(self)(self.d,
                          self.x * other.x + self.y * other.y * self.d,
                          self.x * other.y + self.y * other.x)
    
    def __truediv__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        dividend = self
        divisor = other
        if isinstance(divisor, Rational):
            divisor = type(dividend)(dividend.d, divisor, 0)
        if dividend.d != divisor.d:
            raise ValueError("Dividend and divisor must be members of same field 𝐐(√d).")
        divisor_norm = divisor.norm
        if divisor_norm == 0:
            raise ZeroDivisionError(f"Divisor must have nonzero norm: {divisor_norm=}")
        return type(dividend)(dividend.d,
                          Fraction(dividend.x * divisor.x - dividend.y * divisor.y * dividend.d, divisor_norm),
                          Fraction(-dividend.x * divisor.y + dividend.y * divisor.x, divisor_norm))
    
    def __rtruediv__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        dividend = other
        divisor = self
        if isinstance(dividend, Rational):
            dividend = type(divisor)(divisor.d, dividend, 0)
        if dividend.d != divisor.d:
            raise ValueError("Dividend and divisor must be members of same field 𝐐(√d).")
        divisor_norm = divisor.norm
        if divisor_norm == 0:
            raise ZeroDivisionError(f"Divisor must have nonzero norm: {divisor_norm=}")
        return type(dividend)(dividend.d,
                          Fraction(dividend.x * divisor.x - dividend.y * divisor.y * dividend.d, divisor_norm),
                          Fraction(-dividend.x * divisor.y + dividend.y * divisor.x, divisor_norm))
    
    def __pow__(self: typing.Self, exponent: int) -> typing.Self:
        if not isinstance(exponent, int):
            raise TypeError(f"{exponent=} must be an integer.")
        if exponent <= 0 and self.norm == 0:
            raise ZeroDivisionError(f"{exponent=} must be positive for non-invertible elements of 𝐐(√d).")
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
    
    def __lt__(self: typing.Self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        return float(self) < float(other)
    
    def __le__(self: typing.Self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        return float(self) < float(other) or self == other
    
    def __gt__(self: typing.Self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        return float(self) > float(other)
    
    def __ge__(self: typing.Self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(self.d, other, 0)
        return float(self) > float(other) or self == other
    
    def __int__(self: typing.Self) -> int:
        return int(float(self))
    
    def identity(self: typing.Self) -> typing.Self:
        return type(self)(self.d, self.x, self.y)

    def conjugate(self: typing.Self) -> typing.Self:
        return type(self)(self.d, self.x, -self.y)
    
    def __abs__(self: typing.Self) -> float:
        return math.sqrt(self * self.conjugate())
    
    def __str__(self: typing.Self) -> str:
        return f"{self.x} {self.y:+} * √{self.d}"

    @property
    def norm(self: typing.Self) -> Fraction:
        return self.x ** 2 - self.d * self.y ** 2
    
    @property
    def trace(self: typing.Self) -> Fraction:
        return 2 * self.x
    
    @property
    def is_integral(self: typing.Self) -> bool:
        return self.norm == int(self.norm) and self.trace == int(self.trace)
    
    def minimal_polynomial(self: typing.Self) -> RationalQuadraticPolynomial:
        """
        When d is not a perfect square,
        Gal(𝐐(√d) / 𝐐) = {identity, conjugate}
        When d is a perfect square,
        Gal(𝐐(√d) / 𝐐) = {identity}
        
        m_a(x)
        = (x-a)(x-a.conjugate())
        = x ** 2 - (a + a.conjugate()) * x + a * a.conjugate()
        = x ** 2 - 2 * a.x * x + abs(a) ** 2
        """
        return RationalQuadraticPolynomial(1, -2 * self.x, (self * self.conjugate()).x)
    
    def lft_GL2Z(self: typing.Self, matrix: gl2z.GL2Z) -> typing.Self:
        """
        Linear fractional transformation.
        """
        return (matrix.alpha * self + matrix.beta) / (matrix.gamma * self + matrix.delta)
    
    def GL2Q_left_regular_representation(self: typing.Self) -> gl2z.GL2Q:
        """
        Left regular representation of 𝐐(√d) in GL₂(𝐐).
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        § 4.2.3 "The Matrix (or Regular) Representation of an Algebraic Number", p. 160.
        """
        alpha = self.x
        beta = self.y * self.d
        gamma = self.y
        delta = self.x
        return gl2z.GL2Q(alpha, beta, gamma, delta)
    
    @classmethod
    def from_GL2Q_left_regular_representation(cls, matrix: gl2z.GL2Q) -> typing.Self:
        if matrix.gamma == 0:
            raise ValueError(f"{matrix.gamma=} must not be 0.")
        d = Fraction(matrix.beta, matrix.gamma)
        if d != int(d):
            raise ValueError(f"{matrix} must belong to image of left regular representation of 𝐐(√d) in GL₂(𝐐).")
        d = int(d)
        x = matrix.alpha
        y = matrix.gamma
        return cls(d, x, y)
    
    @classmethod
    def discriminant(cls, alpha1: typing.Self, alpha2: typing.Self) -> Rational:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Proposition 4.4.1, p. 165.
        """
        return (alpha1 * alpha2.conjugate() - alpha2 * alpha1.conjugate()) ** 2
    
    @classmethod
    def discriminant_by_trace(cls, alpha1: typing.Self, alpha2: typing.Self) -> Rational:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Proposition 4.4.1, p. 165.
        """
        alpha11 = (alpha1 * alpha1).trace
        alpha12 = (alpha1 * alpha2).trace
        alpha21 = (alpha2 * alpha1).trace
        alpha22 = (alpha2 * alpha2).trace
        return alpha11 * alpha22 - alpha12 * alpha21
        

def eigenvalues(m: gl2z.M2Z) -> tuple[RealQuadraticNumber, RealQuadraticNumber] | tuple[Fraction, Fraction]:
    """
    x^2 - tr(m)x + det(m) = 0
    a = 1, b = -tr(m), c = det(m)
    D = b ** 2 - 4ac
    x1 = -b/2 - sqrt(D)/2
    x2 = -b/2 + sqrt(D)/2
    """
    a = 1
    b = -m.trace
    c = m.det
    D = b ** 2 - 4 * a * c
    if D != math.isqrt(D) ** 2:
        x1 = RealQuadraticNumber(D, Fraction(-b, 2), -Fraction(1, 2))
        x2 = RealQuadraticNumber(D, Fraction(-b, 2), Fraction(1, 2))
        return x1, x2
    else:
        x1 = Fraction(-b - math.isqrt(D), 2)
        x2 = Fraction(-b + math.isqrt(D), 2)
        return x1, x2


class PureQuadraticSurd:
    def __init__(self, d: int, coefficient_fraction: Rational) -> None:
        coefficient_fraction = Fraction(coefficient_fraction)
        d_squarefull, d_squarefree = prime_numbers.squarefull_and_squarefree_parts(d)
        radicand_squarefull_sqrt_fraction = Fraction(math.isqrt(d_squarefull))

        self.d = d_squarefree
        self.coefficient_fraction = coefficient_fraction * radicand_squarefull_sqrt_fraction

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.d}, {self.coefficient_fraction})"
    
    def __eq__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other), Fraction(1))
        return self.d == other.d and self.coefficient_fraction == other.coefficient_fraction
    
    def __neg__(self) -> typing.Self:
        return type(self)(self.d, -self.coefficient_fraction)
    
    def __hash__(self) -> int:
        return hash((self.d, self.coefficient_fraction))
    
    def __str__(self) -> str:
        return f"{self.coefficient_fraction} * √{self.d}"

    @classmethod
    def from_coefficient_fraction_and_radicand_fraction(cls, coefficient_fraction: Rational, radicand_fraction: Rational) -> typing.Self:
        coefficient_fraction, radicand_fraction = Fraction(coefficient_fraction), Fraction(radicand_fraction)
        numerator, denominator = radicand_fraction.numerator, radicand_fraction.denominator
        numerator_squarefull, numerator_squarefree = prime_numbers.squarefull_and_squarefree_parts(numerator)
        denominator_squarefull, denominator_squarefree = prime_numbers.squarefull_and_squarefree_parts(denominator)
        radicand_squarefull_sqrt_fraction = Fraction(math.isqrt(numerator_squarefull), math.isqrt(denominator_squarefull))
        # sqrt(p/q) = sqrt(pq) / sqrt(qq) = sqrt(pq) / q = 1/q * sqrt(pq) = 1/q * d, d = pq
        rationalize_denominator_fraction = Fraction(1, denominator_squarefree)

        d = numerator_squarefree * denominator_squarefree
        coefficient_fraction = coefficient_fraction * radicand_squarefull_sqrt_fraction * rationalize_denominator_fraction
        return cls(d, coefficient_fraction)
    
    @classmethod
    def sqrt_rational(cls, rational_number: Rational) -> typing.Self:
        if isinstance(rational_number, Rational):
            rational_number = Fraction(rational_number)
        return cls.from_coefficient_fraction_and_radicand_fraction(1, rational_number)


    def __mul__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(1, other)
        d1, d2 = self.d, other.d
        d = d1 * d2
        return type(self)(d, self.coefficient_fraction * other.coefficient_fraction)
    
    __rmul__ = __mul__

    def __truediv__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(1, Fraction(other))
        d1, coefficient_fraction1, d2, coefficient_fraction2 = self.d, self.coefficient_fraction, other.d, other.coefficient_fraction
        d_fraction = Fraction(d1, d2)
        coefficient_fraction = Fraction(coefficient_fraction1, coefficient_fraction2)
        return type(self).from_coefficient_fraction_and_radicand_fraction(coefficient_fraction, d_fraction)

    def __rtruediv__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(1, Fraction(other))
        d1, coefficient_fraction1, d2, coefficient_fraction2 = other.d, other.coefficient_fraction, self.d, self.coefficient_fraction
        d_fraction = Fraction(d1, d2)
        coefficient_fraction = Fraction(coefficient_fraction1, coefficient_fraction2)
        return type(self).from_coefficient_fraction_and_radicand_fraction(coefficient_fraction, d_fraction)
    
    def __float__(self) -> float:
        return self.coefficient_fraction * math.sqrt(self.d)
    
    def __lt__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other), Fraction(1))
        return float(self) < float(other)
    
    def __le__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other), Fraction(1))
        return float(self) < float(other) or self == other
    
    def __gt__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other), Fraction(1))
        return float(self) > float(other)
    
    def __ge__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other), Fraction(1))
        return float(self) > float(other) or self == other
    
    def __abs__(self) -> typing.Self:
        return self if self >= 0 else -self
    
    def __pow__(self, exponent: int) -> typing.Self:
        if self == 0 and exponent <= 0:
            return ValueError(f"For non-invertible element {self=}, {exponent=} must be positive integer.")
        if exponent == 0:
            return type(self)(1, 1)
        elif exponent > 0:
            return ft.reduce(operator.mul, exponent * [self], type(self)(1, 1))
        elif exponent < 0:
            return ft.reduce(operator.mul, -exponent * [1 / self], type(self)(1, 1))
        
    def as_real_quadratic_number(self: typing.Self) -> RealQuadraticNumber:
        return RealQuadraticNumber(self.d, 0, self.coefficient_fraction)

    @property
    def norm(self) -> Fraction:
        """
        For α = c√d  (d square–free):

            N(α) =  c²        if d = 1   (α is rational)
            N(α) = –c² d      if d > 1   (two Galois conjugates differ by sign)
        """
        if self.d == 1:                       # rational element
            return self.coefficient_fraction ** 2
        return -(self.coefficient_fraction) ** 2 * self.d
    
    @property
    def trace(self) -> Fraction:
        """
        self + conjugate(self) = 0 unless d = 1
        """
        if self.d == 1:
            return self.coefficient_fraction
        else:
            return 0
    


class RealQuadraticCompositum:
    """
    Element of compositum field K = 𝐐(√p₁, √p₂, …)
    This is an infinite Galois extension. Any particular element belongs to a compositum field
    that is a finite Galois extension.
        rational_part: Rational
        surd_terms : dict[frozenset[int] | int, PureQuadraticSurd]
    """

    def __init__(self,
                 rational_part: Rational = 0,
                 surd_terms: dict[frozenset[int] | int, PureQuadraticSurd] | None = None) -> None:

        self.rational_part: Fraction = Fraction(rational_part)
        self.surd_terms: dict[frozenset[int], PureQuadraticSurd] = {}

        if not surd_terms:
            return

        for key, surd in surd_terms.items():

            # ── 1. coerce the key into the canonical frozenset[int] form ──
            if isinstance(key, int):                       # allow plain prime-product
                radset: frozenset[int] = frozenset({key})
            elif isinstance(key, frozenset) and all(isinstance(p, int) for p in key):
                radset = key                               # already canonical
            else:
                raise TypeError("surd_terms keys must be frozenset[int] or int")

            # ── 2. value-type check (unchanged) ───────────────────────────
            if not isinstance(surd, PureQuadraticSurd):
                raise TypeError("surd_terms values must be PureQuadraticSurd")

            coeff = surd.coefficient_fraction
            if coeff == 0:
                continue

            # ── 3. √1–terms fold into the rational part ───────────────────
            if math.prod(radset) == 1:
                 self.rational_part += coeff
                 continue

            # ── 4. merge with any existing coefficient for the same radset ─
            old   = self.surd_terms.get(radset)
            new_c = coeff + (old.coefficient_fraction if old else 0)
            if new_c:
                self.surd_terms[radset] = PureQuadraticSurd(math.prod(radset), new_c)
            elif radset in self.surd_terms:                # cancellation to zero
                del self.surd_terms[radset]

    def __repr__(self) -> str:
        return (f"{type(self).__name__}(rational_part={self.rational_part!r}, "
                f"surd_terms={self.surd_terms!r})")

    def __str__(self) -> str:
        if not self.surd_terms:
            return str(self.rational_part)
        parts: list[str] = [str(self.rational_part)] if self.rational_part else []
        for rs in sorted(self.surd_terms, key=lambda k: math.prod(k)):
            coeff = self.surd_terms[rs].coefficient_fraction
            sign  = "+" if coeff >= 0 else "-"
            mag   = abs(coeff)
            mag_s = "" if mag == 1 else str(mag) + "·"
            parts.append(f"{sign} {mag_s}√{math.prod(rs)}")
        return " ".join(parts).lstrip("+ ").replace("+ -", "- ")

    def __eq__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(other)
        return (self.rational_part == other.rational_part and
                self.surd_terms   == other.surd_terms)
    
    def __float__(self) -> float:
        return self.rational_part + math.fsum(surd for surd in self.surd_terms.values())
    
    def __lt__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other))
        return float(self) < float(other)
    
    def __le__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other))
        return float(self) < float(other) or self == other
    
    def __gt__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other))
        return float(self) > float(other)
    
    def __ge__(self, other: typing.Self | Rational) -> bool:
        if isinstance(other, Rational):
            other = type(self)(Fraction(other))
        return float(self) > float(other) or self == other
    
    def __abs__(self) -> typing.Self:
        surd_terms = self.surd_terms
        abs_surd_terms = {k: abs(v) for k, v in surd_terms.items()}
        return type(self)(abs(self.rational_part), abs_surd_terms)

    def __neg__(self) -> typing.Self:
        neg = {rs: PureQuadraticSurd(term.d, -term.coefficient_fraction)
               for rs, term in self.surd_terms.items()}
        return type(self)(-self.rational_part, neg)

    def __add__(self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(other)

        new_r = self.rational_part + other.rational_part
        new_s: dict[frozenset[int], PureQuadraticSurd] = {}

        for rs in self.surd_terms.keys() | other.surd_terms.keys():
            c = (self.surd_terms.get(rs, PureQuadraticSurd(math.prod(rs), 0))
                 .coefficient_fraction
               + other.surd_terms.get(rs, PureQuadraticSurd(math.prod(rs), 0))
                 .coefficient_fraction)
            if c:
                new_s[rs] = PureQuadraticSurd(math.prod(rs), c)

        return type(self)(new_r, new_s)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return type(self)(other) - self

    def __mul__(self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            other = type(self)(other)

        new_r = self.rational_part * other.rational_part
        new_s: dict[frozenset[int], Fraction] = {}

        for rs, term in self.surd_terms.items():
            if other.rational_part:
                new_s[rs] = new_s.get(rs, 0) + term.coefficient_fraction * other.rational_part
        for rs, term in other.surd_terms.items():
            if self.rational_part:
                new_s[rs] = new_s.get(rs, 0) + term.coefficient_fraction * self.rational_part

        for rs1, t1 in self.surd_terms.items():
            for rs2, t2 in other.surd_terms.items():
                common = rs1 & rs2
                rs_new = rs1 ^ rs2
                coeff  = t1.coefficient_fraction * t2.coefficient_fraction * math.prod(common)
                if not rs_new:
                    new_r += coeff
                else:
                    new_s[rs_new] = new_s.get(rs_new, 0) + coeff

        sdict = {rs: PureQuadraticSurd(math.prod(rs), c) for rs, c in new_s.items() if c}
        return type(self)(new_r, sdict)

    __rmul__ = __mul__

    def conjugates(self) -> list[typing.Self]:
        if not self.surd_terms:
            return [self]
        primes = sorted({p for rs in self.surd_terms for p in rs})
        conj_list: list[typing.Self] = []
        for signs in it.product([-1, 1], repeat=len(primes)):
            sign_map = dict(zip(primes, signs))
            conj_s: dict[frozenset[int], PureQuadraticSurd] = {}
            for rs, term in self.surd_terms.items():
                sgn = math.prod(sign_map[p] for p in rs)
                coeff = term.coefficient_fraction * sgn
                if coeff:
                    conj_s[rs] = PureQuadraticSurd(term.d, coeff)
            conj_list.append(type(self)(self.rational_part, conj_s))
        return conj_list

    @ft.cached_property
    def trace(self) -> Fraction:
        if not self.surd_terms:
            return self.rational_part
        n = len({p for rs in self.surd_terms for p in rs})
        return self.rational_part * (1 << n)

    @ft.cached_property
    def norm(self) -> Fraction:
        if not self.surd_terms:
            return self.rational_part
        prod = type(self)(1)
        for sigma in self.conjugates():
            prod *= sigma
        if prod.surd_terms:
            raise ArithmeticError("norm didn’t clear surds")
        return prod.rational_part

    def inverse(self) -> typing.Self:
        if self.rational_part == 0 and not self.surd_terms:
            raise ZeroDivisionError("division by zero")
        if not self.surd_terms:
            return type(self)(Fraction(1, self.rational_part))

        adj = type(self)(1)
        for sigma in self.conjugates():
            if sigma != self:
                adj *= sigma

        nrm = self.norm
        inv_r = adj.rational_part / nrm
        inv_s = {rs: PureQuadraticSurd(term.d, term.coefficient_fraction / nrm)
                 for rs, term in adj.surd_terms.items()}
        return type(self)(inv_r, inv_s)

    def __truediv__(self, other):  return self * (type(self)(other) if isinstance(other, Rational) else other).inverse()

    def __rtruediv__(self, other): return type(self)(other) * self.inverse()

    def __pow__(self, exponent: int) -> typing.Self:
        if self == 0 and exponent <= 0:
            return ValueError(f"For non-invertible element {self=}, {exponent=} must be positive integer.")
        if exponent == 0:
            return type(self)(1)
        elif exponent > 0:
            return ft.reduce(operator.mul, exponent * [self], type(self)(1))
        elif exponent < 0:
            return ft.reduce(operator.mul, -exponent * [1 / self], type(self)(1))
        
    @classmethod
    def from_pure_quadratic_surd(cls, surd: PureQuadraticSurd) -> typing.Self:
        return cls(0, {frozenset([surd.d]): surd})
    
    @classmethod
    def sqrt_rational(cls, q: Rational | typing.Self) -> typing.Self:
        if isinstance(q, Rational):
            q = Fraction(q)
        if isinstance(q, cls):
            if q == q.rational_part:
                q = q.rational_part
            else:
                raise TypeError(f"{q=} must be rational.")
        if q == 0:
            return cls(0)
        pure_quadratic_surd = PureQuadraticSurd.sqrt_rational(q)
        compositum_surd = cls.from_pure_quadratic_surd(pure_quadratic_surd)
        return compositum_surd


class RealQuadraticField(abc.Container):
    """
    Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
    Proposition 4.4.1, p. 165, for definition of discriminant using integral basis.
    Proposition 5.1.1, p. 223, for integral basis and discriminant formula for real quadratic fields.

    Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
    Chapter I, Section 6, "Quadratic Number Fields and Their Units", pp. 35-38.
    """

    def __init__(self: typing.Self, d: int) -> None:
        if not isinstance(d, int):
            raise TypeError("d must be integer")
        if not d > 1:
            raise ValueError("d must be > 1")
        _, squarefree_part_d = prime_numbers.squarefull_and_squarefree_parts(d)
        self.d = squarefree_part_d

    def __repr__(self: typing.Self) -> str:
        return f"{type(self).__name__}({d})"

    def __eq__(self: typing.Self, other: typing.Self) -> bool:
        return self.d == other.d
    
    def __hash__(self: typing.Self) -> int:
        return hash(self.d)
    
    def __contains__(self: typing.Self, item: Number) -> bool:
        if isinstance(item, Rational):
            return True
        if isinstance(item, RealQuadraticNumber):
            return self.d == item.d

    @property
    def D(self: typing.Self) -> int: # discriminant
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        p. 167:

        Proposition 4.4.5. The algebraic numbers 𝛼1, ..., 𝛼n form an integral basis
        if and only if they are algebraic integers and if d(𝛼1, ..., 𝛼n) = d(K), where
        d(K) is the discriminant of K.
        """
        return self.d if self.d % 4 == 1 else 4 * self.d
    
    @ft.cached_property
    def fundamental_unit(self: typing.Self) -> RealQuadraticNumber:
        x, y = pell.solve_pell_equation(self.d)
        return RealQuadraticNumber(self.d, x, y)

    @ft.cached_property
    def regulator_float(self: typing.Self) -> float:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Definition 4.9.8, p. 211.
        """
        return math.log(float(self.fundamental_unit))

    def regulator_decimal(self: typing.Self, precision: int = 10) -> Decimal:
        decimal.getcontext().prec = precision
        fundamental_unit = self.fundamental_unit
        x = int(fundamental_unit.x)
        y = int(fundamental_unit.y)
        d = self.d
        return (Decimal(x) + Decimal(y) * Decimal(d).sqrt()).ln()

    def prime_decomposition_type(self, p: int) -> int:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        p. 198: "Theorem 4.8.13 immediately shows how prime numbers decompose in a quadratic field"

        Corollary 4.8.12. The decomposition type of a prime number p in a quadratic
        field K of discriminant D is the following: if (D / p) = -1 then p is inert. If
        (D / p) = 0 then p is ramified (i.e. p𝐙_K = 𝖕²). Finally, if (D / p) = +1, then p splits
        (completely), i.e. p𝐙_K = 𝖕1𝖕2.

        Proposition 5.1.4, p. 224.
        """
        if not isinstance(p, int):
            raise TypeError(f"{p=} must be an integer.")
        if not prime_numbers.isprime(p):
            raise ValueError(f"{p=} must be prime.")
        return prime_numbers.kronecker_symbol(self.D, p)

    @property
    def sqrtd(self: typing.Self) -> RealQuadraticNumber:
        return RealQuadraticNumber(self.d, 0, 1)
    
    @property
    def omega(self):
        # 𝓞_𝐐(√d) = 𝐙[ω]
        d = self.d
        if d % 4 == 1:                # integral basis (1, (1+√d)/2)
            return RealQuadraticNumber(d, Fraction(1, 2), Fraction(1, 2))
        else:                         # basis (1, √d)
            return RealQuadraticNumber(d, 0, 1)
        
    @property
    def omega_order(self: typing.Self) -> RealQuadraticNumber:
        # 𝓞_𝐐(√d) = 𝐙[ω]
        D = self.D
        return RealQuadraticNumber(D, Fraction(D, 2), Fraction(1, 2))


    @property
    def delta(self):
        """
        δ = ω – Tr(ω)/2
        • For d ≡ 1 (mod 4) we have ω = (1+√d)/2          ⇒ δ = (√d–1)/2
        • For d ≡ 2,3 (mod 4) we have ω = √d               ⇒ δ = √d
        𝓞_𝐐(√d) = 𝐙[𝛿]
        """
        d = self.d
        if d % 4 == 1:
            return RealQuadraticNumber(d, Fraction(1, 2), -Fraction(1, 2))
        else:                           # d ≡ 2,3  →  ω = √d,  Tr ω = 0
            return RealQuadraticNumber(d, 0, 1)


    @property
    def integral_basis(self: typing.Self) -> tuple[RealQuadraticNumber, RealQuadraticNumber]:
        b1 = RealQuadraticNumber(self.d, 1, 0)
        b2 = self.omega
        return b1, b2

    def __str__(self: typing.Self) -> str:
        omega = self.omega
        fundamental_unit = self.fundamental_unit
        return f"𝐐(√{self.d}):\tdiscriminant D={self.D}, ring of integers 𝓞_𝐐(√d)=𝐙[{omega}], fundamental unit {fundamental_unit}"
    

def hilbert_theorem_90(d: int, t: Fraction) -> RealQuadraticNumber:
    return RealQuadraticNumber(d, 1, t) / RealQuadraticNumber(d, 1, -t)


def hilbert_theorem_90_inverse(alpha: RealQuadraticNumber) -> Fraction:
    """
    x  =  (1 + dt²) / (1 − dt²)
    y  =  2t      / (1 − dt²)
    """
    x, y = alpha.x, alpha.y
    if x == -1:
        raise ValueError(f"{x=} must not be -1.")
    return y / (x + 1)


class QuadraticIntegerRing(abc.Container):
    """
    Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
    Proposition 4.4.1, p. 165, for definition of discriminant using integral basis.
    Proposition 5.1.1, p. 223, for integral basis and discriminant formula for real quadratic fields.

    Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
    Chapter I, Section 6, "Quadratic Number Fields and Their Units", pp. 35-38.
    """

    __slots__ = ["d"]

    def __init__(self: typing.Self, d: int) -> None:
        if not isinstance(d, int):
            raise TypeError("d must be integer")
        if not d > 1:
            raise ValueError("d must be > 1")
        _, squarefree_part_d = prime_numbers.squarefull_and_squarefree_parts(d)
        self.d = squarefree_part_d

    def __repr__(self: typing.Self) -> str:
        return f"{type(self).__name__}({d})"

    def __eq__(self: typing.Self, other: typing.Self) -> bool:
        return self.d == other.d
    
    def __hash__(self: typing.Self) -> int:
        return hash(self.d)
    
    def __contains__(self: typing.Self, item: RealQuadraticNumber | Rational) -> bool:
        d = self.d
        if isinstance(item, Rational):
            item = RealQuadraticNumber(d, item, 0)
        if not isinstance(item, RealQuadraticNumber):
            raise TypeError(f"{item=} must be RealQuadraticNumber or Rational.")
        if item.d != self.d:
            return False
        return item.is_integral


if __name__ == "__main__":
    d = 67

    if d % 4 == 1:
        b1, b2 = (RealQuadraticNumber(d, 1, 0), RealQuadraticNumber(d, Fraction(1, 2), Fraction(1, 2)))
    elif d % 4 in [2, 3]:
        b1, b2 = (RealQuadraticNumber(d, 1, 0), RealQuadraticNumber(d, 0, 1))
    assert b1.is_integral and b2.is_integral
    assert b1.minimal_polynomial().is_integral and b2.minimal_polynomial().is_integral
    determinant = b1 * b2.conjugate() - b2 * b1.conjugate()
    discriminant = (determinant ** 2).x
    assert discriminant == RealQuadraticField(d).D

    d = 48

    assert math.isclose(RealQuadraticField(d).regulator_float, math.log(float(RealQuadraticField(d).fundamental_unit)))
    assert math.isclose(RealQuadraticField(d).regulator_decimal(), math.log(float(RealQuadraticField(d).fundamental_unit)))

    assert RealQuadraticField(d).fundamental_unit.norm == 1

    assert float(RealQuadraticField(d).fundamental_unit) > 1

    d = 43
    x = 2
    y = Fraction(2, 7)
    assert RealQuadraticNumber(d, x, y) in RealQuadraticField(d)

    d = 43
    x = 2
    y = Fraction(2, 7)
    assert RealQuadraticNumber(d, x, y).is_integral and RealQuadraticNumber(d, x, y).minimal_polynomial().is_integral() or \
        (not RealQuadraticNumber(d, x, y).is_integral and not RealQuadraticNumber(d, x, y).minimal_polynomial().is_integral())

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

    real_quadratic_number = RealQuadraticNumber(d, Fraction(11, 2), Fraction(1, 25))
    m = real_quadratic_number.GL2Q_left_regular_representation()
    assert m.det == real_quadratic_number.norm
    assert m.trace == real_quadratic_number.trace
    assert RealQuadraticNumber.from_GL2Q_left_regular_representation(m) == real_quadratic_number

    d = 15
    source_number = RealQuadraticNumber(d, Fraction(11, 2), Fraction(1, 25))
    target_number = source_number.lft_GL2Z(gl2z.P)
    divided_number = target_number / source_number
    m1 = divided_number.GL2Q_left_regular_representation()
    m2 = target_number.GL2Q_left_regular_representation() / source_number.GL2Q_left_regular_representation()
    assert m1 == m2

    d = 13
    number = RealQuadraticNumber(d, Fraction(11, 2), Fraction(1, 25))
    matrix = gl2z.S
    transformed_number = number.lft_GL2Z(matrix)
    image = number.GL2Q_left_regular_representation()
    transformed_image = image.action_GL2Z_on_real_quadratic_field_image(matrix)
    preimage_transformed_image = RealQuadraticNumber.from_GL2Q_left_regular_representation(transformed_image)
    assert preimage_transformed_image == transformed_number

    d = 20
    power_check_height = 10
    fundamental_unit_sqrtd = RealQuadraticField(d).fundamental_unit
    assert all(fundamental_unit_sqrtd ** (-k) == (fundamental_unit_sqrtd ** k).conjugate() for k in range(power_check_height))

    d = 17
    omega_d = RealQuadraticField(d).omega
    assert RealQuadraticNumber.discriminant(RealQuadraticNumber(d, 1, 0), omega_d) == RealQuadraticField(d).D

    d = 30
    assert RealQuadraticNumber.discriminant(RealQuadraticNumber(d, 1, 0), RealQuadraticField(d).fundamental_unit) == \
        RealQuadraticNumber.discriminant_by_trace(RealQuadraticNumber(d, 1, 0), RealQuadraticField(d).fundamental_unit)
    
    d = 26
    t = Fraction(-2, 7) # t != 1
    x = hilbert_theorem_90(d, t)
    assert x.norm == 1
    assert hilbert_theorem_90_inverse(x) == t

    d = 26
    epsilon = RealQuadraticField(d).fundamental_unit
    t = hilbert_theorem_90_inverse(epsilon)
    assert hilbert_theorem_90(d, t) == epsilon

    m = gl2z.M2Z(4, 1, 1, 1)
    eigenvalues = eigenvalues(m)
    assert m.trace == sum(eigenvalues)
    assert m.det == math.prod(eigenvalues)

    d = 34
    K = RealQuadraticField(d)
    O_K = QuadraticIntegerRing(d)
    omega_d = K.omega
    sqrtd = K.sqrtd
    epsilon_d = K.fundamental_unit
    assert omega_d in O_K
    assert sqrtd in O_K
    assert epsilon_d in O_K

    # 1. square-factor is pulled out:  √8  →  2√2  so coefficient triples (3*2)
    s = PureQuadraticSurd(8, 3)
    assert s.d == 2
    assert s.coefficient_fraction == 6

    # 2. value-equality (same d & coefficient) but different Python objects
    assert PureQuadraticSurd(2, Fraction(3, 2)) == PureQuadraticSurd(2, Fraction(3, 2))

    # 3. multiplication by rational keeps d and rescales coefficient
    t = PureQuadraticSurd(2, Fraction(3, 2)) * 2
    assert t == PureQuadraticSurd(2, 3)

    # 4. mixing radicands multiplies them:  (3/2)√2 · 5√3  = (15/2)√6
    u = PureQuadraticSurd(2, Fraction(3, 2)) * PureQuadraticSurd(3, 5)
    assert u == PureQuadraticSurd(6, Fraction(15, 2))

    u = PureQuadraticSurd(2, Fraction(3, 2)) * PureQuadraticSurd(3, 5)
    assert u * u ** (-1) == 1
    assert u ** 5 * u ** (-5) == 1
    assert u.norm ** 2 == (u ** 2).norm

    x1 = RealQuadraticCompositum(1, {frozenset([2]): PureQuadraticSurd(2, 5)})
    x2 = RealQuadraticCompositum(1, {frozenset([3]): PureQuadraticSurd(3, 5)})
    x1 * x2
    assert (x2 * x2.inverse()).surd_terms == {}       # product is rational
    assert (x1 / x2) * x2 == x1                       # division round-trip

    a = RealQuadraticCompositum(0, {frozenset([2]): PureQuadraticSurd(2, 3)})   # 3√2
    b = RealQuadraticCompositum(0, {frozenset([2]): PureQuadraticSurd(2, 4)})   # 4√2
    c = RealQuadraticCompositum(1, {frozenset([2]): PureQuadraticSurd(2, 3)})   # 1 + 3√2
    d = RealQuadraticCompositum(1, {frozenset([3]): PureQuadraticSurd(3, 5)})   # 1 + 5√3

    # 1. surds with same d merge under addition
    assert (a + b) == RealQuadraticCompositum(0, {frozenset([2]): PureQuadraticSurd(2, 7)})

    # 2. (3√2)(4√2) = 24  – product collapses to rational
    assert (a * b) == RealQuadraticCompositum(24)

    # 3. trace  of 1 + 3√2   is   2  (because two conjugates: ±√2)
    assert c.trace == 2

    # 4. norm  of 1 + 3√2   is   1 − 18 = −17
    assert c.norm == -17

    # 5. inverse round-trip
    assert c * c.inverse() == RealQuadraticCompositum(1)

    # 6. division round-trip (non-trivial radicands 2 & 3)
    assert (c / d) * d == c

    # 1. Adding rational numbers (within the field):
    assert RealQuadraticCompositum(5) + RealQuadraticCompositum(2) == RealQuadraticCompositum(7)  # 5 + 2 = 7

    # 2. Rational promotion and surd addition:
    sqrt2 = RealQuadraticCompositum(0, {2: PureQuadraticSurd(2, 1)})

    assert sqrt2 + 2 == RealQuadraticCompositum(2, {2: PureQuadraticSurd(2, 1)})  # 2 + √2
    # 3. Multiplying surds (√2 * √2 = 2, a rational):
    assert sqrt2 * sqrt2 == RealQuadraticCompositum(2)  # √2 * √2 = 2

    # 5. Galois conjugates, trace, and norm in Q(√2, √3):
    elem = RealQuadraticCompositum(1, {2: PureQuadraticSurd(2, 3), 3: PureQuadraticSurd(3, 1)})  # 1 + 3√2 + 1√3
    conjugates = elem.conjugates()
    # There should be 4 conjugates for two independent primes:
    assert len(conjugates) == 4
    # Sum of conjugates equals 2^2 * (rational_part) = 4:
    total = RealQuadraticCompositum(0)
    for conj in conjugates:
        total += conj
    assert total == RealQuadraticCompositum(4)           # trace = 4
    # Product of conjugates is rational (norm):
    prod = RealQuadraticCompositum(1)
    for conj in conjugates:
        prod *= conj
    assert prod.surd_terms == {}                        # all surds cancel out
    assert prod.rational_part == elem.norm            # norm = 184 in this case

    a  = RealQuadraticCompositum(1, {2: PureQuadraticSurd(2, 1)})      # 1 + √2
    b  = RealQuadraticCompositum(0, {3: PureQuadraticSurd(3, 2)})  # 2√3
    assert a.trace == 2                 # two conjugates (±√2)
    assert b.norm  == -12               # N(2√3)=−4·3
    assert (a * b).inverse() * (a * b) == RealQuadraticCompositum(1)
