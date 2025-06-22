import math
import decimal
from decimal import Decimal
from fractions import Fraction
from numbers import Number
from numbers import Rational
from collections import abc
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
    
    def lft_GL2Z(self: typing.Self, matrix: gl2z.GL2Z) -> typing.Self: # linear fractional transformation
        return (matrix.alpha * self + matrix.beta) / (matrix.gamma * self + matrix.delta)
    
    def GL2Q_left_regular_representation(self: typing.Self) -> gl2z.GL2Q: # Left regular representation of 𝐐(√d) in GL₂(𝐐)
        """
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

    def GL2Q_integral_basis_representation(self: typing.Self) -> gl2z.GL2Q:
        d = self.d
        x = self.x
        y = self.y

        if d % 4 in [2, 3]:
            """
            α = x + y√d
            ω = √d
            α = x + yω
            α * 1 = x * 1 + y * ω
            α * ω  = (x + yω) * ω = x√d + y * d = (d * y) + x√d = (d * y) * 1 + x * ω

            [ x      y * d ]
            [ y      x     ]
            """
            return gl2z.GL2Q(x, d * y, y, x)
        elif d % 4 == 1:
            """
            α = x + y√d
            ω = (1 + √d) / 2
            √d = 2ω − 1
            
            √d * ω  
            = √d * (1 + √d)/2  
            = (√d + d)/2  
            = d/2 + (1/2) * √d  
            = d/2 + (1/2) * (2ω − 1)  
            = d/2 + ω − 1/2  
            = (d − 1)/2 + ω 

            α * 1  
            = x + y * √d  
            = x + y(2ω − 1)  
            = (x − y) + (2y) * ω  
            = (x − y) * 1 + (2y) * ω  

            α * ω  
            = (x + y√d) * ω  
            = xω + y(√d * ω)
            = xω + y((d − 1)/2 + ω)
            = (d − 1)/2 * y * 1 + (x + y) * ω

            [ x - y      (d-1)/2 * y ]
            [ 2y         x + y       ]
            """
            return gl2z.GL2Q(x - y, Fraction(d - 1, 2) * y, 2 * y, x + y)
        else:
            raise ValueError(f"{d=} must not be a multiple of 4.")
        

    @classmethod
    def from_GL2Q_integral_basis_representation(cls, matrix: gl2z.GL2Q) -> typing.Self:
        if matrix.gamma == 0:
            raise ValueError(f"{matrix.gamma=} must not be 0.")
        x = Fraction(matrix.alpha + matrix.delta, 2)
        y = Fraction(matrix.delta - matrix.alpha, 2)
        d = Fraction(matrix.beta, y) * 2 + 1
        if d != int(d):
            raise ValueError(f"{matrix} must belong to image of integral basis representation of 𝐐(√d) in GL₂(𝐐).")
        d = int(d)
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
    
    @property
    def integer_coefficient_tuple(self: typing.Self) -> tuple[int,int]:
        d = self.d
        if not self.is_integral:
            raise TypeError(f"{self} must be an integral element of 𝐐(√d).")
        if d % 4 in (2, 3):       # ω = √d
            return (int(self.x), int(self.y))
        else:                    # d % 4 == 1,  ω = (1+√d)/2
            return (int(self.x - self.y), int(2 * self.y))
        

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


class RealQuadraticField(abc.Container):
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
    
    @property
    def fundamental_unit(self: typing.Self) -> RealQuadraticNumber:
        x, y = pell.solve_pell_equation(self.d)
        return RealQuadraticNumber(self.d, x, y)

    @property
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
    def omega(self: typing.Self) -> RealQuadraticNumber: # 𝓞_𝐐(√d) = 𝐙[ω]
        d = self.d
        if d % 4 == 1:
            return RealQuadraticNumber(d, Fraction(1, 2), Fraction(1, 2))
        elif d % 4 in [2, 3]:
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


class NonzeroIdeal(abc.Container):
    """
    Nonzero ideal in the ring of integers 𝓞_K of K=𝐐(√d).

    Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.

    p. 166:
    Theorem 4.4.2. The ring 𝐙_K is a free Z-module of rank n=deg(K). This is true more generally for any non-zero ideal of
    𝐙_K.

    p. 166:
    Definition 4.4.3. A 𝐙-basis of the free module 𝐙_K will be called an integral
    basis of K. The discriminant of an integral basis is independent of the choice
    of that basis, and is called the discriminant of the field K and is denoted by
    d(K).

    Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
    Chapter I, Section 7, "Relationship of Quadratic Forms to Ideals", pp. 38-50.

    p. 44:
    Lemma 1.19.
    (a) If a ≠ 0 and b' are integers such that a divides N(b' + δ) in 𝐙, then
    (a, b' + δ) = ⟨a, b' + δ⟩ in the sense that the free abelian subgroup of R generated
    by a and b' + δ coincides with the ideal generated by a and b' + δ.
    (b) If I is any nonzero ideal in R, then I is of the form I = ⟨a, r⟩ for some
    integer a > 0 and some r in R.
    """

    __slots__ = ["a", "r"]


    @staticmethod
    def orientation(r1: RealQuadraticNumber, r2: RealQuadraticNumber) -> RealQuadraticNumber:
        if not (r1.is_integral and r2.is_integral):
            raise ValueError(f"{r1=} and {r2=} must be integral elements of 𝐐(√d).")
        if r1.d != r2.d:
            raise ValueError(f"{r1=} and {r2=} must belong to same ring of integers 𝓞_𝐐(√d).")
        return r1 * r2.conjugate() - r1.conjugate() * r2

    def __init__(self, a: RealQuadraticNumber | int, r: RealQuadraticNumber) -> None:
        if not isinstance(r, RealQuadraticNumber):
            raise TypeError(f"{r=} must be in 𝐐(√d).")
        d = r.d
        if isinstance(a, int):
            a = RealQuadraticNumber(d, a, 0)
        if not isinstance(a, RealQuadraticNumber):
            raise TypeError(f"{a=} must be in 𝐐(√d).")
        if not (a.is_integral and r.is_integral):
            raise ValueError(f"{a=} and {r=} must be integral elements of 𝐐(√d).")
        if a.d != d or r.d != d:
            raise ValueError("Generators must belong to the same ring of integers 𝓞_𝐐(√d).")
        oriented_volume = type(self).orientation(a, r)
        if oriented_volume == 0:
            raise ValueError("A nonzero ideal in the ring of integers of 𝐐(√d) has rank 2.")

        # ----- clear denominators so that we can work with an integer 2×2 matrix
        # coordinates of a, r in the naïve {1, √d} basis
        x1, y1 = a.x, a.y  # Fractions
        x2, y2 = r.x, r.y
        den_lcm = math.lcm(x1.denominator, y1.denominator, x2.denominator, y2.denominator)
        X1, Y1 = int(x1 * den_lcm), int(y1 * den_lcm)
        X2, Y2 = int(x2 * den_lcm), int(y2 * den_lcm)

        # ----- column‑HNF on the integer matrix [[X1,X2],[Y1,Y2]]
        M = gl2z.M2Z(X1, X2,
                     Y1, Y2)
        # We need H = M·V with V ∈ GL₂(ℤ).  Obtain V via row‑HNF on Mᵗ.
        Mt = gl2z.M2Z(X1, Y1,
                      X2, Y2)
        U, Ht = gl2z.hnf_2x2(Mt)        # Ht = U·Mt  (row‑HNF)
        V = U.transpose()               # V = Uᵗ  ∈ GL₂(ℤ)
        H = M * V                       # column‑HNF

        # columns of H give the *scaled* generators
        a_scaled = RealQuadraticNumber(d, Fraction(H.a11, den_lcm), Fraction(H.a21, den_lcm))
        r_scaled = RealQuadraticNumber(d, Fraction(H.a12, den_lcm), Fraction(H.a22, den_lcm))

        # ----- make a positive
        if a_scaled < RealQuadraticNumber(d, 0, 0):
            a_scaled = -a_scaled
            r_scaled = -r_scaled

        # ----- reduce r modulo a so that 0 ≤ r.y < a
        # r' ≡ r  (mod a)  with y‑coordinate in that interval
        q = math.floor(float(r_scaled.y / a_scaled.y)) if a_scaled.y != 0 else 0
        r_reduced = r_scaled - q * a_scaled
        while r_reduced.y < 0:
            r_reduced += a_scaled
        while r_reduced.y >= a_scaled.y and a_scaled.y != 0:
            r_reduced -= a_scaled

        self.a, self.r = a_scaled, r_reduced

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.a}, {self.r})"
    
    def __eq__(self, other: typing.Self) -> bool:
        return self.a == other.a and self.r == other.r
    
    @property
    def d(self) -> int:
        return self.r.d
    
    @property
    def oriented_volume(self) -> RealQuadraticNumber:
        return type(self).orientation(self.a, self.r)
    
    @property
    def volume(self) -> RealQuadraticNumber:
        if self.oriented_volume < 0:
            return -self.oriented_volume
        return self.oriented_volume

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

    @property
    def norm(self) -> int:
        d = self.d
        if d % 4 == 1: # D = d
            D_sqrt = RealQuadraticNumber(d, 0, 1)
            return int(self.volume / D_sqrt)
        elif d % 4 in [2, 3]: # D = 4d
            D_sqrt = RealQuadraticNumber(d, 0, 2)
            return int(self.volume / D_sqrt)

    def __contains__(self, item: RealQuadraticNumber | Rational) -> bool:
        """
        a1 * r1 + a2 * r2 = r
        a1(x1 + y1√d) + a2(x2 + y2√d) = x + y√d
        a1x1 + a2x2 + (a1y1 + a2y2)√d = x + y√d
        x1a1 + x2a2 = x
        y1a1 + y2a2 = y

        Cramer's rule:
        a1 = (x * y2 - x2 * y) / (x1y2 - x2y1) 
        a2 = (x1 * y - x * y1) / (x1y2 - x2y1)
        """
        if isinstance(item, Rational):
            item = RealQuadraticNumber(self.d, item, 0)
        if not item.is_integral:
            raise ValueError(f"{item=} must be an integral element of 𝐐(√d).")
        a, r = self.a, self.r
        ax, ay, rx, ry, x, y = a.x, a.y, r.x, r.y, item.x, item.y
        D = ax * ry - rx * ay # nonzero by nonzero orientation

        a1 = Fraction(x * ry - rx * y, D)
        a2 = Fraction(ax * y - x * ay, D)
        if a1.denominator == 1 and a2.denominator == 1:
            return True
        return False
    
    def __le__(self, other: typing.Self) -> typing.Self:
        if self.a in other and self.r in other:
            return True
        return False
    
    def __lt__(self, other: typing.Self) -> typing.Self:
        if self.a in other and self.r in other and self != other:
            return True
        return False
    
    def __ge__(self, other: typing.Self) -> typing.Self:
        if other.a in self and other.r in self:
            return True
        return False
    
    def __gt__(self, other: typing.Self) -> typing.Self:
        if other.a in self and other.r in self and self != other:
            return True
        return False


    def __mul__(self: typing.Self, other: typing.Self | Rational) -> typing.Self:
        if isinstance(other, Rational):
            # Multiply ideal by a rational (integer/fraction) – scale generators
            return type(self)(self.a * other, self.r * other)
        if not isinstance(other, type(self)):
            return NotImplemented
        if self.d != other.d:
            raise ValueError("Ideals must belong to the same quadratic field to be multiplied.")
        # Extract coordinates of generators in the naive {1, √d} basis
        d = self.d
        x1, y1 = self.a.x, self.a.y
        x2, y2 = self.r.x, self.r.y
        x3, y3 = other.a.x, other.a.y
        x4, y4 = other.r.x, other.r.y
        # Compute common denominator to clear all fractions
        den_lcm = math.lcm(x1.denominator, y1.denominator, 
                            x2.denominator, y2.denominator,
                            x3.denominator, y3.denominator,
                            x4.denominator, y4.denominator)
        # Convert all coordinates to integers
        X1, Y1 = int(x1 * den_lcm), int(y1 * den_lcm)
        X2, Y2 = int(x2 * den_lcm), int(y2 * den_lcm)
        X3, Y3 = int(x3 * den_lcm), int(y3 * den_lcm)
        X4, Y4 = int(x4 * den_lcm), int(y4 * den_lcm)
        # Form the 2×4 integer matrix of combined generators
        M = gl2z.M2x4Z(X1, X2, X3, X4,
                       Y1, Y2, Y3, Y4)
        # Compute column Hermite normal form H = M * V  (rank ≤ 2)
        U, H = gl2z.hnf_2x4(M)
        # Identify pivot columns (with nonzero entries)
        c1 = next(c for c in range(1, 5) if H.entry(1, c) != 0)
        c2 = next(c for c in range(c1 + 1, 5) if H.entry(2, c) != 0)
        # Scaled ideal generators from HNF columns (divide back by den_lcm)
        a_scaled = RealQuadraticNumber(d, Fraction(H.entry(1, c1), den_lcm),
                                          Fraction(H.entry(2, c1), den_lcm))
        r_scaled = RealQuadraticNumber(d, Fraction(H.entry(1, c2), den_lcm),
                                          Fraction(H.entry(2, c2), den_lcm))
        # Ensure a_scaled is positive (choose canonical orientation)
        if a_scaled < RealQuadraticNumber(d, 0, 0):
            a_scaled = -a_scaled
            r_scaled = -r_scaled
        # Reduce r_scaled modulo a_scaled (adjust second generator)
        if a_scaled.y != 0:
            q = math.floor(float(r_scaled.y / a_scaled.y))
        else:
            q = 0
        r_reduced = r_scaled - q * a_scaled
        while r_reduced.y < 0:
            r_reduced += a_scaled
        while a_scaled.y != 0 and r_reduced.y >= a_scaled.y:
            r_reduced -= a_scaled
        # If product ideal is principal, simplify generators to standard form
        norm_product = self.norm * other.norm
        p = math.isqrt(norm_product)
        if p * p == norm_product:
            K = RealQuadraticField(d)
            a_final = RealQuadraticNumber(d, p, 0)
            r_final = p * K.omega
            return type(self)(a_final, r_final)
        # Otherwise, return the ideal with the computed generators
        return type(self)(a_scaled, r_reduced)



    @classmethod
    def prime_ideal(cls, d: int, p: int, t_sgn: int = 1) -> typing.Self:
        """
        Henri Cohen, A Course in Computational Algebraic Number Theory, Prop. 5.1.4, p. 224.

        Let K = Q(√D), 𝓞_K = Z[ω], with ω = (D + √D)/2.

        (1) If (D/p)=0 then p ramifies and
            pZ_K = 𝔭²,  𝔭 = ⟨p, ω⟩  (unless p=2 and D≡12 mod 16, then 𝔭=⟨2,1+ω⟩).
        (2) If (D/p)=–1 then p is inert and
            pZ_K itself is prime,  i.e.  (p) = ⟨p, p·ω⟩  (norm p²).
        (3) If (D/p)=+1 then p splits and
            pZ_K = 𝔭₁·𝔭₂,  
            𝔭₁ = ⟨p, ω – (D + b)/2⟩,  
            𝔭₂ = ⟨p, ω – (D – b)/2⟩,  
            where b² ≡ D (mod 4p) and you pick the sign of b via t_sgn=±1.
        """
        K = RealQuadraticField(d)
        D = K.D
        typ = K.prime_decomposition_type(p)

        # (1) ramified
        if typ == 0:
            if p == 2 and D % 16 == 12:
                return cls(p, 1 + K.omega)
            return cls(p, K.omega)
        # (2) inert: principal ideal (p) of norm p^2
        elif typ == -1:
            return cls(p, p * K.omega)
        # (3) split
        # solve b^2 ≡ D (mod 4p), pick ± via t_sgn
        elif typ == 1:
            b = prime_numbers.solve_quadratic_congruence(D, 4 * p)
            r = K.omega - Fraction(D + t_sgn * b, 2)
            return cls(p, r)
        else:
            raise ArithmeticError("Supporting ingredients are broken.")


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

    d = 13
    number = RealQuadraticNumber(d, Fraction(11, 2), Fraction(1, 25))
    mq = number.GL2Q_integral_basis_representation()
    assert RealQuadraticNumber.from_GL2Q_integral_basis_representation(mq) == number

    d = 19
    u = RealQuadraticNumber(d, 2, 3)
    assert u.norm == u.GL2Q_integral_basis_representation().det
    assert u.trace == u.GL2Q_integral_basis_representation().trace

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
    t = Fraction(-2, 7)
    x = hilbert_theorem_90(d, t)
    assert x.norm == 1
    assert hilbert_theorem_90_inverse(x) == t

    # Ideals

    d = 13
    a = RealQuadraticNumber(d, 1, 0) # 1
    r = RealQuadraticNumber(d, 0, 1) # √d
    ideal1 = NonzeroIdeal(a, r)
    ideal2 = NonzeroIdeal(r, a)
    assert ideal1 == ideal2

    d = 13
    u = RealQuadraticNumber(d, 1, 0) # 1
    v = RealQuadraticNumber(d, 0, 1) # √d
    ideal1 = NonzeroIdeal(u, v)
    ideal2 = NonzeroIdeal(v, u)
    assert ideal1 == ideal2
    assert 3 * u + 5 * v in ideal1

    # unit 1 not in ⟨2, v⟩
    d = 13
    u = RealQuadraticNumber(d, 1, 0) # 1
    v = RealQuadraticNumber(d, 0, 1) # √d
    ideal = NonzeroIdeal(2 * u, v)
    assert u not in ideal

    d = 17
    omega_d = RealQuadraticField(d).omega
    assert NonzeroIdeal(RealQuadraticNumber(d, 1, 0), omega_d).norm == 1

    # ⟨2, √13⟩ is a proper sub-ideal of ⟨1, √13⟩
    d  = 13
    u  = RealQuadraticNumber(d, 1, 0)   # 1
    v  = RealQuadraticNumber(d, 0, 1)   # √13
    ideal1  = NonzeroIdeal(2*u, v)      # ⟨2, √13⟩
    ideal2  = NonzeroIdeal(u,   v)
    assert ideal1 < ideal2
    assert ideal2 > ideal1
    assert not (ideal2 < ideal1)

    m = gl2z.M2Z(4, 1, 1, 1)
    eigenvalues = eigenvalues(m)
    assert m.trace == sum(eigenvalues)
    assert m.det == math.prod(eigenvalues)

    d = 13
    u = RealQuadraticNumber(d, 1, 0)    # 1
    v = RealQuadraticNumber(d, 0, 1)    # √13
    ideal1 = NonzeroIdeal(2 * u, v)         # ⟨2, √13⟩
    ideal2 = NonzeroIdeal(u, u + v)         # ⟨1, 1 + √13⟩
    ideal_product = ideal1 * ideal2
    print(ideal1, ideal2, ideal_product)
    assert ideal1.norm * ideal2.norm == ideal_product.norm
    assert ideal_product <= ideal1 and ideal_product <= ideal2

    # ---------------------------------------------------------------------------
    # ❶  Split prime – Q(√5), p = 11
    #     p𝓞_K = 𝔭₁ · 𝔭₂   with 𝔭₁ ≠ 𝔭₂
    # ---------------------------------------------------------------------------
    d, p = 5, 11                          # 5 is quadratic residue mod 11
    𝔭1 = NonzeroIdeal.prime_ideal(d, p,  1)   # ⟨p, p + t√d⟩
    𝔭2 = NonzeroIdeal.prime_ideal(d, p, -1)   # conjugate prime
    assert 𝔭1.norm == 𝔭2.norm == p

    # π = 𝔭1 * 𝔭2                            # should be the (principal) ideal (p)
    # assert π.norm == p ** 2                # N(𝔭₁𝔭₂) = N(𝔭₁)·N(𝔭₂)
    # assert π <= 𝔭1 and π <= 𝔭2            # product ideal is contained in each factor
    # assert RealQuadraticNumber(d, p, 0) in π   # the integer p itself lies in (p)


    # # ---------------------------------------------------------------------------
    # # ❷  Ramified prime – Q(√5), p = 5
    # #     p𝓞_K = 𝔭²
    # # ---------------------------------------------------------------------------
    # d, p = 5, 5
    # 𝔭  = NonzeroIdeal.prime_ideal(d, p)   # ⟨p, ω⟩  (ramified)

    # π  = 𝔭 * 𝔭                            # (p)  principal

    # assert 𝔭.norm == p                    # N(𝔭) = p
    # assert π.norm == p ** 2               # N(𝔭²) = p²
    # assert π <= 𝔭                         # 𝔭² ⊂ 𝔭
    # assert RealQuadraticNumber(d, p, 0) in π   # p ∈ (p)


    # ---------------------------------------------------------------------------
    # ❸  Inert prime – Q(√5), p = 3
    #     single prime ideal of norm p²
    # ---------------------------------------------------------------------------
    # d, p = 5, 3
    # 𝔮   = NonzeroIdeal.prime_ideal(d, p)  # ⟨p, √d⟩  (inert)

    # assert 𝔮.norm == p ** 2               # inert ⇒ N(𝔮) = p²
    # assert (𝔮 * 𝔮).norm == (p ** 2) ** 2  # norm multiplicativity
    # assert (𝔮 * 𝔮) <= 𝔮                  # product is contained in its factor


    # ---------------------------------------------------------------------------
    # ❹  Conjugate product in a different field – Q(√13), p = 17 (splits)
    # ---------------------------------------------------------------------------
    # d, p = 13, 17                         # 13 is a quadratic residue mod 17
    # 𝔭_plus  = NonzeroIdeal.prime_ideal(d, p,  1)
    # 𝔭_minus = NonzeroIdeal.prime_ideal(d, p, -1)

    # π = 𝔭_plus * 𝔭_minus                  # should be (p)

    # assert π.norm == p ** 2
    # assert RealQuadraticNumber(d, p, 0) in π
    # assert π <= 𝔭_plus and π <= 𝔭_minus
