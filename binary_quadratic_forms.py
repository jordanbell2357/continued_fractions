import math
from fractions import Fraction
from numbers import Rational
from collections import abc
import typing

import gl2z
import quadratic_fields


class IndefiniteBQF(abc.Hashable):
    """
    Primitive indefinite binary quadratic forms.

    Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
    Definition 5.2.3, p. 225, for binary quadratic forms.
    Definition 5.6.2, p. 262, for reduced indefinite binary quadratic forms.
    Section 5.6, p. 263, for (a, b, c) being reduced if and only if 0 < (-b + √D) / (2|a|) < 1 and (b + √D) / (2|a|) > 1.
    Definition 5.6.4, p. 263 for reduction operator and helper function r.
    Algorithm 5.6.5, p. 263 for reduction algorithm for indefinite quadratic forms.

    Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
    Chapter I, Sections 3, "Equivalence and Reduction of Quadratic Forms", pp. 12-24.
    Chapter I, Section 4, "Composition of Forms, Class Group", pp. 24-31.
    Chapter I, Section 5, "Genera", pp. 31-34.
    Chapter I, Section 7, "Relationship of Quadratic Forms to Ideals", pp. 38-50.
    """

    __slots__ = ["a", "b", "c"]

    def __init__(self: typing.Self, a: int, b: int, c: int) -> None:
        if not all(isinstance(x, int) for x in [a, b, c]):
            raise TypeError("a, b, c must all be integers.")
        if a==0 or c == 0:
            raise ValueError(f"For indefinite binary quadratic form ax²+bx+c, a and c must be nonzero: {a=}, {c=}.")
        if math.gcd(a, b, c) > 1:
            raise ValueError(f"Primitive requires gcd(a, b, c) = 1: {math.gcd(a, b, c)=}")
        D = b ** 2 - 4 * a * c
        if D <= 0:
            raise ValueError(f"Discriminant {D=} must be positive.")
        if math.sqrt(D) == math.isqrt(D):
            raise ValueError(f"Discriminant {D=} of must not be a perfect square.")
        self.a = a
        self.b = b
        self.c = c

    def __repr__(self: typing.Self) -> str:
        return f"{type(self).__name__}({self.a}, {self.b}, {self.c})"

    @property # discriminant
    def D(self: typing.Self) -> int:
        return self.b ** 2 - 4 * self.a * self.c

    def __eq__(self: typing.Self, other: typing.Self) -> bool:
        return self.a == other.a and self.b == other.b and self.c == other.c
    
    def __hash__(self: typing.Self) -> int:
        return hash((self.a, self.b, self.c))
    
    def __iter__(self: typing.Self) -> abc.Iterator:
        return iter([self.a, self.b, self.c])
    
    def __str__(self) -> str:
        return f"{self.a}x²\t{self.b:+}xy\t{self.c:+}y².\tD={self.D}"
    
    def GL2Z_action(self: typing.Self, matrix: gl2z.GL2Z) -> typing.Self:
        """
        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms", p. 12.

        axx + bxy + cyy
        == a(alpha x + beta y)(alpha x + beta y) + b(alpha x + beta y)(gamma x + delta y) + c(gamma x + delta y)(gamma x + delta y)
        == (a*alpha**2 + b*alpha*gamma + c*gamma**2) xx
           + (2*a*alpha*beta + b**alpha*delta + b*beta*gamma + 2*c*gamma*delta) xy
           + (a*beta**2 + b*beta*delta + c*delta**2) yy
        == Axx + Bxy + Cyy
        """
        if matrix.det not in [-1, 1]:
            raise TypeError(f"{matrix=} must belong to GL₂(𝐙).")
        a, b, c = self.a, self.b, self.c
        alpha, beta, gamma, delta = matrix.alpha, matrix.beta, matrix.gamma, matrix.delta
        A = a * alpha ** 2 + b * alpha * gamma + c * gamma ** 2
        B = 2 * a * alpha * beta + b * alpha * delta + b * beta * gamma + 2 * c * gamma * delta
        C = a * beta ** 2 + b * beta * delta + c * delta ** 2
        return type(self)(A, B, C)


    @property
    def is_reduced(self: typing.Self) -> bool:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Definition 5.6.2, p. 262: an indefinite binary quadratic form (a,b,c) is reduced when |√D - 2|a|| < b < √D.
        """
        sqrtD = math.sqrt(self.D)
        return abs(sqrtD - 2 * abs(self.a)) < self.b < sqrtD

    
    def reduce_with_exponent(self: typing.Self) -> tuple[typing.Self, int]:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Definition 5.6.4 and Algorithm 5.6.5, p. 263.
        """

        def r(D: int, b: int, a: int) -> int:
            if a == 0:
                raise ValueError(f"{a=} must be nonzero")
            if abs(a) > math.sqrt(D):
                for r0 in range(-abs(a) + 1, abs(a) + 1):
                    if (r0 - b) % (2 * abs(a)) == 0:
                        return r0
            elif abs(a) < math.sqrt(D):
                for r0 in range(math.isqrt(D) - 2 * abs(a) + 1, math.isqrt(D) + 1):
                    if (r0 - b) % (2 * abs(a)) == 0:
                        return r0

        _, b, c, D = self.a, self.b, self.c, self.D
        r0 = r(D, -b, c)
        m = (b + r0) // (2 * c)

        a1 = c
        b1 = r0
        c1 = (r0 * r0 - D) // (4 * c)
        return type(self)(a1, b1, c1), m

    def reduced_with_exponent_list(self: typing.Self) -> tuple[typing.Self, list[int]]:
        bqf = self
        exponent_list = []
        while not bqf.is_reduced:
            bqf, mi = bqf.reduce_with_exponent()
            exponent_list.append(mi)
        return bqf, exponent_list

    @property
    def real_quadratic_number_associate(self: typing.Self) -> quadratic_fields.RealQuadraticNumber:
        return quadratic_fields.RealQuadraticNumber(self.D, Fraction(-self.b, 2 * abs(self.a)), Fraction(1, 2 * abs(self.a)))
    
    def integral_quadratic_polynomial_associate(self: typing.Self) -> quadratic_fields.IntegralQuadraticPolynomial:
        return self.real_quadratic_number_associate.minimal_polynomial()

    def evaluate(self: typing.Self, x: Rational | quadratic_fields.RealQuadraticNumber, y: Rational | quadratic_fields.RealQuadraticNumber) -> int:
        return self.a * x ** 2 + self.b * x * y + self.c * y **2
    
    @staticmethod
    def exponent_list_to_word(exponent_list: list[int]) -> list[str]:
        word_list = []
        for m in exponent_list:
            # gl2z.gl2z.S and gl2z.gl2z.T
            word = "S" + "T" * m
            word_list.append(word)
        return "".join(word_list)
    
    def in_GL2Q(self) -> gl2z.GL2Q:
        return gl2z.GL2Q(2 * self.a, self.b, self.b, 2 * self.c)
    
    def compose(self: typing.Self, other: typing.Self) -> typing.Self:
        pass


def class_number(D: int) -> int:
    """
    Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
    Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms", Theorem 1.6, p. 14:

    "Fix a nonsquare discriminant D.

    "(a) The Dirichlet class number h(D) is finite. In fact, any form of discriminant
    D is properly equivalent to a form (a, b, c) with |b| ≤ |a| ≤ |c| and therefore
    has 3|ac| ≤ |D|, and the number of forms of discriminant D satisfying all these
    inequalities is finite."

    "(b) An odd prime p with GCD(D, p) = 1 is primitively representable by some
    form (a, b, c) of discriminant D if and only if (D/p) = ± 1. In this case the number
    of proper equivalence classes of forms primitively representing p is either 1 or 2,
    and these classes are carried to one another by GL(2, Z). In fact, if (D/p) = ± 1,
    then b² ≡ D mod 4p for some integer b, and representatives of these classes may
    be taken to be (p, ±b, (b² - D) / 4p)."
    """


if __name__ == "__main__":
    bqf = IndefiniteBQF(2, 0, -5) # primitive indefinite BQF

    assert bqf.is_reduced and (0 < float(bqf.real_quadratic_number_associate) < 1 and -float(bqf.real_quadratic_number_associate.conjugate()) > 1) or \
        not bqf.is_reduced and not (0 < float(bqf.real_quadratic_number_associate) < 1 and -float(bqf.real_quadratic_number_associate.conjugate()) > 1)

    d = 5  # positive integer that is not a perfect square
    bqf = IndefiniteBQF(1, 0, -d)
    tau = bqf.real_quadratic_number_associate
    # The BQF is equal to 0 when (x,y) = (τ,1).
    assert bqf.evaluate(tau, 1) == quadratic_fields.RealQuadraticNumber(bqf.D, 0, 0)

    bqf = IndefiniteBQF(3, 11, 2)
    reduced_bqf, exponent_list = bqf.reduced_with_exponent_list()
    word = IndefiniteBQF.exponent_list_to_word(exponent_list)
    product_matrix = gl2z.word_to_matrix(word)
    assert bqf.GL2Z_action(product_matrix) == reduced_bqf

    m = gl2z.T
    bqf = IndefiniteBQF(1, 17, -14)
    bqf_transformed = bqf.GL2Z_action(m)
    assert bqf_transformed.D == bqf.D

    m = gl2z.S
    bqf = IndefiniteBQF(1, 0, -14)
    bqf_transformed = bqf.GL2Z_action(m)
    bqf_transformed_image = bqf_transformed.in_GL2Q()
    bqf_image = bqf.in_GL2Q()
    bqf_image_transformed = bqf_image.transpose_action_GL2Z(m)
    assert bqf_image_transformed  ==  bqf_transformed_image

    m = gl2z.P
    bqf = IndefiniteBQF(1, 0, -14)
    bqf_image = bqf.in_GL2Q()    
    assert bqf.D == -bqf_image.det



