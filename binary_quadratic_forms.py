import math
from fractions import Fraction
from numbers import Rational
from collections import abc
import itertools as it
import typing

import cflib
import gl2z
import quadratic_fields
import prime_numbers


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
    Chapter I, Section 6, "Quadratic Number Fields and Their Units", pp. 35-38.
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
        Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms", pp. 12-24.
        
        p. 12:
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

        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms", pp. 12-24.

        p. 21:
        "Let us call a primitive form (a, b, c) of discriminant D > 0 reduced when it satisfies the conditions
        0 < b < √D and √D − b < 2|a| < √D + b."
        """
        sqrtD = math.sqrt(self.D)
        return abs(sqrtD - 2 * abs(self.a)) < self.b < sqrtD

    
    def reduce_with_exponent(self: typing.Self) -> tuple[typing.Self, int]:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Definition 5.6.4 and Algorithm 5.6.5, p. 263.

        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms", pp. 12-24.

        p. 22:
        Theorem 1.8. Fix a positive nonsquare discriminant D.
        (a) Each form of discriminant D is properly equivalent to some reduced form
        of discriminant D.
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
    
    def reduced(self: typing.Self) -> typing.Self:
        bqf_reduced, _ = self.reduced_with_exponent_list()
        return bqf_reduced
    

    @staticmethod
    def are_equivalent(bqf1: typing.Self, bqf2: typing.Self) -> bool:
        """
        Proper (SL2Z) equivalence.
        """
        bqf1_reduced, _ = bqf1.reduced_with_exponent_list()
        bqf2_reduced, _ = bqf2.reduced_with_exponent_list()
        return bqf1_reduced == bqf2_reduced


    def equivalent_bqf_with_word(self: typing.Self) -> tuple[typing.Self, str]:
        """
        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms", pp. 12-24.
        
        p. 14:
        Theorem 1.6. Fix a nonsquare discriminant D.
        (a) The Dirichlet class number h(D) is finite. In fact, any form of discriminant
        D is properly equivalent to a form (a, b, c) with |b| ≤ |a| ≤ |c| and therefore
        has 3|ac| ≤ |D|, and the number of forms of discriminant D satisfying all these
        inequalities is finite.
        """
        word = ""
        bqf = self
        while not abs(bqf.b) <= abs(bqf.a) <= abs(bqf.c):
            if abs(bqf.a) > abs(bqf.c):
                bqf = bqf.GL2Z_action(gl2z.S)
                word += "S"
            n = round(Fraction(-bqf.b, 2 * bqf.a))
            if n > 0:
                bqf = bqf.GL2Z_action(gl2z.GL2Z.N(n))
                word += "T" * n
            elif n < 0:
                bqf = bqf.GL2Z_action(gl2z.GL2Z.N(n))
                word += "U" * n # U == T ** (-1)
            else:
                pass # in this case, next iteration will at most add "S" to word and will be last iteration
        return bqf, word


    @property
    def real_quadratic_number_associate(self: typing.Self) -> quadratic_fields.RealQuadraticNumber:
        return quadratic_fields.RealQuadraticNumber(self.D, Fraction(-self.b, 2 * abs(self.a)), Fraction(1, 2 * abs(self.a)))
    
    def integral_quadratic_polynomial_associate(self: typing.Self) -> quadratic_fields.IntegralQuadraticPolynomial:
        return self.real_quadratic_number_associate.minimal_polynomial()

    def evaluate(self: typing.Self, x: Rational | quadratic_fields.RealQuadraticNumber, y: Rational | quadratic_fields.RealQuadraticNumber) -> int:
        return self.a * x ** 2 + self.b * x * y + self.c * y **2
    
    def equivalent_bqf_to_evaluate(self, x0: int, y0: int) -> typing.Self:
        if math.gcd(x0, y0) != 1:
            raise ValueError(f"{x0=}, {y0=} must be relatively prime.")
        eea = cflib.EEA(x0, y0)
        a = x0
        b = -eea.bezout_y
        c = y0
        d = eea.bezout_x
        m = gl2z.GL2Z(a, b, c, d)
        bqf = self.GL2Z_action(m)
        return bqf
    
    @staticmethod
    def exponent_list_to_word(exponent_list: list[int]) -> list[str]:
        word_list = []
        for m in exponent_list:
            # gl2z.GL2Z.S and gl2z.GL2Z.T
            word = "S" + "T" * m
            word_list.append(word)
        return "".join(word_list)
    
    def in_GL2Q(self) -> gl2z.GL2Q:
        return gl2z.GL2Q(2 * self.a, self.b, self.b, 2 * self.c)

    def is_fundamental_discriminant(D: int) -> bool:
        squarefull_D, squarefree_D = prime_numbers.squarefull_and_squarefree_parts(D)
        while squarefull_D % 4 == 0:
            squarefull_D = squarefull_D // 4
        if squarefull_D == 1:
            if squarefree_D % 2 == 1 or squarefree_D % 16 in [8, 12]:
                return True
        return False
    
    @property
    def to_fundamental_discriminant(self) -> int:
        """
        Hua Loo Keng, Introduction to Number Theory, Translated from the Chinese by Peter Shiu, Springer, 1982.
        Chapter 12, "Binary Quadratic Forms", p. 322.
        Theorem 11.1. Each discriminant d is uniquely expressible as f * m ** 2 where f is a
        fundamental discriminant.
        """
        D = self.D
        _, q = prime_numbers.squarefull_and_squarefree_parts(D)
        if D % 2 == 0:
            f = q
        else:
            if q % 4 == 1:
                f = q
            elif q % 4 in [2, 3]:
                f = 4 * q
        return f


    
    @classmethod
    def principal_bqf_for_fundamental_discriminant(cls, D: int) -> typing.Self:
        if not cls.is_fundamental_discriminant(D):
            raise ValueError(f"{D=} must be a fundamental discriminant.")
        if D % 4 == 1:
            a = 1
            b = 1
            c = (1 - D) // 4
            return cls(a, b, c)
        else: # D % 4 == 0 since D is a fundamental discriminant
            d = D // 4
            a = 1
            b = 0
            c = -d
            return cls(a, b, c)

    @classmethod
    def primitively_represent_odd_prime(cls, D: int, p: int) -> list[typing.Self]:
        """
        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms", pp. 12-24.
        
        p. 14:
        Theorem 1.6. Fix a nonsquare discriminant D.
        (b) An odd prime p with GCD(D, p) = 1 is primitively representable by some
        form (a, b, c) of discriminant D if and only if (D/p) = ± 1. In this case the number
        of proper equivalence classes of forms primitively representing p is either 1 or 2,
        and these classes are carried to one another by GL(2, Z). In fact, if (D/p) = ± 1,
        then b² ≡ D mod 4p for some integer b, and representatives of these classes may
        be taken to be (p, ±b, (b² - D) / 4p).
        """
        if not cls.is_fundamental_discriminant(D):
            raise ValueError(f"{D=} must be a fundamental discriminant.")
        if p == 2 or not prime_numbers.isprime(p):
            raise ValueError(f"{p=} must be an odd prime.")
        if prime_numbers.kronecker_symbol(D, 4 * p) != 1:
            return None
        b = prime_numbers.solve_quadratic_congruence(D, 4 * p)
        bqf1 = cls(p, b, (b ** 2 - D) // (4 * p))
        bqf2 = cls(p, -b, (b ** 2 - D) // (4 * p))
        if IndefiniteBQF.are_equivalent(bqf1, bqf2):
            return [bqf1]
        return [bqf1, bqf2]

    def reduced_left_neighbor(self: typing.Self) -> typing.Self:
        if not self.is_reduced:
            raise ValueError(f"{self=!r} must be reduced.")
        D = self.D
        c = self.a
        for b0 in range(math.isqrt(D) - 2 * abs(c) + 1, math.isqrt(D) + 1):
            if (self.b + b0) % (2 * c) == 0:
                b = b0
                break
        else:
            # impossible for a reduced form: if this occurs there is error in implementation
            raise ArithmeticError(f"Could not find a valid left-neighbor for {self=!r}.")
        a = (b**2 - D) // (4 * c)
        return type(self)(a, b, c)


    def reduced_right_neighbor(self: typing.Self) -> typing.Self:
        """
        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms",
        
        p. 21:
        Two forms (a, b, c) and (a', b', c') of discriminant D > 0 will be said to be
        neighbors if c = a' and b + b' ≡ 0 mod 2c. More precisely we say in this
        case that (a', b', c') is a neighbor on the right of (a, b, c) and that (a, b, c) is
        a neighbor on the left of (a', b', c').

        p. 22:
        Theorem 1.8. Fix a positive nonsquare discriminant D.
        (b) Each reduced form of discriminant D is a neighbor on the left of one and
        only one reduced form of discriminant D and is a neighbor on the right of one
        and only one reduced form of discriminant D.
        (c) The reduced forms of discriminant D occur in uniquely determined cycles,
        each one of even length, such that each member of a cycle is an iterated neighbor
        on the right to all members of the cycle and consequently is properly equivalent
        to all other members of the cycle.
        (d) Two reduced forms of discriminant D are properly equivalent if and only
        if they lie in the same cycle in the sense of (c).
        """
        if not self.is_reduced:
            raise ValueError(f"{self=!r} must be reduced.")
        D = self.D
        a = self.c
        for b0 in range(math.isqrt(D) - 2 * abs(self.c) + 1, math.isqrt(D) + 1):
            if (self.b + b0) % (2 * self.c) == 0:
                b = b0
                break
        else:
            # impossible for a reduced form: if this occurs there is error in implementation
            raise ArithmeticError(f"Could not find a valid left‐neighbor for {self=!r}.")

        c = (b ** 2 - D) // (4 * a)
        return type(self)(a, b, c)


    @classmethod
    def count_reduced_bqf_fundamental_discriminant(cls, D: int) -> int:
        if not cls.is_fundamental_discriminant(D):
            raise ValueError(f"{D=} must be a fundamental discriminant.")
        reduced_bqf_count = 0
        bqf = cls.principal_bqf_for_fundamental_discriminant(D)
        reduced_bqf, _ = bqf.reduced_with_exponent_list()
        reduced_bqf_count += 1
        right_neighbor = reduced_bqf.reduced_right_neighbor()
        while right_neighbor != reduced_bqf:
            right_neighbor = right_neighbor.reduced_right_neighbor()
            reduced_bqf_count += 1
        return reduced_bqf_count

    def primitively_represent(self, m: int) -> tuple[int, int, int]:
        """
        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 4, "Composition of Forms, Class Group", pp. 24-31.

        p. 26:
        Lemma 1.10. If (a, b, c) is a primitive form of nonsquare discriminant and if
        m ≠ 0 is an integer, then (a, b, c) primitively represents some integer relatively
        prime to m.
        """
        if m == 0:
            raise ValueError(f"{m=} must be nonzero integer.")
        a, _, c = self
        x0 = math.prod(p for p in prime_numbers.make_prime_factor_list(math.gcd(a, m)) if c % p != 0)
        y0 = math.prod(p for p in prime_numbers.make_prime_factor_list(m) if a % p != 0)
        l = self.evaluate(x0, y0)
        return l, x0, y0
    

    @classmethod
    def compose(cls, bqf1: typing.Self, bqf2: typing.Self) -> typing.Self:
        """
        Return the reduced composition of two reduced primitive
        indefinite forms with the same positive discriminant.

        Cohen, A Course in Computational Algebraic Number Theory
        Algorithm 5.6.8
        
        Keeps all coefficients ≤ √D during the whole computation.
        """
        if bqf1.D != bqf2.D:
            raise ValueError("Discriminants must match.")
        if not (bqf1.is_reduced and bqf2.is_reduced):
            raise ValueError("Both inputs must be reduced forms.")

        a1, b1, _ = bqf1
        a2, b2, _ = bqf2
        D = bqf1.D # common discriminant

        # Step 1. g = gcd(a₁, a₂, (b₁+b₂)/2)
        r  = (b1 + b2) // 2 # b1, b2 have same parity
        g  = math.gcd(a1, a2, r)

        # Step 2. Divide the inputs by g
        a1_, a2_, r_ = a1 // g, a2 // g, r // g # pairwise coprime

        # Step 3. Solve  a1_·s + a2_·t = r_  via Bézout
        eea = cflib.EEA(a1_, a2_)                 # a1_·x + a2_·y = 1
        s0, _ = eea.bezout_x, eea.bezout_y

        s  = (s0 * r_) % a2_                      # 0 ≤ s < a2_
        t  = (r_ - s * a1_) // a2_                # ensures equality

        # Step 4. Compose
        a3 = a1_ * a2_ + s * t * g
        b3 = b1 + 2 * a1 * s
        # Keep b₃ within (−a₃, a₃] to ease the final reduction
        b3 = (b3 + a3) % (2 * a3) - a3

        c3 = (b3 * b3 - D) // (4 * a3)

        composed = cls(a3, b3, c3)

        # Step 5.  Final reduction (≤ 2 iterations)
        reduced = composed.reduced()

        # Step 6.  Canonical choice: ensure a > 0
        if reduced.a < 0:
            reduced = reduced.reduced_right_neighbor()

        return reduced
    
    def image_mod_D(self: typing.Self) -> set[int]:
        D = self.D
        domain_DxD = it.product(range(D), repeat=2)
        image_set = {self.evaluate(x, y) % D for x, y in domain_DxD}
        image_set_relatively_prime = {k for k in image_set if math.gcd(D, k) == 1}
        return image_set_relatively_prime




class ProperEquivalenceClass(abc.Container):
    def __init__(self, bqf: IndefiniteBQF) -> None:
        reduced_bqf, _ = bqf.reduced_with_exponent_list()
        self.bqf = reduced_bqf

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.bqf})"

    def __eq__(self, other: typing.Self) -> typing.Self:
        return self.bqf == other.bqf
    
    def __hash__(self) -> int:
        return hash(self.bqf)
    
    def __contains__(self, item: IndefiniteBQF) -> bool:
        return IndefiniteBQF.are_equivalent(self.bqf, item)

    @property
    def D(self) -> int:
        return self.bqf.D

    @classmethod
    def identity_class(cls, D: int) -> typing.Self:
        if not IndefiniteBQF.is_fundamental_discriminant(D):
            raise ValueError(f"{D=} must be a fundamental discriminant.")
        if D % 4 == 0:
            bqf = IndefiniteBQF(1, 0, -D // 4)
            return cls(bqf)
        elif D % 4 == 1:
            bqf = IndefiniteBQF(1, 1, -(D - 1) // 4)
            return cls(bqf)
        
    def inverse(self) -> typing.Self:
        bqf_inverse = IndefiniteBQF(self.bqf.a, -self.bqf.b, self.bqf.c)
        return type(self)(bqf_inverse)
    
    def __mul__(self, other: typing.Self) -> typing.Self:
        bqf_composed = IndefiniteBQF.compose(self.bqf, other.bqf)
        return type(self)(bqf_composed)

    def __truediv__(self, other: typing.Self) -> typing.Self:
        if self.D != other.D:
            raise ValueError(f"{self} and {other} must be associated with the same discriminant.")
        return self * other.inverse()
    

def classnumber_h(D: int) -> int:
    """
    Number of proper equivalence classes of primitive *indefinite*
    binary quadratic forms of fundamental discriminant ``D > 0``.

    The routine enumerates **all Knapp-reduced** forms

        (a, b, c)   with
            0 < b < √D                  and
            √D − b < 2|a| < √D + b      and
            4ac = b² − D ,  gcd(a,b,c)=1 ,

    then partitions them into cycles under the right-neighbour map,
    each cycle corresponding to one proper class.

    Parameters
    ----------
    D : int
        Positive, non-square, *fundamental* discriminant.

    Returns
    -------
    int
        Dirichlet class number h(D) (proper classes).

    Raises
    ------
    ValueError
        If D is not a positive fundamental discriminant.
    """
    if not IndefiniteBQF.is_fundamental_discriminant(D) or D <= 0:
        raise ValueError(f"{D} must be a positive fundamental discriminant.")

    sqrtD  = math.isqrt(D)                    # ⌊√D⌋  (integer)
    forms  : list[IndefiniteBQF] = []

    # Enumerate every Knapp-reduced primitive form exactly once
    for b in range(1, sqrtD + 1):             # 0 < b < √D  ⇒  b ≤ ⌊√D⌋
        lower = sqrtD - b                     # left inequality bound
        upper = sqrtD + b                     # right inequality bound
        # a must be positive;  2a > lower  and  2a < upper
        a_min = lower // 2 + 1
        a_max = (upper - 1) // 2              # strict inequality
        for a in range(a_min, a_max + 1):
            n = b * b - D
            den = 4 * a
            if n % den:                       # 4a ∤ b² − D  ⇒  not a form
                continue
            c = n // den
            if math.gcd(a, b, c) != 1:        # primitive only
                continue
            forms.append(IndefiniteBQF(a, b, c))

    # Group reduced forms into neighbour-cycles
    visited: set[IndefiniteBQF] = set()
    cycles = 0
    for f in forms:
        if f in visited:
            continue
        cycles += 1
        g = f
        while True:
            visited.add(g)
            g = g.reduced_right_neighbor()
            if g == f:
                break
    return cycles




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


    bqf = IndefiniteBQF(1, 0, -14)
    bqf_transformed = bqf.GL2Z_action(gl2z.S)
    assert bqf_transformed.a == bqf.c and bqf_transformed.b == - bqf.b and bqf_transformed.c == bqf.a

    bqf = IndefiniteBQF(1, 0, -14)
    n = 5
    bqf_transformed = bqf.GL2Z_action(gl2z.T ** n)
    assert bqf_transformed.a == bqf.a and bqf_transformed.b == 2 * bqf.a * n + bqf.b

    bqf = IndefiniteBQF(3, 11, 2)
    equivalent_bqf, word = bqf.equivalent_bqf_with_word()
    assert abs(equivalent_bqf.b) <= abs(equivalent_bqf.a) <= abs(equivalent_bqf.c)
    assert 3 * abs(equivalent_bqf.a * equivalent_bqf.c) <= abs(equivalent_bqf.D)
    assert abs(equivalent_bqf.a) <= abs(equivalent_bqf.D) and \
        abs(equivalent_bqf.b) <= abs(equivalent_bqf.D) and \
        abs(equivalent_bqf.c) <= abs(equivalent_bqf.D)
    
    bqf = IndefiniteBQF(3, 11, 2)
    reduced_bqf, _ = bqf.reduced_with_exponent_list()
    right_neighbor = reduced_bqf.reduced_right_neighbor()
    m = gl2z.GL2Z(0, -1, 1, (reduced_bqf.b + right_neighbor.b) // (2 * reduced_bqf.c))
    reduced_bqf_transformed = reduced_bqf.GL2Z_action(m)
    assert reduced_bqf_transformed == right_neighbor

    D = 97
    assert IndefiniteBQF.count_reduced_bqf_fundamental_discriminant(D) % 2 == 0

    bqf = IndefiniteBQF(2, 8, -5)
    bqf_reduced, _ = bqf.reduced_with_exponent_list()
    assert bqf_reduced.is_reduced

    bqf = IndefiniteBQF(1, 0, -26)
    bqf_reduced, _ = bqf.reduced_with_exponent_list()
    assert IndefiniteBQF.are_equivalent(bqf, bqf_reduced)

    m = 87
    bqf = IndefiniteBQF(2,8,-5)
    bqf_reduced, _ = bqf.reduced_with_exponent_list()
    _, x0, y0 = bqf_reduced.primitively_represent(m)
    assert math.gcd(x0, y0) == 1

    m = 87
    bqf = IndefiniteBQF(2, 8, -5)
    bqf_reduced, _ = bqf.reduced_with_exponent_list()
    _, x0, y0 = bqf_reduced.primitively_represent(m)
    bqf_equivalent = bqf_reduced.equivalent_bqf_to_evaluate(x0, y0)
    assert bqf_equivalent.evaluate(1, 0) == bqf_equivalent.a
    assert bqf_equivalent.evaluate(1, 0) == bqf_reduced.evaluate(x0, y0)

    bqf = IndefiniteBQF(2, 8, -5)
    C = ProperEquivalenceClass(bqf)
    assert bqf in C

    bqf1 = IndefiniteBQF(1, 3, -2)   # D = 17, reduced
    bqf2 = IndefiniteBQF(2, 1, -2)   # D = 17, reduced
    bqf_composed = IndefiniteBQF.compose(bqf1, bqf2)
    C1 = ProperEquivalenceClass(bqf1)
    C2 = ProperEquivalenceClass(bqf2)
    assert C1 * C2 == ProperEquivalenceClass(bqf_composed)

    # for D in range(100):
    #     if IndefiniteBQF.is_fundamental_discriminant(D):
    #         print(f"h({D})", "=", classnumber_h(D), sep="\t")


    bqf = IndefiniteBQF(1, 0, -14)
    bqf_image = bqf.in_GL2Q()    
    assert bqf.D == -bqf_image.det

    m = gl2z.S
    bqf = IndefiniteBQF(1, 0, -14)
    bqf_transformed = bqf.GL2Z_action(m)
    bqf_transformed_image = bqf_transformed.in_GL2Q()
    bqf_image = bqf.in_GL2Q()
    bqf_image_transformed = bqf_image.transpose_action_GL2Z(m)
    assert bqf_image_transformed  ==  bqf_transformed_image


    bqf1 = IndefiniteBQF(1, 3, -2)   # D = 17, reduced
    bqf2 = IndefiniteBQF(2, 1, -2)   # D = 17, reduced
    # same genus
    assert bqf1.image_mod_D() == bqf2.image_mod_D()

    bqf1 = IndefiniteBQF(1, 3, -2)   # D = 17, reduced
    bqf2 = IndefiniteBQF(2, 1, -2)   # D = 17, reduced
    assert len(bqf1.image_mod_D()) == prime_numbers.euler_totient(bqf1.D) // 2

    bqf1 = IndefiniteBQF(1, 4, -8)
    bqf2 = IndefiniteBQF(-1, 4, 8)
    bqf3 = IndefiniteBQF(3, 3, -5)
    bqf4 = IndefiniteBQF(-3, 3, 5)
    print(bqf1.D, bqf2.D, bqf3.D, bqf4.D)
    bqf1_reduced = bqf1.reduced()
    print(IndefiniteBQF.is_fundamental_discriminant(bqf1.D))

    print(IndefiniteBQF.is_fundamental_discriminant(bqf1_reduced.to_fundamental_discriminant))

    for d in range(30):
        print(d, IndefiniteBQF.is_fundamental_discriminant(d))
