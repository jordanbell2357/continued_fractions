import math
from fractions import Fraction
from numbers import Rational
from collections import Counter
from collections import abc
import functools as ft
import operator
import itertools as it
import typing

import cflib
import gl2z
import pell
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

    Daniel E. Flath, Introduction to Number Theory, Wiley, 1989.

    Franz Halter-Koch, Quadratic Irrationals: An Introduction to Classical Number Theory, CRC Press, 2013.

    Primes of the Form 𝑥²+𝑛𝑦²: Fermat, Class Field Theory, and Complex Multiplication. Third Edition with Solutions
    David A. Cox, with contributions by Roger Lipsett
    AMS Chelsea Publishing: An Imprint of the American Mathematical Society, 2022
    § 1.2.A, "Lagrange, Legendre and Quadratic Forms", pp. 20-22
    § 1.2.C, "Elementary Genus Theory", pp. 27-31
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

    @property
    def D(self: typing.Self) -> int:
        """
        Discriminant.
        D = b^2 - 4ac
        D ≡ 0, 1 (mod 4) because D ≡ b^2 (mod 4) and b^2 ≡ 0, 1 (mod 4)
        """
        return self.b ** 2 - 4 * self.a * self.c

    def __eq__(self: typing.Self, other: typing.Self) -> bool:
        return self.a == other.a and self.b == other.b and self.c == other.c
    
    def __hash__(self: typing.Self) -> int:
        return hash((self.a, self.b, self.c))

    def __iter__(self: typing.Self) -> abc.Iterator:
        return iter([self.a, self.b, self.c])
    
    def __abs__(self: typing.Self) -> int:
        return max(abs(x) for x in self)
    
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
    

    def stabilizer_GL2Z(self) -> gl2z.GL2Z:
        """
        Duncan A. Buel, Binary Quadratic Forms: Classical Theory and Modern Computations 
        Springer, 1989

        Section 3.2 Automorphs, Pell's Equation, p. 31

        We recall that an
        automorph of a binary quadratic form is a nontrivial transformation
        (1.1) of determinant +1 under which the form is equivalent to itself. 

        Theorem 3.9. If Δ is any discriminant of binary quadratic forms,
        then there exists a solution (x, y) to the Pell equation
        x² -Δy² = 4. (3.2)
        There is a one-to-one correspondence between automorphs of (definite
        or indefinite) forms (a, b, c) of discriminant Δ and solutions of the
        Pell equation (3.2).

        https://mathoverflow.net/questions/344562/infinite-cyclic-subgroups-of-textsl-2-mathbbz
        """

        a, b, c = self
        D = self.D
        x, y = pell.solve_pell_equation(D)
        t = 2 * x
        u = 2 * y
        g = gl2z.GL2Z((t + b * u) // 2, -c * u, a * u, (t - b * u) // 2)
        return g


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
        a, b, c = self
        D = self.D
        sqrtD = math.sqrt(D)
        # return abs(sqrtD - 2 * abs(self.a)) < self.b < sqrtD
        return 0 < b < sqrtD and sqrtD - b < 2 * abs(a) < sqrtD + b

    
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
        r1 = bqf1.reduced()
        r2 = bqf2.reduced()
        if r1 == r2:
            return True                    # same representative

        rn = r1.reduced_right_neighbor()  # step once round the cycle
        while rn != r1:                   # finite because cycles are finite
            if rn == r2:
                return True
            rn = rn.reduced_right_neighbor()
        return False                       # walked full cycle – different class


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
        return self.a * x ** 2 + self.b * x * y + self.c * y ** 2


    def equivalent_bqf_to_evaluate(self: typing.Self, x0: int, y0: int) -> typing.Self:
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
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Definition 5.1.2, p. 224.
        """
        if D == 1:
            return False
        if prime_numbers.is_squarefree(D) and D % 4 == 1:
            return True
        elif D % 4 == 0:
            m = D // 4
            if prime_numbers.is_squarefree(m) and m % 4 in [2, 3]:
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
        if D % 2 == 1:
            return q
        else: # D % 2 == 0
            if q % 4 == 1:
                return q
            elif q % 4 in [2, 3]:
                return 4 * q


    @classmethod
    def principal_bqf_for_discriminant(cls, D: int) -> typing.Self:
        if D % 4 not in [0, 1]:
            raise ValueError(f"{D=} must be a discriminant.")
        if D % 4 == 1:
            a = 1
            b = 1
            c = (1 - D) // 4
            return cls(a, b, c)
        else: # D % 4 == 0 since D is a discriminant
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
        if D % 4 not in [0, 1]:
            raise ValueError(f"{D=} must be a discriminant.")
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
        """
        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 3, "Equivalence and Reduction of Quadratic Forms", pp. 12-24.
        
        p. 21:
        Two forms (a, b, c) and (a', b', c') of discriminant D > 0 will be said to be
        neighbors if c = a' and b + b' ≡ 0 mod 2c. More precisely we say in this
        case that (a', b', c') is a neighbor on the right of (a, b, c) and that (a, b, c) is
        a neighbor on the left of (a', b', c').

        p. 22:
        Theorem 1.8. Fix a positive nonsquare discriminant D.
        (a) Each form of discriminant D is properly equivalent to some reduced form
        of discriminant D.
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
    def list_reduced_bqf_fundamental_discriminant(cls, D: int) -> int:
        if not cls.is_fundamental_discriminant(D):
            raise ValueError(f"{D=} must be a fundamental discriminant.")
        reduced_bqf_list = []
        bqf = cls.principal_bqf_for_discriminant(D)
        reduced_bqf= bqf.reduced()
        reduced_bqf_list.append(reduced_bqf)
        right_neighbor = reduced_bqf.reduced_right_neighbor()
        while right_neighbor != reduced_bqf:
            reduced_bqf_list.append(right_neighbor)
            right_neighbor = right_neighbor.reduced_right_neighbor()
        return reduced_bqf_list


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


    def inverse(self: typing.Self) -> typing.Self:
        """
        Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
        Chapter I, Section 4, "Composition of Forms, Class Group", pp. 24-31.

        Theorem 1.12, pp. 27-28.
        (c) Under the product operation in (b), the set of proper equivalence classes
        of primitive forms of discriminant D is a finite abelian group. The identity is the
        class of (1, 0, −D/4) if D ≡ 0 mod 4 and is the class of (1, 1, −(D − 1)/4) if
        D ≡ 1 mod 4. The group inverse of the class of (a, b, c) is the class of (a, −b, c).   
        """
        return type(self)(self.a, -self.b, self.c)


    def image_mod_D(self: typing.Self) -> set[int]:
        D = self.D
        domain_DxD = it.product(range(D), repeat=2)
        image_set = {self.evaluate(x, y) % D for x, y in domain_DxD}
        image_set_relatively_prime = {k for k in image_set if math.gcd(D, k) == 1}
        return image_set_relatively_prime


    def min_abs_image(self: typing.Self) -> set[int]:
        """
        m(f)=inf|f(x,y)|,  x,y∈𝐙²,  (x,y)≠(0,0)
        """
        D = self.D
        domain_DxD = it.product(range(-math.isqrt(D) - 1, math.isqrt(D) + 1), repeat=2)
        abs_image_set = {abs(self.evaluate(x, y)) for x, y in domain_DxD}
        abs_image_set.remove(0)
        return min(abs_image_set)


    def min_abs_image_by_reduction(self: typing.Self) -> int:
        """
        m(f) = min { |f(x,y)| : |f(x,y)| != 0, x,y in Z }

        Johannes Buchmann and Ulrich Vollmer
        Binary Quadratic Forms: An Algorithmic Approach
        Springer, 2007

        p. 139
        Corollary 6.16.2. The minimum of an integral indefinite irreducible form f
        is the absolute smallest integer that appears as an a in a form (a, b, c) in the
        proper cycle of f.
        """
        # Step 1: one reduction
        g = self.reduced()

        # Step 2: walk the cycle, take the smallest |a|
        min_abs = abs(g.a)
        neighbour = g.reduced_right_neighbor()
        while neighbour != g:
            min_abs = min(min_abs, abs(neighbour.a))
            if min_abs == 1:      # early exit – cannot do better
                break
            neighbour = neighbour.reduced_right_neighbor()

        return min_abs


    @property
    def conductor(self) -> int:
        """
        Johannes Buchmann and Ulrich Vollmer
        Binary Quadratic Forms: An Algorithmic Approach
        Springer, 2007

        Definition 3.3.2, p. 37:
        The conductor of a discriminant ∆ is the largest positive integer f such
        that ∆/f^2 is a discriminant. We denote it by f(∆).
        """
        D = self.D
        factor_D = prime_numbers.make_prime_factor_counter(D)
        factor_f  = Counter()

        for p, v in factor_D.items():
            if p == 2:
                # ----------- subtle 2-adic cases -------------------------
                t = D // (2 ** v)               # t ≡ 1 or 3 (mod 4)
                if t % 4 == 1:                  # Case 1
                    if v % 2 == 0:              #   v even
                        e = v // 2
                    else:                       #   v odd
                        e = 0 if v < 3 else (v - 3) // 2
                else:                           # t % 4 == 3   (Case 2)
                    e = (v - v % 2 - 2) // 2
                if e > 0:
                    factor_f[2] = e
            else:
                # ----------- odd primes: floor(v/2) ----------------------
                if v >= 2:
                    factor_f[p] = v // 2

        f = math.prod(p ** e for p, e in factor_f.items())
        return f if f > 0 else 1


    def lift(self, f: int) -> typing.Self:
        """
        Return a primitive form whose discriminant is f²·D.

        Strategy
        --------
        Two lifts are possible:

        *  y ↦ f·y  →  (a,  f·b,  f²·c)
        *  x ↦ f·x  →  (f²·a,  f·b,   c)

        Exactly one of them is guaranteed to stay primitive
        (because a primitive triple (a,b,c) cannot have *both*
        a and c divisible by f² when f > 1).

        Parameters
        ----------
        f : positive integer  (the desired conductor multiplier)

        Raises
        ------
        ValueError
            if *f* ≤ 0 or if neither orientation yields a primitive form.
        """
        if not isinstance(f, int) or f <= 0:
            raise ValueError(f"{f=} must be a positive integer.")

        if f == 1:
            return self

        a, b, c = self.a, self.b, self.c

        # orientation 1 : y  ->  f·y
        a1, b1, c1 = a, b * f, c * f * f
        if math.gcd(a1, b1, c1) == 1:
            return type(self)(a1, b1, c1)

        # orientation 2 : x  ->  f·x
        a2, b2, c2 = a * f * f, b * f, c
        if math.gcd(a2, b2, c2) == 1:
            return type(self)(a2, b2, c2)

        # occurs when lifted form is not primitive
        raise ValueError(f"Lift by conductor {f=} does not preserve primitivity.")


    def descend(self) -> typing.Self:
        """
        Use SL₂(𝐙) shear to drop the conductor f of the
        discriminant Δ = f²·Δ₀ to obtain a primitive form of the
        fundamental discriminant Δ₀.

        The procedure tries, for n = 0,...,f−1:

        1.  Right shear Tⁿ = (1 n 0 1):
            succeeds if c' ≡ 0 (mod f²).
            Return (a , b'/f , c'/f²).

        2.  Left shear Uⁿ = (1 0 n 1): (orientation 2)
            succeeds if a' ≡ 0 (mod f²).
            Return (a'/f² , b'/f , c).
        """
        f = self.conductor
        if f == 1:
            return self

        # Right shear
        for n in range(f):
            g = gl2z.GL2Z.N(n)            # (1 n; 0 1)
            bqf = self.GL2Z_action(g)       # coefficients (a , b', c')
            if bqf.b % f == 0 and bqf.c % f ** 2 == 0:
                return type(self)(bqf.a, bqf.b // f, bqf.c // f ** 2)

        # Left shear
        for n in range(f):
            g = gl2z.GL2Z(1, 0, n, 1)     # (1 0; n 1)
            bqf = self.GL2Z_action(g)       # coefficients (a', b', c )
            if bqf.b % f == 0 and bqf.a % f ** 2 == 0:
                return type(self)(bqf.a // f ** 2, bqf.b // f, bqf.c)

        # If both loops fail (should not happen if implemented correctly, by theory)
        raise ArithmeticError("Descent failed: no suitable SL₂(𝐙) shear found.")

    @classmethod
    def compose_bqf(cls, bqf1: typing.Self, bqf2: typing.Self) -> typing.Self:
        """
        Duncan A. Buel, Binary Quadratic Forms: Classical Theory and Modern Computations 
        Springer, 1989
        Theorem 4.10, p. 62.
        """
        if bqf1.D != bqf2.D:
            raise ValueError(f"{bqf1} and {bqf2} must have the same discriminant.")
        bqf1 = bqf1.reduced()
        bqf2 = bqf2.reduced()
        D = bqf1.D
        a1, b1, _, a2, b2, _ = *bqf1, *bqf2
        beta = (b1 + b2) // 2
        eea = cflib.EEA(a1, a2)
        m = eea.gcd
        t1, t2 = eea.bezout_x, eea.bezout_y
        # t1 a1 + t2 a2 = m
        eea = cflib.EEA(m, beta)
        n = eea.gcd # n = math.gcd(a1, a2, beta)
        u1, u2 = eea.bezout_x, eea.bezout_y
        # u1 m + u2 beta = n
        # u1(t1 a1 + t2 a2) + u2 beta = n
        t, u, v = u1 * t1, u1 * t2, u2
        A = a1 * a2 // n ** 2
        B = (a1 * b2 * t + a2 * b1 * u + v * (b1 * b2 + D) // 2) // n
        # D = B^2 - 4AC => C = (B^2 - D) / (4A)
        C = (B ** 2 - D) // (4 * A)
        return cls(A, B, C).reduced()


    @classmethod
    def compose(cls, bqf1: typing.Self, bqf2: typing.Self) -> typing.Self:
        """
        Henri Cohen, A Course in Computation Algebraic Number Theory, Graduate Texts in Mathematics, Volume 138, Springer, 1996.
        Definition 5.4.6, p. 246.
        """
        if bqf1.D != bqf2.D:
            raise ValueError(f"{bqf1} and {bqf2} must have the same discriminant.")
        bqf1 = bqf1.reduced()
        bqf2 = bqf2.reduced()
        D = bqf1.D
        a1, b1, c1, a2, b2, c2 = *bqf1, *bqf2
        s = (b1 + b2) // 2
        n = (b1 - b2) // 2
        eea = cflib.EEA(a1, a2)
        u1, v1 = eea.bezout_x, eea.bezout_y
        d1 = eea.gcd
        eea = cflib.EEA(d1, s)
        u2, v2 = eea.bezout_x, eea.bezout_y
        d = eea.gcd
        # u1 a1 + v1 a2 = d1
        # u2 d1 + v2 s = d
        # u2(u1 * a1 + v1 * a2) + v2 * s = d
        # u1u2 a1 + v1u2 a2 + v2 s = d
        u, v, w = u1 * u2, v1 * u2, v2
        d0 = math.gcd(d, c1, c2, n)
        a3 = d0 * a1 * a2 // d ** 2
        b3 = b2 + 2 * a2 // d * (v * (s - b2) - w * c2)
        c3 = (b3 ** 2 - D) // (4 * a3)
        return cls(a3, b3, c3).reduced()
    
    
    def __mul__(self, other: typing.Self) -> typing.Self:
        return type(self).compose(self, other)

    def __rmul__(self, other: typing.Self) -> typing.Self:
        return type(self).compose(other, self)
    
    def __truediv__(self, other: typing.Self) -> typing.Self:
        bqf1 = self
        bqf2 = other.inverse()
        return bqf1 * bqf2

    def __rtruediv__(self, other: typing.Self) -> typing.Self:
        bqf1 = other
        bqf2 = self.inverse()
        return bqf1 * bqf2
    
    def __pow__(self, n: int) -> typing.Self:
        if not isinstance(n, int):
            raise TypeError(f"Exponent {n=} must be integer.")
        if n == 0:
            return type(self).principal_bqf_for_discriminant(self.D)
        elif n > 0:
            bqf = self
            return ft.reduce(operator.mul, n * [bqf])
        else: # n < 0
            bqf_inverse = self.inverse()
            return ft.reduce(operator.mul, (-n) * [bqf_inverse])
        

def class_group(D: int) -> list[IndefiniteBQF]:
    """
    For a positive *fundamental* discriminant *D* return **one reduced
    representative per proper equivalence class** of primitive indefinite
    binary quadratic forms of discriminant *D*.

    Theory
    ------
    A primitive form (a,b,c) with D = b²-4ac > 0 is when  
        0 < b < √D  and  |√D − 2a| <  b.

    Every class contains exactly one reduced form, so listing all reduced
    forms – modulo equality – gives the whole class group.

    Algorithm
    ---------
    1.  Loop through  1 ≤ b < √D  with the parity constraint  b² ≡ D (mod 4).
    2.  For each such b, loop through 1 ≤ a < √D, keep those that satisfy  
           √D − b < 2a < √D + b   and   c = (b² − D)/(4a) ∈ ℤ.
    3.  Require gcd(a,b,c)=1 (primitivity) and a>0 (always true for reduced
        forms).  Append IndefiniteBQF(a,b,c).reduced() to set.
    """

    if not IndefiniteBQF.is_fundamental_discriminant(D) or D <= 0:
        raise ValueError(f"{D=} must be a positive fundamental discriminant.")

    sqrtD = math.sqrt(D)
    isqrtD = math.isqrt(D)
    bqf_list: list[IndefiniteBQF] = []

    for b in range(1, isqrtD + 1):
        if (b ** 2 - D) % 4: # b² ≡ D (mod 4)
            continue
        for a in range(1, isqrtD + 1):
            if not (sqrtD - b < 2 * a < sqrtD + b): # 2a < √D + b < 2√D
                continue
            numerator = b * b - D # negative
            denominator = 4 * a # positive
            if numerator % denominator != 0:
                continue
            c = numerator // denominator # c < 0
            if math.gcd(a, b, c) != 1: # continue if not primitive
                continue
            bqf = IndefiniteBQF(a, b, c).reduced()   # canonical triple
            if not any(IndefiniteBQF.are_equivalent(bqf, bqfi) for bqfi in bqf_list):
                bqf_list.append(bqf)

    bqf_list = sorted(bqf_list, key=lambda bqf: (bqf.a, bqf.b, bqf.c)) # lexicographic order
    return bqf_list


def genus_group_order(D: int) -> int:
    """
    Anthony W. Knapp, Advanced Algebra, Digital Second Edition, 2016.
    Chapter I, Section 5, "Genera", pp. 31-34.

    p. 33:
    Theorem 1.14. For a fundamental discriminant D, the principal genus P of
    primitive integer forms is a subgroup of the form class group H, and the cosets
    of P are the various genera. Thus the set G of genera is exactly the set of cosets
    H/P and inherits a group structure from class multiplication. The subgroup P
    coincides with the subgroup of squares in H, and consequently every nontrivial
    element of G has order 2.
    REMARKS. The group G is called the genus group of discriminant D. The
    hypothesis that D is fundamental is needed only for the conclusion that every
    member of P is a square in H. Since every nontrivial element of G has order
    2 when D is fundamental, application of the Fundamental Theorem of Finitely
    Generated Abelian Groups or use of vector-space theory over a 2-element field
    shows that G is the direct sum of cyclic groups of order 2; in particular, the order
    of G is a power of 2. Problems 25–29 at the end of the chapter show that the
    order of G is 2g, where g + 1 is the number of distinct prime factors of D.
    """

    if not IndefiniteBQF.is_fundamental_discriminant(D) or D <= 0:
        raise ValueError(f"{D} must be a positive fundamental discriminant.")
    bqf = IndefiniteBQF.principal_bqf_for_discriminant(D)
    return prime_numbers.euler_totient(D) // len(bqf.image_mod_D())


def genus_group_order_by_divisors(D: int) -> int:
    if not IndefiniteBQF.is_fundamental_discriminant(D) or D <= 0:
        raise ValueError(f"{D} must be a positive fundamental discriminant.")
    s = prime_numbers.prime_little_omega(D)
    if D % 4 == 1 or D % 8 in [0, 4]:
        s += 1
    return 2 ** (s - 1)


def genus_group(D: int) -> list[list[IndefiniteBQF]]:
    """
    For a positive *fundamental* discriminant *D* return the list of genera,
    each genus being given as a list of reduced representatives of the
    SL₂(ℤ)-classes it contains.

    Theory
    ------
    *  The form class group **H** is the set of proper-equivalence classes of
       primitive forms of discriminant *D*.
    *  Every form f has an “image set”
          f.image_mod_D() ⊆ (ℤ/Dℤ)*
       and two classes lie in the **same genus** ⇔ their image sets coincide
       (Knapp I§5).  Thus a genus is determined by that residue set.
    *  The number of distinct genera is
          |G| = genus_group_order(D)

    Algorithm
    ---------
    1.  Get one reduced representative for every class with class_group.
    2.  Bucket those representatives by the key  
           key(f) := frozenset(f.image_mod_D()).
        Each bucket is a genus.
    3.  Sort the lists.

    Returns
    -------
    A list whose length is `genus_group_order(D)`;  each entry is a
    list of `IndefiniteBQF` instances belonging to one genus,
    with the first list always the *principal* genus (it contains the
    principal form (1,1,(1−D)/4) or (1,0,−D/4)).

    """
    if not IndefiniteBQF.is_fundamental_discriminant(D) or D <= 0:
        raise ValueError(f"{D=} must be a positive fundamental discriminant.")

    # 1.  All classes
    classes = class_group(D)

    # 2.  Bucket by genus key
    genus_buckets: dict[frozenset[int], list[IndefiniteBQF]] = {}
    for f in classes:
        key = frozenset(f.image_mod_D())
        genus_buckets.setdefault(key, []).append(f)

    # 3.  Normalise / pretty-print
    genus_list = [sorted(bucket, key=lambda g: (g.a, g.b, g.c))
                  for bucket in genus_buckets.values()]
    genus_list.sort(key=lambda bucket: (bucket[0].a, bucket[0].b, bucket[0].c))
    return genus_list


if __name__ == "__main__":
    bqf = IndefiniteBQF(2, 0, -5) # primitive indefinite BQF
    assert (
        (bqf.is_reduced
            and (0 < float(bqf.real_quadratic_number_associate) < 1
            and -float(bqf.real_quadratic_number_associate.conjugate()) > 1))
            or
        (not bqf.is_reduced
         and not (0 < float(bqf.real_quadratic_number_associate) < 1
         and -float(bqf.real_quadratic_number_associate.conjugate())) > 1)
         )
    
    bqf = IndefiniteBQF(3, 11, 2)
    equivalent_bqf, word = bqf.equivalent_bqf_with_word()
    assert abs(equivalent_bqf.b) <= abs(equivalent_bqf.a) <= abs(equivalent_bqf.c)
    assert 3 * abs(equivalent_bqf.a * equivalent_bqf.c) <= abs(equivalent_bqf.D)
    assert abs(equivalent_bqf.a) <= abs(equivalent_bqf.D) and \
        abs(equivalent_bqf.b) <= abs(equivalent_bqf.D) and \
        abs(equivalent_bqf.c) <= abs(equivalent_bqf.D)

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
    reduced_bqf = bqf.reduced()
    right_neighbor = reduced_bqf.reduced_right_neighbor()
    m = gl2z.GL2Z(0, -1, 1, (reduced_bqf.b + right_neighbor.b) // (2 * reduced_bqf.c))
    reduced_bqf_transformed = reduced_bqf.GL2Z_action(m)
    assert reduced_bqf_transformed == right_neighbor

    D = 97
    assert len(IndefiniteBQF.list_reduced_bqf_fundamental_discriminant(D)) % 2 == 0

    bqf = IndefiniteBQF(2, 8, -5)
    bqf_reduced = bqf.reduced()
    assert bqf_reduced.is_reduced

    bqf = IndefiniteBQF(1, 0, -26)
    bqf_reduced = bqf.reduced()
    assert IndefiniteBQF.are_equivalent(bqf, bqf_reduced)

    m = 87
    bqf = IndefiniteBQF(2,8,-5)
    bqf_reduced = bqf.reduced()
    _, x0, y0 = bqf_reduced.primitively_represent(m)
    assert math.gcd(x0, y0) == 1

    m = 87
    bqf = IndefiniteBQF(2, 8, -5)
    bqf_reduced = bqf.reduced()
    _, x0, y0 = bqf_reduced.primitively_represent(m)
    bqf_equivalent = bqf_reduced.equivalent_bqf_to_evaluate(x0, y0)
    assert bqf_equivalent.evaluate(1, 0) == bqf_equivalent.a
    assert bqf_equivalent.evaluate(1, 0) == bqf_reduced.evaluate(x0, y0)

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

    bqf = IndefiniteBQF(1, 0, -14)
    g = bqf.stabilizer_GL2Z()
    assert bqf.GL2Z_action(g) == bqf

    bqf = IndefiniteBQF(1, 0, -14)
    assert bqf.min_abs_image() == bqf.min_abs_image_by_reduction()

    # D = 205
    bqf1 = IndefiniteBQF(1, 13, -9)
    bqf2 = IndefiniteBQF(-1, -13, 9)
    bqf3 = IndefiniteBQF(3, 13, -3)
    bqf4 = IndefiniteBQF(-3, 13, 3)
    assert bqf1.image_mod_D() == bqf2.image_mod_D()
    assert bqf3.image_mod_D() == bqf4.image_mod_D()
    assert bqf1.image_mod_D() != bqf3.image_mod_D()

    # D = 17, reduced, same genus
    bqf1 = IndefiniteBQF(1, 3, -2)
    bqf2 = IndefiniteBQF(2, 1, -2)
    assert bqf1.image_mod_D() == bqf2.image_mod_D()

    for D in range(200):
        if IndefiniteBQF.is_fundamental_discriminant(D):
            assert genus_group_order(D) == genus_group_order_by_divisors(D)

    bqf1 = IndefiniteBQF.principal_bqf_for_discriminant(45)
    fundamental_D = bqf1.to_fundamental_discriminant
    bqf2 = IndefiniteBQF.principal_bqf_for_discriminant(fundamental_D)
    assert bqf2.conductor == 1
    assert bqf2.D == bqf1.D // bqf1.conductor ** 2

    # 1.  Round-trip sanity:  descend ∘ lift ≡ identity  (checks both directions)
    D = 5
    f = 3 # conductor
    bqf = IndefiniteBQF.principal_bqf_for_discriminant(D)       # fundamental Δ0=5
    bqf_lift = bqf.lift(f)                                            # Δ = 3²·5 = 45
    assert IndefiniteBQF.are_equivalent(bqf_lift.descend(), bqf)      # proper-equivalence test

    # 2.  Discriminant-vs-conductor consistency after a lift
    assert bqf_lift.D == f ** 2 * bqf.D                                    # Δ scales by f²
    assert bqf_lift.conductor == f                                   # f(Δ)=3 after the lift
    assert math.gcd(bqf_lift.a, bqf_lift.b, bqf_lift.c) == 1                     # primitivity preserved

    # 3.  To fundamental discriminant: a non-fundamental form descends to conductor 1
    bqf  = IndefiniteBQF.principal_bqf_for_discriminant(20)      # Δ=20, conductor f=2
    bqf_descend  = bqf.descend()                                           # should drop to Δ0=5
    assert bqf.conductor == 2 and bqf_descend.conductor == 1               # before / after check
    assert bqf_descend.D == bqf.D // (bqf.conductor ** 2)                    # Δ0 = Δ/f²

    # 4.  Even discriminant with f = 2  (checks descent after a 4-factor)
    D = 32
    bqf = IndefiniteBQF.principal_bqf_for_discriminant(D)      # (1,0,−8)  Δ = 32
    assert IndefiniteBQF.are_equivalent(bqf.descend(),
                                        IndefiniteBQF.principal_bqf_for_discriminant(D // bqf.conductor ** 2))   # Δ₀ = 8

    # 5.  Round-trip with a composite conductor 6 = 2·3
    D = 13
    f = 6
    bqf = IndefiniteBQF.principal_bqf_for_discriminant(D)     # fundamental  Δ₀ = 13
    bqf_lift = bqf.lift(f)                                           # Δ = 6²·13 = 468 (orientation chosen automatically)
    assert bqf_lift.conductor == f and bqf_lift.D == f ** 2 * bqf.D
    assert IndefiniteBQF.are_equivalent(bqf_lift.descend(), bqf)       # back to fundamental form

    # 6.  A lift that *must* fail (a and c already share the prime factor 2)
    bqf = IndefiniteBQF(2, 1, -2)    # primitive, reduced, Δ = 17 (fundamental)
    f_invalid = 2
    try:
        bqf.lift(f_invalid)
        assert False, "Must raise ValueError for invalid conductor (non-primitive lift)."
    except ValueError: # correct behaviour
        pass

    assert IndefiniteBQF.primitively_represent_odd_prime(5, 3) is None

    bqf_list = IndefiniteBQF.primitively_represent_odd_prime(205, 3)
    assert len(bqf_list) == 2                 # both inequivalent orbits returned
    assert not IndefiniteBQF.are_equivalent(*bqf_list)
    for g in bqf_list:                        # sanity: coefficient a=p and disc=D
        assert g.a == 3 and g.D == 205

    D = 17
    print(class_group(D))