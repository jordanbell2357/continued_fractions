import math
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


class IndefiniteBQF:
    """Henri Cohen, A Course in Computation Algebraic Number Theory, Chapter 5:
    Definition 5.2.3, p. 225, for binary quadratic forms.
    Definition 5.6.2, p. 262, for reduced indefinite binary quadratic forms.
    Definition 5.6.4, p. 263 for reduction operator.
    Algorithm 5.6.5, p. 263 for reduction algorithm for indefinite quadratic forms.
    """

    @classmethod
    def reduced(cls, bqf: typing.Self) -> typing.Self:
        while not bqf.is_reduced:
            bqf = bqf.reduce()
        return bqf

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
        self.gcd = math.gcd(self.a, self.b, self.c)
        self.is_primitive = self.gcd == 1
        self.is_reduced = self.D > 0 and abs(math.sqrt(self.D) - 2 * abs(self.a)) < self.b < math.sqrt(self.D)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.a}, {self.b}, {self.c})"

    def __eq__(self, other: typing.Self) -> bool:
        return self.a == other.a and self.b == other.b and self.c == other.c
    
    def __float__(self) -> float:
        return (-self.b + math.sqrt(self.D)) / (2 * abs(self.a))
    
    def __rmul__(self, other: SL2Z) -> typing.Self:
        a = self.a * other.alpha ** 2 + self.b * other.alpha * other.gamma + self.c * other.gamma ** 2
        b = 2 * self.a * other.alpha * other.beta + self.b * other.alpha * other.delta + self.b * other.beta * other.gamma + 2 * self.c * other.gamma * other.delta
        c = self.a * other.beta ** 2 + self.b * other.beta * other.delta + self.c * other.delta ** 2
        return type(self)(a, b, c)
    
    def __str__(self) -> str:
        return f"{self.a}x² {self.b:+}xy {self.c:+}y². Discriminant: {self.D}"
    
    def reduce(self: typing.Self) -> typing.Self:
        def r(D: int, b: int, a: int) -> int:
            if abs(a) > math.sqrt(D):
                for k in range(-abs(a), abs(a) + 1):
                    if (k - b) % (2 * a) == 0:
                        return k
            elif abs(a) < math.sqrt(D):
                for k in range(math.isqrt(D) - 2 * abs(a), math.isqrt(D) + 1):
                    if (k - b) % (2 * a) == 0:
                        return k
        a = self.c
        b = r(self.D, -self.b, self.c)
        c = (b ** 2 - self.D) // (4 * self.c)
        return type(self)(a, b, c)

    def evaluate(self, x: int, y: int) -> int:
        return self.a * x ** 2 + self.b * x * y + self.c * y **2



bqf1 = IndefiniteBQF(2, 0, -6)
print(bqf1, bqf1.is_reduced, bqf1.is_primitive, float(bqf1))
bqf2 = bqf1.reduce()
print(bqf2, bqf2.is_reduced, bqf2.is_primitive, float(bqf2))
bqf3 = bqf2.reduce()
print(bqf3, bqf3.is_reduced, bqf3.is_primitive, float(bqf3))
m = SL2Z(-1, -1, 0, -1)
bqfm = m * bqf1
print(bqfm, bqfm.is_reduced, bqfm.is_primitive, float(bqfm))

bqf_reduced = IndefiniteBQF.reduced(bqf1)

print(bqf_reduced, bqf_reduced.is_reduced, bqf_reduced.is_primitive, float(bqf_reduced))

