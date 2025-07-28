import math
from numbers import Rational
from numbers import Number
from fractions import Fraction
from dataclasses import dataclass
from collections import abc
import typing

import gl2z
import quadratic_fields
import farey



@dataclass
class Point:
    x: quadratic_fields.RealQuadraticCompositum
    y: quadratic_fields.RealQuadraticCompositum

    def __str__(self) -> str:
        return f"({self.x}, {self.y})"


@dataclass
class Circle:
    a: Rational
    b: Rational
    r: Rational


def power_of_point_wrt_circle(point_P: Point, circle_O: Circle) -> Number:
    dist_P_O_squared = (point_P.x - circle_O.a) ** 2 + (point_P.y - circle_O.b) ** 2
    return dist_P_O_squared - circle_O.r ** 2


class RadicalLine(abc.Container):
    def __init__(self: typing.Self, c1: Circle, c2: Circle) -> None:
        self.c1 = c1
        self.c2 = c2

    def __contains__(self: typing.Self, point_P: Point) -> bool:
        return power_of_point_wrt_circle(point_P, self.c1) == power_of_point_wrt_circle(point_P, self.c2)


def intersect_circles(c1: Circle, c2: Circle,) -> tuple[Point, Point]:
    """
    Intersection of two circles in the Euclidean plane.
    Euclid, Elements:
    III.9
    III.35
    III.36
    I.47
    """
    a1, b1, r1, a2, b2, r2 = c1.a, c1.b, c1.r, c2.a, c2.b, c2.r

    # Step 1 – centres and distance
    a_delta = a2 - a1
    b_delta = b2 - b1
    d_c1_c2_squared = a_delta * a_delta + b_delta * b_delta
    if d_c1_c2_squared == 0: # concentric
        return None
    d_c1_c2 = quadratic_fields.RealQuadraticCompositum.sqrt_rational(d_c1_c2_squared)

    # Step 2 – distance d₁ of foot D from C₁
    d_foot_c1 = (d_c1_c2_squared + r1 * r1 - r2 * r2) / (2 * d_c1_c2)

    # Step 3 – half-chord length h
    h_squared = r1 * r1 - d_foot_c1 * d_foot_c1
    if h_squared < 0: # no intersection in reals
        return None
    h = quadratic_fields.RealQuadraticCompositum.sqrt_rational(h_squared)

    # Coordinates
    ux, uy = a_delta / d_c1_c2, b_delta / d_c1_c2 # û
    x0 = a1 + d_foot_c1 * ux # D = C₁ + d₁·û
    y0 = b1 + d_foot_c1 * uy
    nx, ny = -b_delta / d_c1_c2, a_delta / d_c1_c2 # n̂ perpendicular to û
    px, py = x0 + h * nx, y0 + h * ny  # P = D + h·n̂
    qx, qy = x0 - h * nx, y0 - h * ny  # Q = D − h·n̂
    p1 = Point(px, py)
    p2 = Point(qx, qy)
    return p1, p2


class FordCircle(abc.Container):
    """
    Ford circles.
    """

    def __init__(self, rational_number: Rational) -> None:
        if isinstance(rational_number, Rational):
            rational_number = Fraction(rational_number)
        else:
            raise TypeError(f"{rational_number=} must be Rational.")
        self.reduced_fraction = rational_number

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.reduced_fraction})"

    @property
    def a(self) -> Fraction:
        return Fraction(self.reduced_fraction.numerator, self.reduced_fraction.denominator)
    
    @property
    def b(self) -> Fraction:
        return Fraction(1, 2 * self.reduced_fraction.denominator ** 2)
    
    @property
    def r(self) -> Fraction:
        return Fraction(1, 2 * self.reduced_fraction.denominator ** 2)
    
    def __str__(self) -> str:
        return f"FordCircle: Center ({self.a}, {self.b}) and radius {self.r}."
    
    def __contains__(self, item: tuple[Rational, Rational]):
        """
        Check for rational  points on circle.
        """
        x, y = item[0], item[1]
        a, b, r = self.a, self.b, self.r
        return (x - a) ** 2 + (y - b) ** 2 == r ** 2
    
    @classmethod
    def are_tangent(cls, circle1: typing.Self, circle2: typing.Self) -> bool:
        a, b = circle1.reduced_fraction.numerator, circle1.reduced_fraction.denominator
        c, d = circle2.reduced_fraction.numerator, circle2.reduced_fraction.denominator
        return (a * d - b * c) ** 2 == 1
    
    @classmethod
    def try_find_point_of_tangency(cls, circle1: typing.Self, circle2: typing.Self) -> tuple[Fraction, Fraction]:
        if not cls.are_tangent(circle1, circle2):
            return None
        p1, p2 = intersect_circles(circle1, circle2)
        assert p1 == p2
        return p1

        
    
    def lft_GL2Z(self: typing.Self, matrix: gl2z.GL2Z) -> typing.Self:
        """
        Linear fractional transformation.
        """ 
        source_fraction = self.reduced_fraction
        alpha, beta, gamma, delta = matrix.alpha, matrix.beta, matrix.gamma, matrix.delta
        target_fraction = (alpha * source_fraction + beta) / (gamma * source_fraction + delta)
        return type(self)(target_fraction)


if __name__ == "__main__":
    x1 = quadratic_fields.RealQuadraticCompositum(1, {frozenset([2]): quadratic_fields.PureQuadraticSurd(2, 3)}) # 1 + 3√2
    x2 = quadratic_fields.RealQuadraticCompositum(1, {frozenset([3]): quadratic_fields.PureQuadraticSurd(3, 5)}) # 1 + 5√3
    c1 = Circle(2, 3, 4)
    c2 = Circle(1, -1, 6)
    p1, p2 = intersect_circles(c1, c2)
    assert power_of_point_wrt_circle(p1, c1) == power_of_point_wrt_circle(p1, c2)
    assert power_of_point_wrt_circle(p2, c1) == power_of_point_wrt_circle(p2, c2)
    assert p1 in RadicalLine(c1, c2) and p2 in RadicalLine(c1, c2)

    print("Two-point :", *intersect_circles(Circle(0, 0, 1), Circle(0, 1, 1)))
    print("Two-point:", *intersect_circles(Circle(Fraction(3, 2), Fraction(-7, 3), 5), Circle(4, 0, 3)))
    print("Tangent :", *intersect_circles(Circle(0, 0, Fraction(1, 2)), Circle(1, 0, Fraction(1, 2))))
    print("Disjoint:", intersect_circles(Circle(0, 0, 2), Circle(10, 0, 2)))

    circle1 = FordCircle(Fraction(1, 2))
    circle2 = FordCircle(Fraction(2, 3))
    assert FordCircle.are_tangent(circle1, circle2) == True

    circle1 = FordCircle(Fraction(1, 2))
    m = gl2z.GL2Z(3, 2, 2, 1)
    circle2 = circle1.lft_GL2Z(m)
    print(circle1, circle2)

    n = 23
    offset = 9
    fraction1, fraction0, fraction2 = farey.farey_list(n)[offset: offset + 3]
    print(FordCircle(fraction1), FordCircle(fraction0), FordCircle(fraction2))

    circle1 = FordCircle(Fraction(1, 2))
    circle2 = FordCircle(Fraction(2, 3))
    print(FordCircle.try_find_point_of_tangency(circle1, circle2))