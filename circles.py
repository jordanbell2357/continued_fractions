import math
from numbers import Rational
from numbers import Number
from fractions import Fraction
from dataclasses import dataclass
from collections import abc
import typing

import quadratic_fields



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
    def __init__(self, c1: Circle, c2: Circle) -> None:
        self.c1 = c1
        self.c2 = c2

    def __contains__(self, point_P: Point) -> bool:
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


# -------------------------------------------------------------------------
# Quick examples
# -------------------------------------------------------------------------
if __name__ == "__main__":
    x1 = quadratic_fields.RealQuadraticCompositum(1, {frozenset([2]): quadratic_fields.PureQuadraticSurd(2, 3)})   # 1 + 3√2
    x2 = quadratic_fields.RealQuadraticCompositum(1, {frozenset([3]): quadratic_fields.PureQuadraticSurd(3, 5)})   # 1 + 5√3
    p1 = Point(x1, x2)
    c1 = Circle(2, 3, 4)
    c2 = Circle(1, -1, 6)
    print(power_of_point_wrt_circle(p1, c1))

    p1, p2 = intersect_circles(Circle(0, 0, 5), Circle(4, 0, 3))
    print(p1, p2)
    print("Two-point:", *intersect_circles(Circle(Fraction(3, 2), Fraction(-7, 3), 5), Circle(4, 0, 3)))
    print("Two-point :", *intersect_circles(Circle(0, 0, 1), Circle(0, 1, 1)))
    print("Tangent :", *intersect_circles(Circle(0, 0, 5), Circle(10, 0, 5)))
    print("Tangent :", *intersect_circles(Circle(0, 0, Fraction(1, 2)), Circle(1, 0, Fraction(1, 2))))
    print("Disjoint:", intersect_circles(Circle(0, 0, 2), Circle(10, 0, 2)))
    c1, c2 = Circle(0, 0, Fraction(1, 2)), Circle(1, 0, Fraction(1, 2))
    p = intersect_circles(c1, c2)[0]
    print(p in RadicalLine(c1, c2))



