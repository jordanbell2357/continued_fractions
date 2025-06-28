import math
from numbers import Rational
from numbers import Number
from fractions import Fraction
from dataclasses import dataclass
from collections import abc
import typing

import prime_numbers



@dataclass
class Point:
    x: Number
    y: Number


@dataclass
class Circle:
    a: Number
    b: Number
    r: Number


def power_of_point_wrt_circle(point_P: Point, circle_O: Circle) -> Number:
    dist_P_O = math.dist((point_P.x, point_P.y), (circle_O.a, circle_O.b))
    return dist_P_O ** 2 - circle_O.r ** 2


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
    d_c1_c2 = math.sqrt(d_c1_c2_squared)

    # Step 2 – distance d₁ of foot D from C₁
    d_foot_c1 = (d_c1_c2_squared + r1 * r1 - r2 * r2) / (2 * d_c1_c2)

    # Step 3 – half-chord length h
    h_squared = r1 * r1 - d_foot_c1 * d_foot_c1
    if h_squared < 0: # no intersection in reals
        return None
    h = math.sqrt(h_squared)

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
    print("Two-point:", intersect_circles(Circle(0, 0, 5), Circle(4, 0, 3)))
    print("Tangent :", intersect_circles(Circle(0, 0, 5), Circle(10, 0, 5)))
    print("Disjoint:", intersect_circles(Circle(0, 0, 2), Circle(10, 0, 2)))

