from fractions import Fraction
import typing
import math

def euclid(x: int, y: int) -> tuple[list[int], list[int], list[int], list[int]]:
    r0 = 1 * x - 0 * y
    r  = 0 * x - (-1) * y

    a0, a = 1, 0
    b0, b = 0, 1

    q_list = []
    r_list = [r0]
    a_list = [a0]
    b_list = [b0]
    
    while r != 0:
        q = r0 // r
        q_list.append(q)
        
        r0, r = r, r0 - q * r
        a0, a = a, a0 - q * a
        b0, b = b, b0 - q * b
        
        r_list.append(r0)
        a_list.append(a0)
        b_list.append(b0)

    r_list.append(r)
    a_list.append(a)
    b_list.append(b)
    
    return q_list, r_list, a_list, b_list


class EEA:
    def __init__(self: typing.Self, x: int, y: int) -> None:
        self.x = x
        self.y = y
        q_list, r_list, a_list, b_list = euclid(self.x, self.y)
        self.q_list = q_list
        self.r_list = r_list
        self.a_list = a_list
        self.b_list = b_list
        self.x_list = [(-1) ** (k + 1) * t for k, t in enumerate(self.b_list)]
        self.y_list = [(-1) ** k * s for k, s in enumerate(self.a_list)]
        self.gcd = self.r_list[-2]
        self.bezout_x = self.a_list[-2]
        self.bezout_y = self.b_list[-2]
        self.convergent_list = list(zip(self.x_list, self.y_list))

    def __repr__(self) -> str:
        return f"{type(self).__name__}(x={self.x}, y={self.y})"
    
    def __len__(self):
        return len(self.q_list)


# Example usage
if __name__ == '__main__':
    x, y = 240, 18

    eea = EEA(x, y)
    q_list, r_list, a_list, b_list = euclid(x, y)
    assert eea.q_list == q_list
    assert eea.r_list == r_list
    assert eea.a_list == a_list
    assert eea.b_list == b_list

    eea = EEA(x, y)
    assert len(eea) == len(eea.q_list)

    eea = EEA(x, y)
    assert eea.gcd == math.gcd(x, y)

    eea = EEA(x, y)
    assert eea.bezout_x * x + eea.bezout_y * y == math.gcd(x, y)

    eea = EEA(x, y)
    assert eea.convergent_list[-1] == Fraction(x, y).as_integer_ratio()
