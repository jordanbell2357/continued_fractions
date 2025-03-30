import decimal
from fractions import Fraction

def euclid(x: int, y: int) -> tuple[list[int], list[int], list[int], list[int]]:
    r0, r = x, y
    
    # s0, s, t0, t now set so that:
    #   r0 = x = 1*x - 0*y
    #   r  = y = 0*x - (-1)*y
    s0, s = 1, 0
    t0, t = 0, -1  # Changed here: t is now -1

    q_list = [None, None]
    r_list = [r0]
    s_list = [s0]
    t_list = [t0]
    
    while r != 0:
        q = r0 // r
        q_list.append(q)
        
        # Standard extended Euclid iteration
        r0, r = r, r0 - q*r
        s0, s = s, s0 - q*s
        t0, t = t, t0 - q*t
        
        r_list.append(r0)
        s_list.append(s0)
        t_list.append(t0)

    r_list.append(r)
    s_list.append(s)
    t_list.append(t)
    
    return q_list, r_list, s_list, t_list

# Example usage
if __name__ == '__main__':
    x, y = 7921, 4050
    decimal.getcontext().prec = 20  # high-precision decimals

    print(f"{x=}", f"{y=}")
    print()
    
    q_list, r_list, s_list, t_list = euclid(x, y)

    # Same convergent construction as before:
    convergents = [Fraction(t, s) for (t, s) in zip(t_list[2:], s_list[2:])]

    print(f"{q_list=}")
    print(f"{r_list=}")
    print(f"{s_list=}")
    print(f"{t_list=}")
    print()

    # The gcd is the next-to-last remainder
    r = r_list[-2]
    print(f"gcd of {x} and {y}:\t", r)
    print()

    # s, t are the next-to-last 's' and 't' values
    s, t = s_list[-2], t_list[-2]
    print(f"s * x - t * y = {s} * {x} - {t} * {y} =", s * x - t * y)
    print()

    cf = q_list[2:]
    print(f"Continued fraction for {x} / {y}:")
    print(cf)
    print()

    print(f"Convergents of continued fraction:")
    for convergent in convergents:
        p, q = convergent.numerator, convergent.denominator
        print(f"Decimal expansion for {p} / {q}:\t", decimal.Decimal(p) / decimal.Decimal(q))
    print()

    print(f"Decimal expansion for {x} / {y}:\t", decimal.Decimal(x) / decimal.Decimal(y))
