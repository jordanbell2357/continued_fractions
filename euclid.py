import decimal
from fractions import Fraction

def euclid(x: int, y: int) -> tuple[list[int], list[int], list[int], list[int]]:
    r0, r = x, y
    s0, s = 1, 0
    t0, t = 0, 1
    
    q_list = [None, None]
    r_list = [r0]
    s_list = [s0]
    t_list = [t0]
    
    while r != 0:
        # Compute the quotient
        q = r0 // r
        
        # Append q to list
        q_list.append(q)
        
        # Iterate
        r0, r = r, r0 - q * r
        s0, s = s, s0 - q * s
        t0, t = t, t0 - q * t
        
        # Append the new r0, s0, t0 values to lists
        r_list.append(r0)
        s_list.append(s0)
        t_list.append(t0)

    r_list.append(r)
    s_list.append(s)
    t_list.append(t)
    
    # Return the lists of q, r, s, t
    return q_list, r_list, s_list, t_list


def get_convergents(cf) -> list[Fraction]:
    N = len(cf)
    if N == 0:
        return []
    
    p_prev2, p_prev1 = 0, 1
    q_prev2, q_prev1 = 1, 0
    convergent_list = []
    
    for k in range(N):
        p_k = cf[k] * p_prev1 + p_prev2
        q_k = cf[k] * q_prev1 + q_prev2
        convergent_list.append(Fraction(p_k, q_k))

        p_prev2, p_prev1 = p_prev1, p_k
        q_prev2, q_prev1 = q_prev1, q_k
    
    return convergent_list
    

# Example usage:
if __name__ == '__main__':
    x, y = 240, 46

    decimal.getcontext().prec = 20

    print(f"{x=}", f"{y=}")

    q_list, r_list, s_list, t_list = euclid(x, y)

    print(f"{q_list=}")
    print(f"{r_list=}")
    print(f"{s_list=}")
    print(f"{t_list=}")

    print()

    r = r_list[-2]
    print(f"gcd of {x} and {y}:\t", r)

    print()

    s, t = s_list[-2], t_list[-2]

    print(f"s * x + t * y = {s} * {x} + {t} * {y} =", s * x + t * y)

    print()

    cf = q_list[2:]
    print(f"Continued fraction for {x} / {y}:\t", cf)

    print()

    convergents = get_convergents(cf)
    print(f"Convergents of continued fraction:")
    for convergent in convergents:
        p, q = convergent.numerator, convergent.denominator
        print(f"Decimal expansion for {p} / {q}:\t", decimal.Decimal(p) / decimal.Decimal(q))

    print()

    print(f"Decimal expansion for {x} / {y}:\t", decimal.Decimal(x) / decimal.Decimal(y))
