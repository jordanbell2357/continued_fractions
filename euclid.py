import decimal

def euclid(x: int, y: int) -> tuple[list[int], list[int], list[int], list[int]]:
    r0, r = x, y
    
    #   r0 = x = 1 * x - 0 * y
    #   r  = y = 0 * x - (-1) * y

    s0, s = 1, 0
    t0, t = 0, 1

    q_list = [None, None]
    r_list = [r0]
    s_list = [s0]
    t_list = [t0]
    
    while r != 0:
        q = r0 // r
        q_list.append(q)
        
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
    x, y = 71755875, 61735500
    decimal.getcontext().prec = 20

    print(f"{x=}", f"{y=}")
    print()
    
    q_list, r_list, s_list, t_list = euclid(x, y)
    print(f"{q_list=}")
    print(f"{r_list=}")
    print(f"{s_list=}")
    print(f"{t_list=}")
    print()

    r = r_list[-2]
    print(f"gcd({x}, {y}) = {r}")
    print()

    s, t = s_list[-2], t_list[-2]
    print("Bézout coefficients")
    print(f"s * {x} + t * {y} = {r}")
    print(f"{s=}, {t=}")
    print()

    cf = q_list[2:]
    print(f"Continued fraction of {x} / {y}")
    print(cf)
    print()

    print(f"Convergents")
    x_list = [(-1) ** (k + 1) * t for k, t in enumerate(t_list)]
    y_list = [(-1) ** k * s for k, s in enumerate(s_list)]
    print(f"{x_list=}")
    print(f"{y_list=}")
    print()

    print("Decimal expansions")
    for k, (xk, yk) in enumerate(zip(x_list[2:], y_list[2:])):
        print(f"x_{k} / y_{k}\t= {decimal.Decimal(xk) / decimal.Decimal(yk)}")
    print(f"x / y\t\t= {decimal.Decimal(x) / decimal.Decimal(y)}")
