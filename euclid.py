def euclid(x, y):
    r0, r = x, y
    
    #   r0 = 1 * x - 0 * y
    #   r  = 0 * x - (-1) * y

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


# Example usage
if __name__ == '__main__':
    x, y = 240, 46

    print(f"{x=}, {y=}")
    print()
    
    q_list, r_list, a_list, b_list = euclid(x, y)
    print(f"{q_list=}")
    print(f"{r_list=}")
    print(f"{a_list=}")
    print(f"{b_list=}")
    print()

    print(80 * "-")
    print()

    print("Greatest common divisor")
    print(f"gcd({x}, {y}) = {r_list[-2] = }")
    print()

    print("Bézout coefficients ax + by = gcd(x, y)")
    print(f"a * {x} + b * {y} = {r_list[-2]}")
    print(f"a = {a_list[-2] = }")
    print(f"b = {b_list[-2] = }")
    print()

    print(f"Continued fraction of {x} / {y}")
    print(f"{q_list = }")
    print()

    print(f"Convergents of continued fraction")
    x_list = [(-1) ** (k + 1) * t for k, t in enumerate(b_list)]
    y_list = [(-1) ** k * s for k, s in enumerate(a_list)]
    convergent_list = list(zip(x_list, y_list))
    print(f"{convergent_list = }")