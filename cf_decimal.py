import decimal

def decimal_continued_fraction(x: decimal.Decimal, max_terms=20):
    """Compute the first `max_terms` partial quotients of the real number x."""

    cf = []
    for _ in range(max_terms):
        a = int(x)          # floor
        cf.append(a)
        frac_part = x - a   # x - floor(x)
        if frac_part == 0:
            # x is an integer, so the continued fraction terminates
            break
        x = 1 / frac_part   # 1 / (x - floor(x))
    return cf


# Example usage:
if __name__ == "__main__":
    # Increase precision to 50 decimal places
    decimal.getcontext().prec = 50
    
    # Example 1: sqrt(7)
    sqrt7 = decimal.Decimal(7).sqrt()
    cf_sqrt7 = decimal_continued_fraction(sqrt7, max_terms=30)
    print("Partial quotients for sqrt(7):", cf_sqrt7)
    
    # Partial quotients for sqrt(7): [2, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1]
