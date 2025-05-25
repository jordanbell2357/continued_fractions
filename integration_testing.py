import prime_numbers
import partitions
import cflib
import farey
import cw
import sb

if __name__ == "__main__":
    # Example: Mertens function
    n = 20

    f = farey.Farey(n)
    assert f.mertens_function == prime_numbers.mertens(n)

    assert farey.mertens_function(n) == prime_numbers.mertens(n)


    # Example: recurrence formula for partition function
    n = 20
    assert partitions.partition_function(n) == sum(prime_numbers.sum_of_divisors(n - k) * partitions.partition_function(k) for k in range(n)) // n