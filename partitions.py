from array import array
from collections import abc

import prime_numbers


def make_partition_list_generator(n: int, max_part: int) -> abc.Generator[list[int]]:
    if n == 0:
        yield []
    else:
        # condition on all cases of max_part
        for k in range(min(n, max_part), 0, -1):
            # find number of partitions of n - k with max_part == k
            for tail in make_partition_list_generator(n - k, k):
                yield [k] + tail


def list_partitions(n: int) -> list[list[int]]:
    return list(make_partition_list_generator(n, n))


def partition_function_array(n: int) -> int:
    p = [0] * (n + 1)
    p[0] = 1
    for i in range(1, n + 1):
        for s in range(i, n + 1):
            p[s] += p[s - i]
    return p[n]


def partition_function(n: int) -> int:
    p = array("L", [0] * (n + 1))
    p[0] = 1
    for i in range(1, n + 1):
        for s in range(i, n + 1):
            p[s] += p[s - i]
    return p[n]


if __name__ == "__main__":
    n = 16
    assert len(list_partitions(n)) == partition_function(n)

    # Recurrence formula for partition function
    n = 20
    assert partition_function(n) == sum(prime_numbers.sum_of_divisors(n - k) * partition_function(k) for k in range(n)) // n

    n = 7
    partition_list_generator = make_partition_list_generator(n, n)
    for partition_list in partition_list_generator:
        print(partition_list)