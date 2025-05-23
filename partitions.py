from array import array
from collections import abc


def partition_list_generator(n: int, max_part: int) -> abc.Generator[list[int]]:
    if n == 0:
        yield []
    else:
        # condition on all cases of max_part
        for k in range(min(n, max_part), 0, -1):
            # find number of partitions of n - k with max_part == k
            for tail in partition_list_generator(n - k, k):
                yield [k] + tail


def list_partitions(n: int) -> list[list[int]]:
    return list(partition_list_generator(n=n, max_part=n))


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
    for n in range(1, 16):
        print(n, partition_function(n))

    n = 20

    assert len(list_partitions(n)) == partition_function(n)
