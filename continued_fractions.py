import itertools as it
import functools as ft
import operator
from fractions import Fraction
import math


def fraction_tuple_to_cf(fraction_tuple: tuple[int, int]) -> list[int]:
    x = Fraction(*fraction_tuple)
    cf = []
    while True:
        a = int(x) # floor
        cf.append(a)
        frac_part = x - a # x - floor(x)
        if frac_part == 0:
            # x is an integer, so the continued fraction terminates
            break
        x = 1 / frac_part # 1 / (x - floor(x))
    return cf

def cf_to_convergent_list(cf: list[int]) -> list[tuple[int, int]]:
    N = len(cf)
    if N == 0:
        return []
    p_prev2, p_prev1 = 0, 1
    q_prev2, q_prev1 = 1, 0
    convergent_list = []
    for k in range(N):
        p_k = cf[k] * p_prev1 + p_prev2
        q_k = cf[k] * q_prev1 + q_prev2
        convergent_list.append((p_k, q_k))
        p_prev2, p_prev1 = p_prev1, p_k
        q_prev2, q_prev1 = q_prev1, q_k
    return convergent_list

def stern_diatomic(n):
    def stern_diatomic_recurse(n):
        if n == 0:
            return 0
        elif n == 1:
            return 1
        else:
            if n % 2 == 0:
                return stern_diatomic_recurse(n // 2)
            else:
                m = (n - 1) // 2
                return stern_diatomic_recurse(m) + stern_diatomic_recurse(m + 1)
    return stern_diatomic_recurse(n)

def mediant(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (x[0] + y[0], x[1] + y[1])

class CalkinWilf:
    ROOT_FRACTION_TUPLE = (1, 1)
    ROOT_MOVE_LIST = []

    @staticmethod
    def bfs_index_to_fraction_tuple_recursive(bfs_index):
        @ft.lru_cache
        def recurse_q(n):
            if n == 0:
                return Fraction(1, 1)
            else:
                return 1 / (2 * int(recurse_q(n - 1)) - recurse_q(n - 1) + 1)
        q = recurse_q(bfs_index - 1)
        return q.as_integer_ratio()
    
    @classmethod
    def bfs_index_to_fraction_tuple(cls, bfs_index) -> tuple[int, int]:
        bit_string = f"{bfs_index:b}"
        runs = [(bit, len(list(g))) for bit, g in it.groupby(reversed(bit_string))]
        move_label_tuple, run_length_tuple = zip(*runs)
        move_label_list, run_length_list = list(move_label_tuple), list(run_length_tuple)
        if move_label_list[0] == "0":
            run_length_list.insert(0, 0)
            move_label_list.insert(0, "1")
        convergent_list = cf_to_convergent_list(run_length_list)
        return convergent_list[-1]

    @staticmethod
    def move_list_to_position(move_list: list[str]) -> int:
        if not move_list:
            return 0
        move_encode_dict = {"L": "0", "R": "1"}
        move_list_encoded = [move_encode_dict[move] for move in move_list]
        move_list_string = "".join(move_list_encoded)
        position = int(move_list_string, base=2)
        return position
    
    @classmethod
    def move_list_to_bfs_index(cls, move_list: list[str]) -> int:
        level = len(move_list)
        position = cls.move_list_to_position(move_list)
        bfs_index = 2 ** level + position
        return bfs_index
    
    @classmethod
    def move_list_to_fraction_tuple(cls, move_list):
        bfs_index = cls.move_list_to_bfs_index(move_list)
        fraction_tuple = cls.bfs_index_to_fraction_tuple(bfs_index)
        return fraction_tuple
    
    @classmethod
    def fraction_tuple_to_move_list(cls, fraction_tuple):
        def fraction_tuple_to_move_list_recurse(cls, fraction_tuple, move_list):
            if fraction_tuple == cls.ROOT_FRACTION_TUPLE:
                return move_list
            else:
                a, b = fraction_tuple
                if a < b:
                    new_move_list = ["L"] + move_list
                    new_fraction_tuple = (a, b - a)
                else:
                    new_move_list = ["R"] + move_list
                    new_fraction_tuple = (a - b, b)
                return fraction_tuple_to_move_list_recurse(cls, new_fraction_tuple, new_move_list)
        return fraction_tuple_to_move_list_recurse(cls, fraction_tuple, [])
    
    @classmethod
    def fraction_tuple_to_bfs_index(cls, fraction_tuple):
        move_list = cls.fraction_tuple_to_move_list(fraction_tuple)
        bfs_index = cls.move_list_to_bfs_index(move_list)
        return bfs_index

    @classmethod
    def bfs_index_to_move_list(cls, bfs_index: int) -> list[str]:
        fraction_tuple = cls.bfs_index_to_fraction_tuple(bfs_index)
        move_list = cls.fraction_tuple_to_move_list(fraction_tuple)
        return move_list
    
    @classmethod
    def move_list_to_node(cls, move_list: list[str]):
        return cls(move_list)        
    
    @classmethod
    def bfs_index_to_node(cls, bfs_index: int):
        move_list = cls.bfs_index_to_move_list(bfs_index)
        return cls(move_list)

    @classmethod
    def fraction_tuple_to_node(cls, fraction_tuple):
        move_list = cls.fraction_tuple_to_move_list(fraction_tuple)
        return cls(move_list)
    
    @classmethod
    def depth_to_bfs_node_list(cls, depth: int):
        cw_tree = []
        N = cls()
        level_node_list = []
        level_node_list.append(N)
        cw_tree.append(level_node_list)
        for n in range(depth):
            previous_level_node_list = level_node_list
            level_node_list = []
            for N in previous_level_node_list:
                level_node_list.append(N.L())
                level_node_list.append(N.R())
            cw_tree.append(level_node_list)
        cw_bfs_list = ft.reduce(operator.iadd, cw_tree, [])
        return cw_bfs_list

    def __init__(self, move_list=None):
        if move_list is None:
            self.move_list = type(self).ROOT_MOVE_LIST
            self.move_string = "".join(self.move_list)
            self.level = len(self.move_list)
            self.position = type(self).move_list_to_position(self.move_list)
            self.bfs_index = 2 ** self.level + self.position
            self.fraction_tuple = type(self).ROOT_FRACTION_TUPLE
            self.cf = fraction_tuple_to_cf(self.fraction_tuple)
        else:
            self.move_list = move_list
            self.move_string = "".join(self.move_list)
            self.level = len(self.move_list)
            self.position = type(self).move_list_to_position(self.move_list)
            self.bfs_index = 2 ** self.level + self.position
            self.fraction_tuple = type(self).move_list_to_fraction_tuple(self.move_list)
            self.cf = fraction_tuple_to_cf(self.fraction_tuple)

    def __repr__(self):
        return f"{type(self).__name__}({self.move_list})"

    def __str__(self):
        return f"Fraction:\t{self.fraction_tuple}\nCW Moves:\t{self.move_string}\nCW BFS index:\t{self.bfs_index}\nCF:\t\t{self.cf}"

    def __eq__(self, other):
        return self.move_list == other.move_list
    
    def __lt__(self, other):
        return len(self.move_list) < len(other.move_list) and other.move_list[:len(self.move_list)] == self.move_list
    
    def L(self):
        new_move_list = self.move_list + ["L"]
        return type(self)(new_move_list)
    
    def R(self):
        new_move_list = self.move_list + ["R"]
        return type(self)(new_move_list)
    
    def P(self):
        move_list = self.move_list
        if not move_list:
            return None
        else:
            new_move_list = move_list[0:-1]
            return type(self)(new_move_list)
    
    def run_tuple_list(self):
        return [(k, len(list(g))) for k, g in it.groupby(self.move_list)]


class SternBrocot:
    ROOT_LEFT_TUPLE, ROOT_MIDDLE_TUPLE, ROOT_RIGHT_TUPLE = (0, 1), (1, 1), (1, 0)

    @classmethod
    def fraction_tuple_to_node(cls, fraction_tuple):
        N = cls([], cls.ROOT_LEFT_TUPLE, cls.ROOT_MIDDLE_TUPLE, cls.ROOT_RIGHT_TUPLE)
        while True:
            fraction_tuple_value = fraction_tuple[0] / fraction_tuple[1]
            node_middle_tuple_fraction_value = N.middle_tuple[0] / N.middle_tuple[1]
            if fraction_tuple_value < node_middle_tuple_fraction_value:
                N = N.L()
            elif fraction_tuple_value > node_middle_tuple_fraction_value:
                N = N.R()
            else:
                return N

    def __init__(self, move_list, left_tuple, middle_tuple, right_tuple):
        self.move_list = move_list
        self.left_tuple = left_tuple
        self.middle_tuple = middle_tuple
        self.right_tuple = right_tuple
        self.fraction_tuple = self.middle_tuple

    def L(self):
        new_move_list = self.move_list + ['L']
        new_middle_tuple = mediant(self.left_tuple, self.middle_tuple)
        return type(self)(new_move_list, self.left_tuple, new_middle_tuple, self.middle_tuple)

    def R(self):
        new_move_list = self.move_list + ['R']
        new_middle_tuple = mediant(self.middle_tuple, self.right_tuple)
        return type(self)(new_move_list, self.middle_tuple, new_middle_tuple, self.right_tuple)
    
    def run_list(self):
        return [(k, len(list(g))) for k, g in it.groupby(self.move_list)]
    
    def cf(self):
        run_list = self.run_list()
        if run_list[0][0] == 'L':
            run_list = [('R', 0)] + run_list
        cf = [run[1] for run in run_list]
        return cf
    
    def __eq__(self, other):
        return self.move_list == other.move_list


if __name__ == "__main__":
    bfs_index = 71

    assert CalkinWilf.bfs_index_to_fraction_tuple_recursive(bfs_index) == \
        CalkinWilf.bfs_index_to_fraction_tuple(bfs_index)

    assert (stern_diatomic(bfs_index), stern_diatomic(bfs_index + 1)) == \
        CalkinWilf.bfs_index_to_fraction_tuple_recursive(bfs_index)
    
    fraction_tuple = CalkinWilf.bfs_index_to_fraction_tuple(bfs_index)
    move_list = CalkinWilf.fraction_tuple_to_move_list(fraction_tuple)
    assert CalkinWilf.bfs_index_to_move_list(bfs_index) == move_list

    N1 = CalkinWilf()
    N2 = CalkinWilf([])
    assert N1 == N2

    depth = int(math.log2(bfs_index))
    cw_bfs_list = CalkinWilf.depth_to_bfs_node_list(depth)
    N1 = CalkinWilf.bfs_index_to_node(bfs_index)
    N2 = cw_bfs_list[bfs_index - 1]
    assert N1 == N2

    N = CalkinWilf.bfs_index_to_node(bfs_index)
    assert N.bfs_index == bfs_index

    N1 = CalkinWilf.bfs_index_to_node(bfs_index)
    N2 = CalkinWilf(N1.move_list)
    assert N1 == N2

    N = CalkinWilf.bfs_index_to_node(bfs_index)
    N1, N2 = N.L(), N.R()
    assert N < N1
    assert N < N2
    assert not (N1 > N2)
    assert not (N2 > N1)

    N = CalkinWilf.bfs_index_to_node(bfs_index)
    assert CalkinWilf.move_list_to_bfs_index(N.move_list) == bfs_index

    N = CalkinWilf.bfs_index_to_node(bfs_index)
    assert CalkinWilf.move_list_to_fraction_tuple(N.move_list) == N.fraction_tuple

    fraction_tuple = (21, 13)

    # N = CalkinWilf.fraction_tuple_to_node(fraction_tuple)
    # print(N)

if __name__ == "__main__":
    fraction_tuple = (23, 196)

    N = SternBrocot.fraction_tuple_to_node(fraction_tuple)
    print(N.fraction_tuple)
