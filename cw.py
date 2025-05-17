import itertools as it
import functools as ft
import operator
from fractions import Fraction
import math
import typing
from pprint import pprint

import cflib

class CalkinWilf:
    ROOT_FRACTION_TUPLE = (1, 1)
    ROOT_MOVE_LIST = []
    ROOT_POSITION = 0
    MOVE_ENCODE_DICT = {"L": "0", "R": "1"}

    @classmethod
    def bfs_index_to_fraction_tuple_recursive(cls, bfs_index: int) -> tuple[int, int]:
        @ft.lru_cache
        def recurse_q(n: int) -> Fraction:
            if n == 0:
                return Fraction(*cls.ROOT_FRACTION_TUPLE)
            else:
                return 1 / (2 * int(recurse_q(n - 1)) - recurse_q(n - 1) + 1)
        q = recurse_q(bfs_index - 1)
        return q.as_integer_ratio()
    
    @classmethod
    def bfs_index_to_fraction_tuple(cls, bfs_index: int) -> tuple[int, int]:
        bit_string = f"{bfs_index:b}"
        runs = [(bit, len(list(g))) for bit, g in it.groupby(reversed(bit_string))]
        move_label_tuple, run_length_tuple = zip(*runs)
        move_label_list, run_length_list = list(move_label_tuple), list(run_length_tuple)
        if move_label_list[0] == "0":
            run_length_list.insert(0, 0)
            move_label_list.insert(0, "1")
        return cflib.cf_to_fraction_tuple(run_length_list)

    @classmethod
    def move_list_to_position(cls, move_list: list[str]) -> int:
        if not move_list:
            return cls.ROOT_POSITION
        move_list_encoded = [cls.MOVE_ENCODE_DICT[move] for move in move_list]
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
    def move_list_to_fraction_tuple(cls, move_list: list[str]) -> tuple[int, int]:
        bfs_index = cls.move_list_to_bfs_index(move_list)
        fraction_tuple = cls.bfs_index_to_fraction_tuple(bfs_index)
        return fraction_tuple
    
    @classmethod
    def fraction_tuple_to_move_list(cls, fraction_tuple: tuple[int, int]) -> list[str]:
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
    def fraction_tuple_to_bfs_index(cls, fraction_tuple: tuple[int, int]) -> int:
        move_list = cls.fraction_tuple_to_move_list(fraction_tuple)
        bfs_index = cls.move_list_to_bfs_index(move_list)
        return bfs_index

    @classmethod
    def bfs_index_to_move_list(cls, bfs_index: int) -> list[str]:
        fraction_tuple = cls.bfs_index_to_fraction_tuple(bfs_index)
        move_list = cls.fraction_tuple_to_move_list(fraction_tuple)
        return move_list
    
    @classmethod
    def move_list_to_node(cls, move_list: list[str]) -> typing.Self:
        return cls(move_list)        
    
    @classmethod
    def bfs_index_to_node(cls, bfs_index: int) -> typing.Self:
        move_list = cls.bfs_index_to_move_list(bfs_index)
        return cls(move_list)

    @classmethod
    def fraction_tuple_to_node(cls, fraction_tuple) -> typing.Self:
        move_list = cls.fraction_tuple_to_move_list(fraction_tuple)
        return cls(move_list)
    
    @classmethod
    def depth_to_cw_tree(cls, depth: int) -> list[list[typing.Self]]:
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
        return cw_tree
    
    @classmethod
    def depth_to_bfs_node_list(cls, depth: int) -> list[typing.Self]:
        cw_tree = cls.depth_to_cw_tree(depth)
        cw_bfs_list = ft.reduce(operator.iadd, cw_tree, [])
        return cw_bfs_list


    def __init__(self, move_list=None) -> None:
        if move_list is None:
            self.move_list = type(self).ROOT_MOVE_LIST
            self.move_string = "".join(self.move_list)
            self.level = len(self.move_list)
            self.position = type(self).ROOT_POSITION
            self.bfs_index = 2 ** self.level + self.position
            self.fraction_tuple = type(self).ROOT_FRACTION_TUPLE
            self.cf = cflib.fraction_tuple_to_cf(self.fraction_tuple)
        else:
            self.move_list = move_list
            self.move_string = "".join(self.move_list)
            self.level = len(self.move_list)
            self.position = type(self).move_list_to_position(self.move_list)
            self.bfs_index = 2 ** self.level + self.position
            self.fraction_tuple = type(self).move_list_to_fraction_tuple(self.move_list)
            self.cf = cflib.fraction_tuple_to_cf(self.fraction_tuple)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.move_list})"

    def __str__(self) -> str:
        return f"Fraction:\t{self.fraction_tuple}\nCF:\t\t{self.cf}\nCW Moves:\t{self.move_string}\nCW BFS index:\t{self.bfs_index}"

    def __eq__(self, other: typing.Self) -> bool:
        return self.move_list == other.move_list
    
    def __lt__(self, other: typing.Self) -> bool:
        return len(self.move_list) < len(other.move_list) and other.move_list[0:len(self.move_list)] == self.move_list
    
    def __len__(self) -> int:
        return len(self.move_list)
    
    def L(self) -> typing.Self:
        new_move_list = self.move_list + ["L"]
        return type(self)(new_move_list)
    
    def R(self) -> typing.Self:
        new_move_list = self.move_list + ["R"]
        return type(self)(new_move_list)
    
    def P(self) -> typing.Self:
        move_list = self.move_list
        if not move_list:
            return None
        else:
            new_move_list = move_list[0:-1]
            return type(self)(new_move_list)
    
    def run_tuple_list(self) -> list[tuple[str, int]]:
        return [(k, len(list(g))) for k, g in it.groupby(self.move_list)]
    

class CalkinWilfTree:
    def __init__(self, depth: int) -> None:
        self.depth = depth
        self.cw_tree = CalkinWilf.depth_to_cw_tree(self.depth)
        self.bfs_node_list = CalkinWilf.depth_to_bfs_node_list(self.depth)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(depth={self.depth})"
    
    def __eq__(self, other: typing.Self) -> bool:
        return self.depth == other.depth
    
    def __lt__(self, other: typing.Self) -> bool:
        return self.depth < other.depth
    
    def __iter__(self):
        return iter(self.cw_tree)


if __name__ == "__main__":
    bfs_index = 71

    assert CalkinWilf.bfs_index_to_fraction_tuple_recursive(bfs_index) == \
        CalkinWilf.bfs_index_to_fraction_tuple(bfs_index)

    assert (cflib.stern_diatomic(bfs_index), cflib.stern_diatomic(bfs_index + 1)) == \
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

    cw_tree = CalkinWilfTree(3)
    pprint(cw_tree.bfs_node_list)