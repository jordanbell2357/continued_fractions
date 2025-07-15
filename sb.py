import itertools as it
from fractions import Fraction
import typing
import random
import functools as ft
import operator
from collections import abc
import textwrap
import heapq
import copy

import cflib


def mediant(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (x[0] + y[0], x[1] + y[1])


class SternBrocot(abc.Sequence):
    """"
    cf. https://mattbaker.blog/2019/01/28/the-stern-brocot-tree-hurwitzs-theorem-and-the-markoff-uniqueness-conjecture/
    and https://rosettacode.org/wiki/Stern-Brocot_sequence    
    """

    ALPHABET = ["L", "R"]
    ROOT_LEFT_TUPLE, ROOT_MEDIANT_TUPLE, ROOT_RIGHT_TUPLE = (0, 1), (1, 1), (1, 0)
    ROOT_MOVE_LIST = []

    @staticmethod
    def sb_left_child_parent_mediant(x: tuple[int, int], m: tuple[int, int]) -> tuple[int, int]:
        return (m[0] - x[0], m[1] - x[1])

    @staticmethod
    def sb_right_child_parent_mediant(y: tuple[int, int], m: tuple[int, int]) -> tuple[int, int]:
        return (m[0] - y[0], m[1] - y[1])
    
    @staticmethod
    def move_string_to_move_list(move_string: str) -> list[str]:
        move_list = [move for move in move_string]
        return move_list

    @classmethod
    def fraction_tuple_to_node(cls, fraction_tuple: tuple[int, int]) -> typing.Self:
        N = cls()
        while True:
            fraction_tuple_value = fraction_tuple[0] / fraction_tuple[1]
            node_mediant_tuple_fraction_value = N.mediant_tuple[0] / N.mediant_tuple[1]
            if fraction_tuple_value < node_mediant_tuple_fraction_value:
                N = N.L()
            elif fraction_tuple_value > node_mediant_tuple_fraction_value:
                N = N.R()
            else:
                return N

    @classmethod
    def move_list_to_node(cls, move_list: list[str]) -> typing.Self:
        N = cls()
        while True:
            if not move_list:
                return N
            elif move_list[0] == "L":
                N = N.L()
                move_list = move_list[1:]
            elif move_list[0] == "R":
                N = N.R()
                move_list = move_list[1:]
            else:
                raise Exception
            
    @classmethod
    def move_list_to_fraction_tuple(cls, move_list: list[str]) -> tuple[int, int]:
        N = cls.move_list_to_node(move_list)
        fraction_tuple = N.mediant_tuple
        return fraction_tuple
    
    @classmethod
    def move_list_to_cf(cls, move_list: list[str]) -> list[int]:
        fraction_tuple = cls.move_list_to_fraction_tuple(move_list)
        cf = cflib.fraction_tuple_to_cf(fraction_tuple)
        return cf
    
    @classmethod
    def level_to_nodes(cls, level: int) -> list[list[typing.Self]]:
        move_string_list = ["".join(moves) for moves in it.product(cls.ALPHABET, repeat=level)]
        move_list_list = [cls.move_string_to_move_list(move_string) for move_string in move_string_list]
        node_list = [cls.move_list_to_node(move_list) for move_list in move_list_list]
        return node_list
    
    @classmethod
    def depth_to_sb_tree(cls, depth: int) -> list[list[typing.Self]]:
        sb_tree = []
        for level in range(depth + 1):
            node_list = cls.level_to_nodes(level)
            sb_tree.append(node_list)
        return sb_tree

    @classmethod
    def depth_to_bfs_node_list(cls, depth: int) -> list[typing.Self]:
        sb_tree = cls.depth_to_sb_tree(depth)
        bfs_node_list = ft.reduce(operator.iadd, sb_tree, [])
        return bfs_node_list


    def __init__(self, move_list=None, left_tuple=None, mediant_tuple=None, right_tuple=None) -> None:
        if move_list is None:
            self.move_list = type(self).ROOT_MOVE_LIST
            self.left_tuple = type(self).ROOT_LEFT_TUPLE
            self.mediant_tuple = type(self).ROOT_MEDIANT_TUPLE
            self.right_tuple = type(self).ROOT_RIGHT_TUPLE
        else:
            self.move_list = move_list
            self.left_tuple = left_tuple
            self.mediant_tuple = mediant_tuple
            self.right_tuple = right_tuple
        self.level = len(self.move_list)
        self.cf = cflib.fraction_tuple_to_cf(self.mediant_tuple)
        self.move_string = "".join(self.move_list)
        self.mediant_string = str(self.mediant_tuple[0]) + "/" + str(self.mediant_tuple[1])
    
    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.move_list}, {self.left_tuple}, {self.mediant_tuple}, {self.right_tuple})"

    def __eq__(self, other: typing.Self) -> bool:
        return self.move_list == other.move_list
    
    def __lt__(self, other: typing.Self) -> bool:
        return len(self.move_list) < len(other.move_list) and other.move_list[0:len(self.move_list)] == self.move_list
    
    def __gt__(self, other: typing.Self) -> bool:
        return len(other.move_list) < len(self.move_list) and self.move_list[0:len(other.move_list)] == other.move_list
    
    def __abs__(self) -> Fraction:
        return Fraction(*self.mediant_tuple)
    
    def __iter__(self) -> abc.Iterator[str]:
        yield from self.move_list
    
    def __getitem__(self, index: int) -> str:
        return self.move_list[index]

    def __len__(self) -> int:
        return len(self.move_list)

    def __str__(self) -> str:
        return f"{self.move_list}, {self.mediant_tuple}"
    
    def L(self) -> typing.Self:
        new_move_list = self.move_list + ['L']
        new_mediant_tuple = mediant(self.left_tuple, self.mediant_tuple)
        return type(self)(new_move_list, self.left_tuple, new_mediant_tuple, self.mediant_tuple)

    def R(self) -> typing.Self:
        new_move_list = self.move_list + ['R']
        new_mediant_tuple = mediant(self.mediant_tuple, self.right_tuple)
        return type(self)(new_move_list, self.mediant_tuple, new_mediant_tuple, self.right_tuple)
    
    def P(self) -> typing.Self:
        move_list = self.move_list
        if not move_list:
            return None
        previous_move = move_list[-1]
        new_move_list = move_list[0:-1]
        left_tuple, mediant_tuple, right_tuple = self.left_tuple, self.mediant_tuple, self.right_tuple
        if previous_move == "L":
            new_left_tuple = left_tuple
            new_mediant_tuple = type(self).sb_left_child_parent_mediant(left_tuple, mediant_tuple)
            new_right_tuple = mediant_tuple
        elif previous_move == "R":
            new_left_tuple = mediant_tuple
            new_mediant_tuple = type(self).sb_right_child_parent_mediant(right_tuple, mediant_tuple)
            new_right_tuple = right_tuple
        return type(self)(new_move_list, new_left_tuple, new_mediant_tuple, new_right_tuple)
    

class SternBrocotTree(abc.Sequence):
    def __init__(self, depth: int) -> None:
        self.depth = depth
        self.sb_tree = SternBrocot.depth_to_sb_tree(self.depth)
        self.bfs_node_list = SternBrocot.depth_to_bfs_node_list(self.depth)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(depth={self.depth})"
    
    def __eq__(self, other: typing.Self) -> bool:
        return self.depth == other.depth

    def __lt__(self, other: typing.Self) -> bool:
        return self.depth < other.depth
    
    def __getitem__(self, index: int) -> SternBrocot:
        return self.bfs_node_list[index]
    
    def __len__(self) -> int:
        return len(self.bfs_node_list)
    
    def __iter__(self) -> abc.Iterator:
        yield from self.bfs_node_list
    
    def __str__(self) -> str:
        tree_string = "\n\n".join(["\n".join(textwrap.wrap("\t".join([node.mediant_string for node in level]))) for level in self.sb_tree])
        return tree_string


if __name__ == "__main__":
    # Example: random move list
    k = 5

    random.seed(444)
    move_list = random.choices(SternBrocot.ALPHABET, k=5)

    N = SternBrocot.move_list_to_node(move_list)
    assert eval(repr(N)) == N

    N = SternBrocot.move_list_to_node(move_list)
    assert SternBrocot.move_list_to_cf(move_list) == N.cf

    N = SternBrocot.move_list_to_node(move_list)
    assert SternBrocot.move_list_to_fraction_tuple(move_list) == N.mediant_tuple


    # Example: fraction tuple to node
    fraction_tuple = (3, 14)

    N = SternBrocot.fraction_tuple_to_node(fraction_tuple)
    assert N.mediant_tuple == fraction_tuple

    N = SternBrocot.fraction_tuple_to_node(fraction_tuple)
    assert N.cf == cflib.fraction_tuple_to_cf(fraction_tuple)    

    assert SternBrocot.move_string_to_move_list("") == []


    # Example: Stern-Brocot tree
    depth = 10

    sb_tree = SternBrocotTree(depth)
    assert eval(repr(sb_tree)) == sb_tree

    sb_tree = SternBrocotTree(depth)
    assert sb_tree.bfs_node_list == SternBrocot.depth_to_bfs_node_list(depth)

    sb_tree = SternBrocotTree(depth)
    assert len(sb_tree) == len(sb_tree.bfs_node_list)


    depth = 3
    sb_tree = SternBrocotTree(depth)
    bfs_node_list = sb_tree.bfs_node_list
    heapq.heapify(bfs_node_list)
    assert bfs_node_list == sb_tree.bfs_node_list
