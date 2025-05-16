import itertools as it
from fractions import Fraction
import typing
import random

import cflib


def mediant(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (x[0] + y[0], x[1] + y[1])


class SternBrocot:
    ALPHABET = ["L", "R"]
    ROOT_LEFT_TUPLE, ROOT_MEDIANT_TUPLE, ROOT_RIGHT_TUPLE = (0, 1), (1, 1), (1, 0)
    ROOT_MOVE_LIST = []
    MOVE_BIT_DICT = {"L": "0", "R": "1"}

    @staticmethod
    def sb_left_child_parent_mediant(x: tuple[int, int], m: tuple[int, int]) -> tuple[int, int]:
        return (m[0] - x[0], m[1] - x[1])

    @staticmethod
    def sb_right_child_parent_mediant(y: tuple[int, int], m: tuple[int, int]) -> tuple[int, int]:
        return (m[0] - y[0], m[1] - y[1])

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
        self.cf = cflib.fraction_tuple_to_cf(self.mediant_tuple)
        self.move_string = "".join(self.move_list)
        L_bit, R_bit = type(self).MOVE_BIT_DICT["L"], type(self).MOVE_BIT_DICT["R"]
        self.binary_move_string = self.move_string.replace("L", L_bit).replace("R", R_bit)

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
        
    def __eq__(self, other: typing.Self) -> bool:
        return self.move_list == other.move_list
    
    def __lt__(self, other: typing.Self) -> bool:
        return len(self.move_list) < len(other.move_list) and other.move_list[0:len(self.move_list)] == self.move_list
    
    def __len__(self) -> int:
        return len(self.move_list)
    
    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.move_list}, {self.left_tuple}, {self.mediant_tuple}, {self.right_tuple})"
    
    def __str__(self) -> str:
        return f"{self.move_list}, {self.mediant_tuple}"


if __name__ == "__main__":
    move_list_length = 5

    random.seed(444)
    move_list = random.choices(SternBrocot.ALPHABET, k=move_list_length)

    N = SternBrocot.move_list_to_node(move_list)
    assert SternBrocot.move_list_to_cf(move_list) == N.cf

    N = SternBrocot.move_list_to_node(move_list)
    assert SternBrocot.move_list_to_fraction_tuple(move_list) == N.mediant_tuple

    numerator = 3
    denominator = 14
    fraction_tuple = Fraction(numerator, denominator).as_integer_ratio()

    N = SternBrocot.fraction_tuple_to_node(fraction_tuple)
    assert N.mediant_tuple == fraction_tuple

    N = SternBrocot.fraction_tuple_to_node(fraction_tuple)
    assert N.cf == cflib.fraction_tuple_to_cf(fraction_tuple)

    print(N.mediant_tuple, N.cf, len(N), N.binary_move_string)
    




