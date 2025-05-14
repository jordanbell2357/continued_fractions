import itertools as it
import cflib


class SternBrocot:
    ROOT_LEFT_TUPLE, ROOT_MIDDLE_TUPLE, ROOT_RIGHT_TUPLE = (0, 1), (1, 1), (1, 0)
    ROOT_MOVE_LIST = []
    ROOT_CF = [1]
    MOVE_BIT_DICT = {"L": "0", "R": "1"}

    @classmethod
    def move_list_to_cf(cls, move_list: list[str]) -> list[int]:
        run_tuple_list = [(k, len(list(g))) for k, g in it.groupby(move_list)]
        if not run_tuple_list:
            return cls.ROOT_CF
        if run_tuple_list[0][0] == 'L':
            run_tuple_list = [('R', 0)] + run_tuple_list
        cf = [run_tuple[1] for run_tuple in run_tuple_list]
        return cf        

    @classmethod
    def fraction_tuple_to_node(cls, fraction_tuple: tuple[int, int]):
        N = cls()
        while True:
            fraction_tuple_value = fraction_tuple[0] / fraction_tuple[1]
            node_middle_tuple_fraction_value = N.middle_tuple[0] / N.middle_tuple[1]
            if fraction_tuple_value < node_middle_tuple_fraction_value:
                N = N.L()
            elif fraction_tuple_value > node_middle_tuple_fraction_value:
                N = N.R()
            else:
                return N
            
    @classmethod
    def move_list_to_node(cls, move_list: list[str]):
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

    def __init__(self, move_list=None, left_tuple=None, middle_tuple=None, right_tuple=None):
        if move_list is None:
            self.move_list = type(self).ROOT_MOVE_LIST
            self.left_tuple = type(self).ROOT_LEFT_TUPLE
            self.middle_tuple = type(self).ROOT_MIDDLE_TUPLE
            self.right_tuple = type(self).ROOT_RIGHT_TUPLE
        else:
            self.move_list = move_list
            self.left_tuple = left_tuple
            self.middle_tuple = middle_tuple
            self.right_tuple = right_tuple
        self.move_string = "".join(self.move_list)
        self.binary_move_string = self.move_string.replace("L", type(self).MOVE_BIT_DICT["L"]) \
            .replace("R", type(self).MOVE_BIT_DICT["R"])
        self.cf = type(self).move_list_to_cf(self.move_list)

    def L(self):
        new_move_list = self.move_list + ['L']
        new_middle_tuple = cflib.mediant(self.left_tuple, self.middle_tuple)
        return type(self)(new_move_list, self.left_tuple, new_middle_tuple, self.middle_tuple)

    def R(self):
        new_move_list = self.move_list + ['R']
        new_middle_tuple = cflib.mediant(self.middle_tuple, self.right_tuple)
        return type(self)(new_move_list, self.middle_tuple, new_middle_tuple, self.right_tuple)
    
    def P(self):
        move_list = self.move_list
        if not move_list:
            return None
        previous_move = move_list[-1]
        new_move_list = move_list[0:-1]
        left_tuple, middle_tuple, right_tuple = self.left_tuple, self.middle_tuple, self.right_tuple
        if previous_move == "L":
            new_left_tuple = left_tuple
            new_middle_tuple = cflib.mediant_left_inverse(left_tuple, middle_tuple)
            new_right_tuple = middle_tuple
        elif previous_move == "R":
            new_left_tuple = middle_tuple
            new_middle_tuple = cflib.mediant_right_inverse(right_tuple, middle_tuple)
            new_right_tuple = right_tuple
        return type(self)(new_move_list, new_left_tuple, new_middle_tuple, new_right_tuple)                    
        
    def __eq__(self, other):
        return self.move_list == other.move_list
    
    def __lt__(self, other):
        return len(self.move_list) < len(other.move_list) and other.move_list[0:len(self.move_list)] == self.move_list
    
    def __repr__(self):
        return f"{type(self).__name__}({self.move_list}, {self.left_tuple}, {self.middle_tuple}, {self.right_tuple})"
    
    def __str__(self):
        return f"{self.move_list}, {self.middle_tuple}, {self.cf}"


if __name__ == "__main__":
    N = SternBrocot.move_list_to_node(["L", "L", "R"])
    print(N)
    N2 = N.L()
    print(N == N2)