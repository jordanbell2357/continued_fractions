from pprint import pprint

import farey
import number_theoretic_functions as ntf
import cw
import sb

if __name__ == "__main__":
    n=20

    f = farey.Farey(n)
    assert f.mertens_function == ntf.mertens(n)

    assert farey.mertens_function(n) == ntf.mertens(n)