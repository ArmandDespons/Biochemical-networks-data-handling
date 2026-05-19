import numpy as np

from numpy.linalg import matrix_rank


def value_to_key(dic: dict, val):

    if not isinstance(dic, dict): raise ValueError(f"First argument must be a dictionary instance")

    if not val in dic.values(): raise ValueError(f"Value '{val}' is not found in the dictionary values")

    return list(dic.keys())[list(dic.values()).index(val)]



def _matrix_parser(matrix: np.ndarray, rows: dict = None):
    """
        Parse a matrix by removing its linearly dependent rows
    """

    output = []

    if rows:
        new_row = {}

    rk = 0
    for i in range(1, matrix.shape[0]+1):
        
        if matrix_rank(matrix[:i, :])==(rk+1):
            output.append( matrix[i-1, :] )

            if rows:
                new_row[value_to_key(rows, i-1)] = rk          

            rk += 1

    if rows:
        return np.array(output), new_row
    
    else: 
        return np.array(output)


def _check_coeffs(coeffs: dict[str, float], reactants, products): 

    for s in coeffs.keys(): 

        if s in [r.formula for r in reactants]:
            coeffs.update({s: -np.abs(coeffs[s])})

        elif s in [p.formula for p in products]:
            coeffs.update({s: np.abs(coeffs[s])})

        else: 
            ValueError(f"Species {s} is not found in the reactants & products")

    


