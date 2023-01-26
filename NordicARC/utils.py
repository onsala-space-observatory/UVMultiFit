import ast
import re

import logging

from typing import List, Any, Union

valid_casa_pos = re.compile(r"J2000 [0-9]{1,2}h[0-9]{2}m[0-9]{2}.[0-9]s -?[0-9]{1,2}d[0-9]{2}m[0-9]{2}.[0-9]s")

def is_valid_python(param: str) -> bool:
    """Check if passed string is valied Python code.

    Args:
        param (str): The string of code to check.

    Returns:
        bool: True if string is valid code.
    """
    logging.debug(f"is_valid_python('{param}')")
    try:
        ast.parse(param)
    except SyntaxError:
        return False
    return True

def is_casa_position(param: str) -> bool:
    """Check if string defines a valid casa position.

    Args:
        param (str): The string to check.

    Returns:
        bool: True if string is a valid casa position.
    """
    logging.debug(f"is_casa_position('{param}')")

    if not isinstance(param, str):
        return False
    return valid_casa_pos.fullmatch(param) is not None

def is_valid_stokes(param: str) -> bool:
    """Check if string is one of the allowed Stokes parameters.

    Args:
        param (str): The string to check.

    Returns:
        bool: True if string is a valid Stokes parameter.
    """
    logging.debug(f"is_valid_stokes({param})")
    allowed = ['PI', 'I', 'Q', 'U', 'V',
               'XX', 'YY', 'XY', 'YX',
               'RR', 'LL', 'LR', 'RL']
    return isinstance(param, str) and param in allowed

def is_list_of_floats(param: Any, length=2) -> bool:
    """Check if passed argument is a list of floats.

    Args:
        param (list): A list of floats.
        length (int, optional): The length of the required list. Defaults to 2.

    Returns:
        bool: True if all elements in the list are floats and the length is correct.
    """
    logging.debug(f"is_list_of_floats({param}) ({length})")
    return isinstance(param, list) and len(param) == length and all(isinstance(elem, float) for elem in param)

def is_list_of_int(param: Any) -> bool:
    """Check if passed argument is a list of integers.

    Args:
        param (list): A list of integers.

    Returns:
        bool: True if the passed argument is a list of integers or a list of lists of integers.
    """
    logging.debug(f"is_list_of_int({param}")
    if not isinstance(param, list):
        return False

    if all(isinstance(elem, int) for elem in param):
        return True

    for p in param:
        if not isinstance(p, list) or not all(isinstance(elem, int) for elem in p):
            return False
    return True

def get_list_of_strings(param: Any) -> List[str]:
    """Turn string argument into list of strings, if necessary.

    Args:
        param (str): A string or list of strings.

    Returns:
        list: A list of strings or None.
    """
    logging.debug(f"get_list_of_strings({param})")
    if isinstance(param, str):
        if len(param) == 0:
            return []
        return [param]

    if isinstance(param, list):
        if len(param) == 0:
            logging.warning("param is list of length zero")
            return []
        for i, p in enumerate(param):
            if not isinstance(p, str):
                logging.warning(f"element {i} is not a string")
                return []
        return param
    return []

def check_proper_motion(pm: Any, model: List[str]) -> List[List[float]]:
    result = [[0.0, 0.0] for i in model]
    if pm != 0.0:
        if len(pm) != len(model):
            logging.error("'proper_motion' must be a list of the same length as model")
        else:
            result = []
            for i in range(len(model)):
                if not is_list_of_floats(pm[i]):
                    logging.error("each element of 'proper_motion' must be a list of two floats")
                else:
                    result.append(pm[i])
    return result


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                        format='%(levelname)s - %(message)s')
    pos = "J2000 12h34m56.0s -01d02m03.0s"
    print(is_casa_position(pos))
    proper_motion = check_proper_motion(0.0, ['delta'])
    print(proper_motion)
    proper_motion = check_proper_motion([[1.0, 2.0], [3.0, 4.0]], ['delta', 'disc'])
    print(proper_motion)
