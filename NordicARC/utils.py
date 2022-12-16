import ast
import re

import logging

valid_casa_pos = re.compile(r"J2000 [0-9]{1,2}h[0-9]{2}m[0-9]{2}.[0-9]s -?[0-9]{1,2}d[0-9]{2}m[0-9]{2}.[0-9]s")

def is_valid_python(param):
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

def is_casa_position(param):
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

def is_valid_stokes(param):
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

def is_list_of_floats(param, length=2):
    """Check if passed argument is a list of floats.

    Args:
        param (list): A list of floats.
        length (int, optional): The length of the required list. Defaults to 2.

    Returns:
        bool: True if all elements in the list are floats and the length is correct.
    """
    logging.debug(f"is_list_of_floats({param}) ({length})")
    return isinstance(param, list) and len(param) == length and all(isinstance(elem, float) for elem in param)

def is_list_of_int(param):
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

def get_list_of_strings(param):
    """Turn string argument into list of strings, if necessary.

    Args:
        param (str): A string or list of strings.

    Returns:
        list: A list of strings or None.
    """
    logging.debug(f"get_list_of_strings({param})")
    if isinstance(param, str):
        if len(param) == 0:
            return None
        return [param]

    if isinstance(param, list):
        if len(param) == 0:
            logging.warning("param is list of length zero")
            return None
        for i, p in enumerate(param):
            if not isinstance(p, str):
                logging.warning(f"element {i} is not a string")
                return None
        return param
    return None

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                        format='%(levelname)s - %(message)s')
    pos = "J2000 12h34m56.0s -01d02m03.0s"
    print(is_casa_position(pos))
