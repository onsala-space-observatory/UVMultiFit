import ast
import re

import logging

valid_casa_pos = re.compile(r"J2000 [0-9]{1,2}h[0-9]{2}m[0-9]{2}.[0-9]s -?[0-9]{1,2}d[0-9]{2}m[0-9]{2}.[0-9]s")

def is_valid_python(param):
    logging.debug(f"is_valid_python('{param}')")
    try:
        ast.parse(param)
    except SyntaxError:
        return False
    return True

def is_casa_position(param):
    logging.debug(f"is_casa_position('{param}')")

    if not isinstance(param, str):
        return False
    return valid_casa_pos.fullmatch(param) is not None

def is_valid_stokes(param):
    logging.debug(f"is_valid_stokes({param})")
    allowed = ['PI', 'I', 'Q', 'U', 'V',
               'XX', 'YY', 'XY', 'YX',
               'RR', 'LL', 'LR', 'RL']
    return isinstance(param, str) and param in allowed

def is_list_of_floats(param, length=2):
    logging.debug(f"is_list_of_floats({param}) ({length})")
    return isinstance(param, list) and len(param) == length and all(isinstance(elem, float) for elem in param)

def is_list_of_int(param):
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
    logging.debug(f"get_list_of_strings({param})")
    if type(param) == str:
        if len(param) == 0:
            return None
        return [param]

    if type(param) == list:
        if len(param) == 0:
            logging.warning("param is list of length zero")
            return None
        for i, p in enumerate(param):
            if type(p) != str:
                logging.warning(f"element {i} is not a string")
                return None
        return param
    return None

def channeler(spw, width=1, maxchans=[3840, 3840, 3840, 3840]):
    """ Function to convert a string with spw selection into lists of channels to select/average."""
    logging.debug(f"channeler({spw})")

    if spw == '':
        spw = ','.join(list(map(str, range(len(maxchans)))))

    entries = spw.split(',')
    output = [[] for i in maxchans]

    for entry in entries:
        check = entry.split(':')
        if check[0] == '*':
            check[0] = f"0~{len(maxchans)-1}"

        spws = list(map(int, check[0].split('~')))
        if len(spws) == 1:
            selspw = [spws[0]]
        else:
            selspw = list(range(spws[0], spws[1] + 1))

        for sp in selspw:
            if sp + 1 > len(maxchans):
                errstr = f"there are only {len(maxchans)} spw in the data, please revise the 'spw' parameter"
                logging.error(errstr)
                return [False, errstr]

            if len(check) == 1:
                channel_ranges = [f"0~{maxchans[sp] - 1}"]
            else:
                chans = check[1]
                channel_ranges = chans.split(';')
            logging.debug(f"channel_range {sp}: {channel_ranges}")
            ranges = []
            for channel_range in channel_ranges:
                ch1, ch2 = list(map(int, channel_range.split('~')))
                if ch1 > ch2:
                    errstr = f"{ch1} is larger than {ch2}, revise channels for spw {sp}"
                    logging.error(errstr)
                    return [False, errstr]
                ch2 = min([ch2, maxchans[sp] - 1])
                for i in range(ch1, ch2 + 1, width):
                    ranges.append(range(i, min([(i + width), ch2 + 1])))

            output[sp] = ranges
    return [True, output]

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                        format='%(levelname)s - %(message)s')
    pos = "J2000 12h34m56.0s -01d02m03.0s"
    print(is_casa_position(pos))

    print(channeler('')[0])
    print(channeler('0')[0])
    print(channeler('1')[0])
    print(channeler('0~3')[0])
    print(channeler('0,2,3')[0])
    # print(channeler('1,3,5')[0])
    print(channeler('*:10~50')[0])
    print(channeler('1~3:30~40')[0])
