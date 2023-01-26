from typing import List

import uvmultimodel as uvmod                      # type: ignore

def QuinnFF(IF: int, refant: int, doModel: bool, doGlobal: bool) -> List:
    """Perform a *global fring fitting* search in delay-rate space.

    Args:
        IF (int): The spectral window to fringe-fit (index refers to the data already read).
        refant (int): The index of the antenna to use as a reference, if this antenna is missing
            or has bad data, another refant will be picked automatically.
        doModel (bool): If True, deconvolve the fixed model before the global finge fitting.
        doGlobal (bool): If True, globalize the solutions, i.e. use all available baselines
            to estimate the antenna gains

    Returns:
        list: a list of 5 elements:
            * An array with the delays (in seconds), one value per antenna.
                (the reference antenna, as well as any other antenna with
                failed solutions, will have null values).
            * An array with the rates (in 1/s), one value per antenna.
            * An array with the phases (in radians), one value per antenna.
            * A double (the delay equivalent to one pixel in the delay/rate matrices).
            * A double (the rate equivalent to one pixel in the delay/rate matrices).
    """

    return uvmod.QuinnFF(IF, refant, doModel, doGlobal)
