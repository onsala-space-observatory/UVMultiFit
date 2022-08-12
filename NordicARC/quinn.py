import uvmultimodel as uvmod

##################################
# FRINGE FITTER:
#
def QuinnFF(IF, refant, doModel, doGlobal):
    """Perform a *Global Fringe Fitting* search in delay-rate space.

    It uses the Quinn estimator for the fringe peak.

    :Parameters:
    ------------
    **IF** : `int`
        The spectral window to fringe-fit (index referred to the data already read).
    **refant** : `int`
        The index of the antenna to use as reference. If this antenna is missing or has
        bad data, another refant will be picked automatically.
    **doModel** : `bool`
        If True, deconvolve the fixed model before the GFF.
    **doGlobal** : `bool`
        If True, globalize the solutions (i.e., use all available baselines to estimate
        the antenna gains).

    :Returns:
    ---------
    Returns a `list` of 5 elements:

    * An array with the delays (in seconds). One value per antenna (the reference
      antenna, as well as any other antenna with failed solutions, will have null
      values).
    * An array with the rates (in 1/s), one value per antenna.
    * An array with the phases (in radians), one value per antenna.
    * A double (the delay equivalent to one pixel in the delay/rate matrices).
    * A double (the rate equivalent to one pixel in the delay/rate matrices).
    """
    return uvmod.QuinnFF(IF, refant, doModel, doGlobal)