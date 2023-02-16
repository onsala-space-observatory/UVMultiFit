# UVMultiFit

A versatile library for fitting models directly to visibility data,
currently implemented in CASA. These models can depend on frequency
and fitting parameters in an arbitrary algebraic way. We have tested
the software with both synthetic data and real observations. In some
cases (e.g., sources with sizes smaller than the diffraction limit of
the interferometer), the results from the fit to the visibilities are
far superior to the output obtained from a mere analysis of the
deconvolved images. We give some illustrative examples in the software
documentation and in Marti-Vidal et al. (2014),
[A&A 563, 136, arXiv:1401.4984](http://arxiv.org/abs/1401.4984).
