from distutils.core import setup, Extension
import numpy as np
import platform

####################################
# CONFIGURATION

######################
# SET THIS TO TRUE IF YOU WANT TO PLAY WITH FRINGE FITTING:
BUILD_QUINN_FITTER = True
######################

# DIRECTORY TO THE GSL LIBRARIES:
if platform.system() == "Linux":
  include_gsl_dir = "/usr/include"
elif platform.system() == "Darwin":
  include_gsl_dir = "/opt/local/include"
else:
  print("Unsupported platform: " + platform.system())
  exit()

# "-export-dynamic" doesn't work on Mac
if platform.system() == "Linux":
  _extra_link_args=["-Xlinker", "-export-dynamic"]
elif platform.system() == "Darwin":
  _extra_link_args=["-Xlinker", "-L/opt/local/lib"]
#####################################

#################
# SCRIPT STARTS #
#################
if BUILD_QUINN_FITTER:
  c_ext = Extension("_uvmultimodel",
                  ["_uvmultimodel.cpp", "_QuinnFringe.cpp"],
                  define_macros = [('QUINN_FITTER', '0')],
                  libraries=['gsl', 'gslcblas', 'fftw3'],
                  extra_compile_args=["-Wno-deprecated", "-O3"],
                  extra_link_args=_extra_link_args)
else:
  c_ext = Extension("_uvmultimodel",
                  ["_uvmultimodel.cpp"],
                  define_macros = [('QUINN_FITTER', '1')],
                  libraries=['gsl', 'gslcblas'],
                  extra_compile_args=["-Wno-deprecated", "-O3"])


setup(
    ext_modules=[c_ext],
    include_dirs=[include_gsl_dir] + [np.get_include()]
)
