from distutils.core import setup, Extension
import numpy as np
import platform
import os
import fnmatch

####################################
# CONFIGURATION

# DEFINE YOUR CASA INSTALLATION
CASA_INSTALLATION = "/home/olberg/Python/casa-pipeline-release-5.6.1-8.el7"

# typically you won't have to change anything below this line
CASA_PYTHON_INCLUDE = CASA_INSTALLATION + "/include/python2.7"
CASA_PYTHON_SITE_PACKAGES = CASA_INSTALLATION + "/lib/python2.7/site-packages"
CASA_NUMPY_INCLUDE = CASA_PYTHON_SITE_PACKAGES + "/numpy/core/include"
np_include = [ CASA_PYTHON_INCLUDE, CASA_NUMPY_INCLUDE ]

if not os.path.exists(CASA_PYTHON_INCLUDE + "/Python.h"):
    print("can't find 'Python.h' in your casa installation")
    exit()

if not os.path.exists(CASA_NUMPY_INCLUDE + "/numpy/arrayobject.h"):
    print("can't find 'numpy/arrayobject.h' in your casa installation")
    exit()

######################
# SET THIS TO TRUE IF YOU WANT TO PLAY WITH FRINGE FITTING:
BUILD_QUINN_FITTER = True

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

# for root, dirnames, filenames in os.walk(CASA_INSTALLATION):
#     for dirname in fnmatch.filter(dirnames, "numpy"):
#         np_include = root.replace("/numpy", "")
#         break

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
    include_dirs=[include_gsl_dir] + np_include
)
