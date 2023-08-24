from setuptools import setup, Extension, find_packages

import sys
import os
import platform

# SET THIS TO TRUE IF YOU WANT TO PLAY WITH FRINGE FITTING:
BUILD_QUINN_FITTER = True

if len(sys.argv) > 2 and sys.argv[2].startswith("--prefix="):
    CASA_INSTALLATION = sys.argv[2].split("=")[1]
else:
    CASA_INSTALLATION = os.getenv("CASA_INSTALLATION")

print("using CASA_INSTALLATION", CASA_INSTALLATION)

# typically you won't have to change anything below this line
CASA_PYTHON_INCLUDE = CASA_INSTALLATION + "/include/python3.6m"
CASA_PYTHON_SITE_PACKAGES = CASA_INSTALLATION + "/lib/python3.6/site-packages"
CASA_NUMPY_INCLUDE = CASA_PYTHON_SITE_PACKAGES + "/numpy/core/include"
np_include = [CASA_PYTHON_INCLUDE, CASA_NUMPY_INCLUDE]

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
    _extra_link_args = ["-Xlinker", "-export-dynamic"]
elif platform.system() == "Darwin":
    _extra_link_args = ["-Xlinker", "-L/opt/local/lib"]

if BUILD_QUINN_FITTER:
    c_ext = Extension("uvmultimodel",
                      ["uvmultimodel.cpp", "QuinnFringe.cpp"],
                      define_macros=[('QUINN_FITTER', '0')],
                      libraries=['gsl', 'gslcblas', 'fftw3'],
                      extra_compile_args=["-Wno-deprecated", "-O3"])
else:
    c_ext = Extension("uvmultimodel",
                      ["uvmultimodel.cpp"],
                      define_macros=[('QUINN_FITTER', '1')],
                      libraries=['gsl', 'gslcblas'],
                      extra_compile_args=["-Wno-deprecated", "-O3"])

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(packages=find_packages(include=["NordicARC"]),
      ext_modules=[c_ext],
      include_dirs=[include_gsl_dir] + np_include,
      long_description=long_description,)
