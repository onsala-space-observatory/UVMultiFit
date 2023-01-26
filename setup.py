from setuptools import setup, Extension, find_packages

import platform

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
      include_dirs=[include_gsl_dir],
      long_description=long_description)
