# from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages

import numpy as np
import platform
import sys
import os
import fnmatch

####################################
# CONFIGURATION

# print("argv:", sys.argv)
# if len(sys.argv) > 2 and sys.argv[2].startswith("--prefix="):
#     CASA_INSTALLATION = sys.argv[2].split("=")[1]
# else:
#     CASA_INSTALLATION = os.getenv("CASA_INSTALLATION")
#
# print("using CASA_INSTALLATION", CASA_INSTALLATION)
#
# # typically you won't have to change anything below this line
# CASA_PYTHON_INCLUDE = CASA_INSTALLATION + "/include/python3.6m"
# CASA_PYTHON_SITE_PACKAGES = CASA_INSTALLATION + "/lib/python3.6/site-packages"
# CASA_NUMPY_INCLUDE = CASA_PYTHON_SITE_PACKAGES + "numpy/core/include"
# np_include = [ CASA_PYTHON_INCLUDE, CASA_NUMPY_INCLUDE ]
#
# if not os.path.exists(CASA_PYTHON_INCLUDE + "/Python.h"):
#     print("can't find 'Python.h' in your casa installation")
#     exit()
#
# if not os.path.exists(CASA_NUMPY_INCLUDE + "/numpy/arrayobject.h"):
#     print("can't find 'numpy/arrayobject.h' in your casa installation")
#     exit()

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
    _extra_link_args = ["-Xlinker", "-export-dynamic"]
elif platform.system() == "Darwin":
    _extra_link_args = ["-Xlinker", "-L/opt/local/lib"]
#####################################


#################
# SCRIPT STARTS #
#################
if BUILD_QUINN_FITTER:
    c_ext = Extension("uvmultimodel",
                      ["uvmultimodel.cpp", "QuinnFringe.cpp"],
                      define_macros=[('QUINN_FITTER', '0')],
                      libraries=['gsl', 'gslcblas', 'fftw3'],
                      extra_compile_args=["-Wno-deprecated", "-O3"],
                      include_dirs=[np.get_include()],
                      extra_link_args=_extra_link_args)
else:
    c_ext = Extension("uvmultimodel",
                      ["uvmultimodel.cpp"],
                      define_macros=[('QUINN_FITTER', '1')],
                      libraries=['gsl', 'gslcblas'],
                      extra_compile_args=["-Wno-deprecated", "-O3"],
                      include_dirs=[np.get_include()])

fh = open("README.md", "r")
long_description = fh.read()
fh.close()

setup(
    name="NordicARC",
    version="3.0.5",
    author="Ivan Marti-Vidal",
    author_email="i.marti-vidal@uv.es",
    description="Fit models directly to visibility data.",
    long_description=long_description,
    packages=find_packages(include=["NordicARC"]),
    url="https://github.com/onsala-space-observatory/UVMultiFit.git",
    ext_modules=[c_ext],
    include_dirs=[include_gsl_dir],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering',
    ],
)
