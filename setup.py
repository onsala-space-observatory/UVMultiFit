# from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages

import platform
import sys
import os

####################################
def locate_osx_lib(library, install_names):
    """Locate library headers on OSX.
    Looks for "library" in common locations, where "library" can be a file or directory
    returns prefix to library.
    *install_names is a tuple containing names to pass to brew, ports, and conda respectively"""
# Homebrew: /usr/local
    prefix_homebrew = "/usr/local"
# Macports: /opt/local
    prefix_macports = "/opt/local"
# Conda: CONDA_PREFIX
    prefix_conda = os.getenv("CONDA_PREFIX")

    if os.path.exists(prefix_homebrew + "/include/" + library):
        return prefix_homebrew
    elif os.path.exists(prefix_macports + "/include/" + library):
        return prefix_macports
    elif prefix_conda and os.path.exists(prefix_conda + "/include/" + library):
        return prefix_conda
    else:
        print("\nError: can't find " + library + " on your system.")
        print("To install with homebrew:")
        print("    brew install " + install_names[0])
        print("To install with Mac Ports:")
        print("    sudo port install " + install_names[1])
        print("To install with Conda:")
        print("    conda install " + install_names[2])
        exit()

####################################
# CONFIGURATION

print("argv:", sys.argv)
if len(sys.argv) > 2 and sys.argv[2].startswith("--prefix="):
    CASA_INSTALLATION = sys.argv[2].split("=")[1]
else:
    CASA_INSTALLATION = os.getenv("CASA_INSTALLATION")

# if platform.system() == "Darwin":
#     CASA_INSTALLATION = CASA_INSTALLATION + "/Python.framework/Versions/2.7"

print("using CASA_INSTALLATION", CASA_INSTALLATION)

# typically you won't have to change anything below this line
CASA_PYTHON_INCLUDE = CASA_INSTALLATION + "/include/python2.7"
CASA_PYTHON_SITE_PACKAGES = CASA_INSTALLATION + "/lib/python2.7/site-packages"
CASA_NUMPY_INCLUDE = CASA_PYTHON_SITE_PACKAGES + "/numpy/core/include"
np_include = [CASA_PYTHON_INCLUDE, CASA_NUMPY_INCLUDE]

if not os.path.exists(CASA_PYTHON_INCLUDE + "/Python.h"):
    print("can't find 'Python.h' in your casa installation")
    exit()

if not os.path.exists(CASA_NUMPY_INCLUDE + "/numpy/arrayobject.h"):
    print("can't find 'numpy/arrayobject.h' in your casa installation")
    exit()

######################
# SET THIS TO TRUE IF YOU WANT TO PLAY WITH FRINGE FITTING:
BUILD_QUINN_FITTER = True

# DIRECTORY TO THE GSL & FFTW3 LIBRARIES:
if platform.system() == "Linux":
    include_gsl_dir = "/usr/include"
    include_fftw3_dir = include_gsl_dir
    if not os.path.exists(include_gsl_dir + "/gsl"):
        print("can't find GSL headers on your system (expected to find " + include_gsl_dir + "/gsl/)")
        exit()
    if not os.path.exists(include_gsl_dir + "/fftw3.h"):
        print("can't find FFTW3 headers on your system (expected to find " + include_gsl_dir + "/fftw3.h)")
        exit()
elif platform.system() == "Darwin":
#Try to locate GSL
        prefix = locate_osx_lib("gsl", ["gsl", "gsl", "gsl"])
        include_gsl_dir = prefix + "/include"
        lib_gsl_dir = prefix + "/lib"
#Try to locate FFTW3
        prefix = locate_osx_lib("fftw3.h", ["fftw", "fftw-3", "fftw"])
        include_fftw3_dir = prefix + "/include"
        lib_fftw3_dir = prefix + "/lib"
else:
    print("Unsupported platform: " + platform.system())
    exit()

# "-export-dynamic" doesn't work on Mac
if platform.system() == "Linux":
    _extra_link_args = ["-Xlinker", "-export-dynamic"]
    _extra_compile_args=["-Wno-deprecated", "-O3"]
elif platform.system() == "Darwin":
    _extra_link_args=["-Xlinker", "-L" + lib_gsl_dir, "-Wl,-rpath," + lib_gsl_dir, "-L" + lib_fftw3_dir, "-Wl,-rpath," + lib_fftw3_dir]
    _extra_compile_args=["-Wno-deprecated", "-O3", "-stdlib=libc++", "-isysroot/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"]
#####################################

#################
# SCRIPT STARTS #
#################
if BUILD_QUINN_FITTER:
    c_ext = Extension("_uvmultimodel",
                      ["_uvmultimodel.cpp", "_QuinnFringe.cpp"],
                      define_macros=[('QUINN_FITTER', '0')],
                      libraries=['gsl', 'gslcblas', 'fftw3'],
                      extra_compile_args=_extra_compile_args,
                      extra_link_args=_extra_link_args)
else:
    c_ext = Extension("_uvmultimodel",
                      ["_uvmultimodel.cpp"],
                      define_macros=[('QUINN_FITTER', '1')],
                      libraries=['gsl', 'gslcblas'],
                      extra_compile_args=_extra_compile_args)

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
    include_dirs=[include_gsl_dir, include_fftw3_dir] + np_include,
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering',
    ],
)
