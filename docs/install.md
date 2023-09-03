# Installation

Steps to install the `UVMultiFit` package for a **CASA 6.x** environment using **python 3.x**:

## Install dependencies

Our C++ extension makes use of the `fftw` and `gsl` libraries, so we
need to install two packages for the compilation to succeed:

#### For Linux (Ubuntu-like), in a terminal:

``` bash
$ sudo apt-get install libgsl-dev libfftw3-dev
```

#### For Mac OS, in a terminal:

If using Mac Ports:

``` bash
$ sudo port install gsl fftw-3
```

If using Homebrew:

``` bash
$ brew install gsl fftw
```

## Clone the github repository or update your installation

> NOTE THAT WE WANT THE **casa6** BRANCH!

Pick a directory where you want to install `UVMultiFit`, e.g.

``` bash
$ cd ~
$ mkdir -p ARC
$ cd ARC
$ git clone --branch casa6 https://github.com/onsala-space-observatory/UVMultiFit.git
$ cd UVMultiFit
```

If you already had a `git` based version installed, pull in the latest changes:

``` bash
$ git pull origin/casa6
```

## Build and install the python package

These instructions assume you have a CASA version 6.x installed,
either the normal release or with the ALMA and/or VLA pipeline
included.

The instructions below use an installation of `casa-6.4.0-16` (ALMA
pipeline) as an example.

You may want to install the package inside your CASA environment,
which would have the advantage that all users of that environment will
be able to use UVMultiFit. But it will require that you have write
permission to that environment.

In that case, **find out, where `casa` is installed** on your computer
and adjust the paths and version numbers below appropriately.  In the
example below this is in `/opt/casa-6.4.0-16`.

Then set the following two environment variables, adjusting
the paths and version numbers to your needs. For an installation into
your personal python environment, simply skip this step:

``` bash
$ PATH=/opt/casa-6.4.0-16/bin:$PATH
$ PYTHONPATH=/opt/casa-6.4.0-16/lib/py/lib/python3.8/site-packages
```

Now, in order to build the package, which will include compilation of
the C++ extension, run the following command:

``` bash
$ python3 -m pip install .
```

If successful, you may check that you can load the C++ module as follows:

``` bash
$ python3
Python 3.8.10 (default, Jun 25 2021, 20:16:58)
[GCC 5.3.1 20160406 (Red Hat 5.3.1-6)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import uvmultimodel as uvmod
>>> dir(uvmod)
['QuinnFF', '__doc__', '__file__', '__loader__', '__name__', '__package__', '__spec__',
 'clearPointers', 'modelcomp', 'setData', 'setModel', 'setNCPU', 'setNspw', 'setWork']
```

## Running the test suite

## Running your own model:

Any feedback and bug report should be sent either to the
[Nordic ARC Node](mailto:contact@nordic-alma.se) or to the source maintainer
[Michael Olberg](mailto:michael.olberg@chalmers.se).
