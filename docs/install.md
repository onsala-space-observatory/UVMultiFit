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

> NOTE THAT WE WANT THE **develop** BRANCH!

Pick a directory where you want to install `UVMultiFit`, e.g.

``` bash
$ cd ~
$ mkdir -p ARC
$ cd ARC
$ git clone --branch develop https://github.com/onsala-space-observatory/UVMultiFit.git
$ cd UVMultiFit
```

If you already had a `git` based version installed, pull in the latest changes:

``` bash
$ git pull origin/develop
```

## Build and install the python package

These instructions assume you have a CASA version 6.x installed,
either the normal release or with the ALMA and/or VLA pipeline
included.

The instructions below use an installation of `casa-6.5.6-22` as an example.

You may want to install the package inside your CASA environment,
which would have the advantage that all users of that environment will
be able to use UVMultiFit. But it will require that you have write
permission to that environment.

In that case, **find out, where `casa` is installed** on your computer
and adjust the paths and version numbers below appropriately.  In the
example below this is in `/opt/casa-6.5.6-22-py3.8.el7`.

Open up the Makefile in an editor and make sure the first line matches **your** casa installation


```
CASADIR=/opt/casa-6.5.6-22-py3.8.el7
```

Now, in order to build the package, which will include compilation of
the C++ extension, run either one of the following commands:

``` bash
$ make install     # this would install into the python3 environment of your casa installation
$ make user        # this will install into your local, personal python environment
```

Also edit the script `activate` and change the line

```
VIRTUAL_ENV="/opt/casa-6.5.6-22-py3.8.el7"
```

to fit **your** casa installation. Then, source this file.

``` bash
source ./activate
```

which should change your bash prompt to `(casa6)`. You may now test
that you can successfully import the compiled C++ extension:

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

Move into the `test-cases` subdirectory and perform any of the
tests. You must however start with `TEST1_disc.py`, the other tests
can then be run in any order. Note, when you run these tests for the
first time they will first generate simulated model data, and
therefore take some time. If the simulation data is already present in
that subdirectory, only the fit will be performed, which is much
faster. Each test will generate an output file `<test-name>.out`,
which can be compared with a file of that same name in subdirectory
`outputs`.


``` bash
$ cd test_cases
$ python3 TEST1_disc.py
$ diff TEST1_disc.out ../outputs/TEST1_disc.out
$ ...
```


## Running your own model:

Any feedback and bug report should be sent either to the
[Nordic ARC Node](mailto:contact@nordic-alma.se) or to the source maintainer
[Michael Olberg](mailto:michael.olberg@chalmers.se).
