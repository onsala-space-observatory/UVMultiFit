# Unit tests

Run the unit tests with the casatools from casa-6.4.0-16. Flag `-rP`
captures print statements and includes their output in the test
output.

``` shell
PYTHONPATH=/home/data/casa-6.4.0-16/lib/py/lib/python3.8/site-packages pytest-3 -rP
```

If modifications to C++ module are needed:

``` shell
edit uvmultimodel.cpp
make user    # make sure Makefile refers to same casa version as above.
```
