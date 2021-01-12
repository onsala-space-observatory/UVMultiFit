CASADIR=/Applications/CASA.app
#/home/olberg/Python/casa-pipeline-release-5.6.1-8.el7

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CASADIR := $(CASADIR)/Contents/Frameworks/Python.framework/Versions/2.7
endif
PYTHON=$(CASADIR)/bin/python

install:
	$(PYTHON) setup.py install --prefix=$(CASADIR)

user:
	export CASA_INSTALLATION=$(CASADIR); $(PYTHON) setup.py install --user

clean:
	rm -fvR $(CASADIR)/lib/python2.7/site-packages/NordicARC*.egg-info
