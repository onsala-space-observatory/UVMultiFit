CASADIR=/home/michael/Python/casa-6.1.2-7-pipeline-2020.1.0.36/lib/py
PYTHON=$(CASADIR)/bin/python3

install:
	$(PYTHON) setup.py install --prefix=$(CASADIR)

user:
	$(PYTHON) setup.py install --user

#clean:
#	rm -fvR $(CASADIR)/lib/python2.7/site-packages/NordicARC*.egg-info
