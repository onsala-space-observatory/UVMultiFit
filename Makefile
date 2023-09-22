CASADIR=/opt/casa-6.5.6-22-py3.8.el7
PYTHON=$(CASADIR)/bin/python3
LIBRARYPATH=$(CASADIR)/lib
PYTHONPATH=$(CASADIR)/lib/py/lib/python3.8/site-packages

install:
	$(PYTHON) setup.py install --prefix=$(CASADIR)

user:
	LD_LIBRARY_PATH=$(LIBRARYPATH) PYTHONPATH=$(PYTHONPATH) $(PYTHON) setup.py install --user

clean:
	rm -fvR build
