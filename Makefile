CASADIR=/media/olberg/d1c658f8-95c1-49a9-ab3e-b77c48f5dbe3/data/casa-6.1.2-7-pipeline-2020.1.0.36
PYTHON=$(CASADIR)/bin/python3

install:
	export PYTHONPATH=$(CASADIR)/lib/py/lib/python3.6/site-packages
	$(PYTHON) setup.py install     #  --prefix=$(CASADIR)

build:
	export PYTHONPATH=$(CASADIR)/lib/py/lib/python3.6/site-packages
	$(PYTHON) setup.py build

clean:
	rm -fvR build

#  	rm -fvR dist
#  	find . -name '*~' | xargs rm -fv
#  	find . -name 'aux.py' | xargs rm -fv
#  	find . -name 'dur.dat' | xargs rm -fv
#  	find . -name 'test*.dat' | xargs rm -fv
#  	find . -name 'STEP*.py' | xargs rm -fv
#  	find . -name '*.png' | xargs rm -fv
#  	find . -name '*.log' | xargs rm -fv
#  	find . -name '*.last' | xargs rm -fv
#  	rm -fvR test/TEST1-2_PROFILES_ELLIPTICITY/TEST1.*
#  	rm -fvR test/TEST1-2_PROFILES_ELLIPTICITY/Disc
#  	rm -fvR test/TEST1-2_PROFILES_ELLIPTICITY/Disc.model
#  	rm -fvR test/TEST1-2_PROFILES_ELLIPTICITY/modelfit.dat
#  	rm -fvR test/TEST3_MULTICOMP/TEST3.*
#  	rm -fvR test/TEST3_MULTICOMP/Discs_Gaussian
#  	rm -fvR test/TEST3_MULTICOMP/Discs_Gaussian.model
#  	rm -fvR test/TEST3_MULTICOMP/modelfit.dat
#  	rm -fvR test/TEST4_WIDEFIELD_MOSAIC/WideField
#  	rm -fvR test/TEST4_WIDEFIELD_MOSAIC/WideField.model
