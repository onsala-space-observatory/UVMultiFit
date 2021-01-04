import os
import glob
import time

DoSimObs = True
DoFit = True

# CASA EXECUTION COMMAND (MAY BE 'casapy --nologger',
# DEPENDING ON YOUR INSTALATION)
casaexe = "casa --nologger"

tic = time.time()

test = "TEST1-2_PROFILES_ELLIPTICITY"

currdir = os.getcwd()
print(currdir)

print("ENTERING %s" % (test))
os.chdir(test)
os.system("%s -c RUNTEST1.py" % casaexe)

tac = time.time()
print("\n THE TEST LASTED %.1f SECONDS.\n\n" % (tac-tic))
