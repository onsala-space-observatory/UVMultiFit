import os
import glob
import time

DoSimObs = True
DoFit = True

# CASA EXECUTION COMMAND (MAY BE 'casapy --nologger',
# DEPENDING ON YOUR INSTALATION)
casaexe = "casa --nologger"

tic = time.time()

test = "TEST4_WIDEFIELD_MOSAIC"

currdir = os.getcwd()
print(currdir)

print("ENTERING %s" % (test))
os.chdir(test)
src = open(glob.glob("MASTER*.py")[0])
aux = open("aux.py", 'w')
print("DoSimObs = %s" % str(DoSimObs), file=aux)
print("DoFit = %s" % str(DoFit), file=aux)
# print("casaexe = '%s'" % str(casaexe), file=aux)
for line in src.readlines():
    print(line[:-1], file=aux)
aux.close()
src.close()
os.system("%s -c aux.py" % casaexe)
# os.system("rm *.log *.last aux.py")

tac = time.time()
print("\n THE TEST LASTED %.1f SECONDS.\n\n" % (tac-tic))
os.chdir(currdir)
