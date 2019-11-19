import os
import glob
import time

# REDO SIMULATIONS? (SET TO TRUE THE FIRST TIME YOU
# RUN THE TEST SUIT)
DoSimObs = True

# DO THE FITS? (OF COURSE!)
DoFit = True

# CASA EXECUTION COMMAND (MAY BE 'casapy --nologger',
# DEPENDING ON YOUR INSTALATION)
casaexe = "casa --nologger"

#################
# SCRIPT STARTS #
#################

tic = time.time()

RESDIR = "ALLTESTS.RESULTS"
SCRDIR = "FIT_SCRIPTS"

alltest = [f[:-2] for f in os.popen("ls -d TEST*/")]

os.system("rm -rf %s %s" % (RESDIR, SCRDIR))
os.system("mkdir %s" % (RESDIR))
os.system("mkdir %s" % (SCRDIR))

currdir = os.getcwd()
for test in alltest:
    print("ENTERING %s" % (test))
    os.chdir(test)
    scr = open(glob.glob("MASTER*.py")[0])
    scr2 = open("aux.py", 'w')
    print >> scr2, "DoSimObs = %s" % str(DoSimObs)
    print >> scr2, "DoFit = %s" % str(DoFit)
    print >> scr2, "casaexe = '%s'" % str(casaexe)
    for line in scr.readlines():
        print >> scr2, line[:-1]
    scr2.close()
    scr.close()
    os.system("%s -c aux.py" % casaexe)
    os.system("mv test*.dat ../%s" % (RESDIR))
    os.system("mv *.png ../%s" % (RESDIR))
    os.system("mv STEP*.py ../%s" % (SCRDIR))

    os.system("rm *.log *.last aux.py")
    os.chdir(currdir)

tac = time.time()
off = open("dur.dat", 'w')
print >> off, "\n\n THE FULL TEST LASTED %.1f SECONDS.\n\n" % (tac-tic)
off.close()

os.system("cat %s/test*.dat dur.dat > %s/FIT_RESULTS.dat" % (RESDIR, RESDIR))
os.system("rm %s/test*.dat dur.dat" % (RESDIR))
