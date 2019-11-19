import numpy as np
import scipy as sp
import pylab as pl
import os
from simutil import *

# TEST 8: DISPERSIVE FRINGE FITTING

# What to do:
# DoSimObs = False # Regenerate the synthetic data
# DoFit = True # Fit them
# casaexe='casa --nologger' # To run CASA in parallel!

# Fringe-Fit only these scans:
DOSCANs = [1]  # range(1,13)

# Array to use:
VLBAarr = """# observatory=VLBA\n# coordsys=XYZ
 -2112065.169063  3705356.505980  4726813.681152  25.  BR
 -1324009.296014  5332181.961648  3231962.393635  25.  FD
 -1995678.809977  5037317.705238  3357328.030390  25.  KP
 -1449752.549265  4975298.582200  3709123.846349  25.  LA
 -5464075.148965  2495248.295195  2148297.249841  25.  MK
  -130872.458404  4762317.103627  4226850.996905  25.  NL
 -2409150.357597  4478573.142538  3838617.335974  25.  OV
 -1640953.908700  5014816.029130  3575411.783269  25.  PT
  2607848.614131  5488069.569639  1932739.681584  25.  SC
"""

# Gains to solve for:
# Tau (ns)  Rate(mHz)   Phase (deg) Ion (TECU)
GAINS = {'BR':   [0.5,        6.,        20.,      1.],
         'FD':   [-0.45,     -3.,       -40.,      2.],
         'KP':   [0.65,      -1.,       170.,      0.5],
         'LA':   [-0.68,      4.,         5.,      1.5],
         'MK':   [0.17,       2.,       -20.,      3.],
         'NL':   [0.18,     -1.5,        10.,      0.25],
         'OV':   [-0.19,      1.,        -4.,      0.75],
         'PT':   [0.34,      1.2,        90.,      5.],
         'SC':   [-0.32,     -3.,       -90.,      6.]}

REFANT = 'LA'

# Name of MS:
imname = 'Compact_Source'

#########################
# SIMULATION PARAMETERS:

# On-source observing time:
TObs = 1.0  # in hr

# integration time (per visib.)
Tint = '2s'

# Duration of duty cycles:
Duty = 500.  # in seconds

# On-source time per duty cycle:
DutyON = 300.  # in seconds

# Source and array:
Coords = '10h00m00.0s 30d00m00.0s'
Hcenter = 0.0  # Hour angle at center of observations
Array = 'VLBA.array'
NU = '5.0GHz'
dNU = '0.512MHz'
Nchan = 128  # IF FREQ AND WIDTH
NUD = 1.44
dNUD = 0.512  # IF FREQ AND WIDTH IN GHz (DISPERSIVE CASE)

TSYS = 60.  # in K.
tau0 = 0.01  # Not important
seed = 42  # Random seed
ADD_NOISE = True

#################
# SCRIPT STARTS #
#################

TECFAC = 1.e16*40.3/3.e8
mas2rad = np.pi/180./3600./1000.
r2d = 180./np.pi
twopi = 2.*np.pi

# Write the array file:
arf = open('VLBA.array', 'w')
print >> arf, VLBAarr
arf.close()

# Compute gains with their right units:
for ant in GAINS.keys():
    GAINS[ant][0] *= 1.e-9
    GAINS[ant][1] *= 1.e-3  # /qa.convertfreq(NU)['value']
    GAINS[ant][2] *= np.pi/180.
    GAINS[ant][3] *= TECFAC

if DoSimObs:
    print 'Generating %s' % imname

    print ' SETTING ANTENNA PARAMETERS'
    util = simutil('')
    stnx, stny, stnz, stnd, padnames, arrname, arrpos = util.readantenna(Array)
    eta_p, eta_s, eta_b, eta_t, eta_q, t_rx = util.noisetemp(telescope='VLBA', freq=NU)
    eta_a = eta_p * eta_s * eta_b * eta_t

    print '\n\n RESETTING THE TSYS TO %.1f K.\n\n' % TSYS
    t_rx = TSYS
    t_sky = 260.
    t_ground = 260.

    # Prepare template measurement set:
    print '\n\n PREPARING MS TEMPLATE\n\n'

    VLBA = me.observatory('VLBA')
    mount = 'alt-az'
    refdate = '2017/01/01/00:00:00'
    usehourangle = True

    os.system('rm -rf %s.ms' % imname)
    sm.open('%s.ms' % imname)

    sm.setconfig(telescopename='VLBA', x=stnx, y=stny, z=stnz,
                 dishdiameter=stnd.tolist(),
                 mount=mount, antname=padnames, padname=padnames,
                 coordsystem='global', referencelocation=VLBA)

    sm.setspwindow(spwname='spw0', freq=NU, deltafreq=dNU,
                   freqresolution=dNU,
                   nchannels=Nchan, refcode="BARY",
                   stokes='RR LL')

    sm.setfield(sourcename="FFTarget", sourcedirection=Coords,
                calcode="TARGET", distance='0m')

    mereftime = me.epoch('TAI', refdate)

    print ' Will shift the date of observation to match the Hour Angle range\n'

    sm.settimes(integrationtime=Tint, usehourangle=usehourangle,
                referencetime=mereftime)

    # Figure out scan distribution:
    NSCAN = int(TObs*3600./DutyON)
    Ttot = NSCAN*Duty
    print '\n\n  Total duration of experiment: %.2f h (%i scans)\n\n' % (Ttot/3600., NSCAN)

    T0s = [Hcenter*3600. + Ttot*(i-NSCAN/2.)/NSCAN for i in range(NSCAN)]

    starttimes = []
    stoptimes = []
    sources = []
    for i in range(NSCAN):
        sttime = T0s[i]
        endtime = (sttime + Duty)
        starttimes.append(str(sttime)+'s')
        stoptimes.append(str(endtime)+'s')
        sources.append("FFTarget")

    for n in range(NSCAN):
        sm.observemany(sourcenames=[sources[n]], spwname='spw0', starttimes=[
                       starttimes[n]], stoptimes=[stoptimes[n]], project='UVFIT_TEST')
    sm.close()

    # Apply the gains to the data:
    print '\n\n WRITING VISIBILITIES \n\n'

    tb.open(imname+'.ms/ANTENNA')
    ANTNAMES = tb.getcol('NAME')
    tb.close()

    tb.open(imname+'.ms/SPECTRAL_WINDOW')
    FREQS = tb.getcol('CHAN_FREQ')[:, 0]
    FREQS -= FREQS[0]
    tb.close()

    tb.open(imname+'.ms', nomodify=False)
    A1 = tb.getcol('ANTENNA1')
    A2 = tb.getcol('ANTENNA2')
    T = tb.getcol('TIME')
    T -= T[0]

    GARR = np.zeros((len(ANTNAMES), 4))
    for ant in range(len(ANTNAMES)):
        GARR[ant, :] = GAINS[ANTNAMES[ant]]

    for i in range(len(T)):
        if not i % 1000:
            sys.stdout.write('\r  Writing vis #%i of %i' % (i, len(T)))
            sys.stdout.flush()
        DELAY = GARR[A1[i]][0] - GARR[A2[i]][0]
        RATE = GARR[A1[i]][1] - GARR[A2[i]][1]
        PHASE = GARR[A1[i]][2] - GARR[A2[i]][2]
        TOTPHASE = 2.*np.pi*(DELAY*FREQS + RATE*T[i]) + PHASE
        DATA = tb.getcell('DATA', i)
        DATA[0, :] = np.exp(1.j*TOTPHASE)
        DATA[1, :] = np.exp(-1.j*TOTPHASE)
        tb.putcell('DATA', i, DATA)
    tb.close()

    # Add noise:
    print '\n\n Corrupting data...\n'

    os.system('rm -rf %s.noisy.ms' % imname)
    os.system('cp -r %s.ms %s.noisy.ms' % (imname, imname))

    sm.openfromms('%s.noisy.ms' % imname)
    sm.setdata(fieldid=0, spwid=0)
    sm.setseed(seed)
    sm.setnoise(spillefficiency=eta_s, correfficiency=eta_q,
                antefficiency=eta_a, trx=t_rx,
                tau=tau0, tatmos=t_sky, tground=t_ground, tcmb=2.725,
                mode="tsys-manual", senscoeff=-1)
    sm.corrupt()
    sm.done()

    clearcal('%s.noisy.ms' % imname, addmodel=True)

    # DISPERSIVE CASE:
    print '\n\n NOW FOR THE DISPERSIVE CASE!\n\n'

    os.system('rm -rf %s.dispersive.noisy.ms' % imname)
    os.system('cp -r %s.noisy.ms %s.dispersive.noisy.ms' % (imname, imname))

    tb.open(imname+'.dispersive.noisy.ms/SPECTRAL_WINDOW', nomodify=False)
    FREQS = tb.getcol('CHAN_FREQ')
    FREQS[:, 0] = 1.e9*np.linspace(NUD, dNUD+NUD, Nchan)
    tb.putcol('CHAN_FREQ', FREQS)
    NUR = tb.getcol('REF_FREQUENCY')
    NUR[:] = NUD*1.e9
    tb.putcol('REF_FREQUENCY', NUR)
    CHW = tb.getcol('CHAN_WIDTH')
    CHW[:] = dNUD*1.e9/Nchan
    tb.putcol('CHAN_WIDTH', CHW)
    CHW = tb.getcol('EFFECTIVE_BW')
    CHW[:] = dNUD*1.e9
    tb.putcol('EFFECTIVE_BW', CHW)
    CHW = tb.getcol('RESOLUTION')
    CHW[:] = dNUD*1.e9/Nchan
    tb.putcol('RESOLUTION', CHW)
    CHW = tb.getcol('TOTAL_BANDWIDTH')
    CHW[:] = dNUD*1.e9
    tb.putcol('TOTAL_BANDWIDTH', CHW)
    tb.close()

    tb.open(imname+'.dispersive.noisy.ms', nomodify=False)
    for i in range(len(T)):
        if not i % 1000:
            sys.stdout.write('\r  Writing vis #%i of %i' % (i, len(T)))
            sys.stdout.flush()
        DELAY = GARR[A1[i]][0] - GARR[A2[i]][0]
        RATE = GARR[A1[i]][1] - GARR[A2[i]][1]
        PHASE = GARR[A1[i]][2] - GARR[A2[i]][2]
        DISPDEL = (GARR[A1[i]][3] - GARR[A2[i]][3])/FREQS[:, 0]**2.
        TOTPHASE = 2.*np.pi*((DISPDEL + DELAY)*FREQS[:, 0] + RATE*T[i]) + PHASE
        DATA = tb.getcell('DATA', i)
        DATA[0, :] = np.exp(1.j*TOTPHASE)
        DATA[1, :] = np.exp(-1.j*TOTPHASE)
        tb.putcell('DATA', i, DATA)
    tb.close()

    # ADD NOISE TO THE DISPERSIVE CASE:
    print '\n\n Corrupting data...\n'

    sm.openfromms('%s.dispersive.noisy.ms' % imname)
    sm.setdata(fieldid=0, spwid=0)
    sm.setseed(seed)
    sm.setnoise(spillefficiency=eta_s, correfficiency=eta_q,
                antefficiency=eta_a, trx=t_rx,
                tau=tau0, tatmos=t_sky, tground=t_ground, tcmb=2.725,
                mode="tsys-manual", senscoeff=-1)
    sm.corrupt()
    sm.done()

if DoFit:
    print("---------------------------------------------")
    print("TEST 8")
    print("---------------------------------------------")
    tempfile = open('STEP8_FIT.py', 'w')
    print '\n\n\n   NON-DISPERSIVE FIT \n\n\n'
    import time

    tempfile.close()

    # Names of the antennas:
    ANTS = GAINS.keys()
    tb.open(imname+'.noisy.ms/ANTENNA')
    ANTNAMES = tb.getcol('NAME')
    tb.close()

    # Gain model for UVMultiFit (only phase gains are solved):
    allgains = {}
    ai = 0
    ri = -1
    for ant in range(len(ANTNAMES)):
        pc = 3*ai
        if ANTNAMES[ant] != REFANT:
            allgains[ant] = '2.*3.1416*(p[%i]*(nu - nu0)*(1.e-9) + p[%i]*t) + p[%i]' % (pc, pc+1, pc+2)
            ai += 1
        else:
            ri = ant

    # Number of parameters (not counting the gains of REFANT!):
    Npar = len(allgains.keys())*3
    message = 'GLOBAL FRINGE FITTING:\n'
    fitTime = 0
    for DOSCAN in DOSCANs:
        tic = time.time()
        message += '\n##############\n   SCAN   %i\n##############\n\n' % DOSCAN
        myfit = uvm.uvmultifit(vis='%s.dispersive.noisy.ms' % imname, spw='0', column='data', scans=[DOSCAN],
                               stokes='RR', model=['delta'], var=['0,0,1'], write='model',
                               p_ini=[0.0 for i in range(Npar)], OneFitPerChannel=False, finetune=True,
                               phase_gains=allgains, LMtune=[1.e-3, 10., 1.e-4, 5, 1.e-3])

        # Fringe fitting based on Quinn estimator:
        QGains = myfit.QuinnFF(0, ri, 0, 1)

        nu = qa.convertfreq(NU)['value']
        pini = []
        bounds = []
        Dtau = 3.*QGains[3]*1.e9
        DRate = 3.*QGains[4]*1.e3
        # Print Quinn estimates:
        if QGains != -1:
            message += '\n\n From Quinn Fringing (though phases are not globalized!): \nAntenna | Quantity     |  True       |  Quinn Est. |   Difference  \n\n'
            for ant in range(len(ANTS)):
                an = ANTNAMES[ant]
                # True gains (i.e., those used in the simulations):
                ph = (GAINS[an][2] - GAINS[REFANT][2])*r2d
                de = (GAINS[an][0] - GAINS[REFANT][0])*1.e9
                rt = (GAINS[an][1] - GAINS[REFANT][1])*1.e3

                # Fitted gains:
                fph = QGains[2][ant]*r2d
                fde = QGains[0][ant]*1.e9
                frt = QGains[1][ant]*1.e3

                # Set fitted gains as a-prioris for the least-squares:
                if an != REFANT:
                    pini += [QGains[0][ant]*1.e9,
                             QGains[1][ant], QGains[2][ant]]
                    bounds += [[pini[-3]-Dtau, pini[-3]+Dtau],
                               [pini[-2]-DRate, pini[-2]+DRate], [-twopi, twopi]]

                # Determine phase ambiguities (w.r.t. true phase gains used in the simul):
                dph = ph-fph
                if dph > 180.:
                    dph -= 360.
                if dph < -180.:
                    dph += 360.

                # Print fit-true gain differences:
                message += '  %s    |  Delay (ns)  |  % 9.3f  |  % 9.3f  |  % 9.3f \n' % (an, de, fde, de-fde)
                message += '  %s    |  Rate (mHz)  |  % 9.3f  |  % 9.3f  |  % 9.3f \n' % (an, rt, frt, rt-frt)
                message += '  %s    |  Phase (deg) |  % 9.3f  |  % 9.3f  |  % 9.3f \n\n' % (an, ph, fph, dph)

        # If the Quinn fitter is not compiled, just set a-prioris manually:
        else:
            message = '\n\n Quinn-based Fringe Fitting is not compiled!\n Will set a-prioris to the TRUE values!'
            for ant in range(len(ANTS)):
                an = ANTNAMES[ant]
                ph = (GAINS[an][2] - GAINS[REFANT][2])
                de = (GAINS[an][0] - GAINS[REFANT][0])
                rt = (GAINS[an][1] - GAINS[REFANT][1])
                pini += [de*1.e9, rt, fph]
        print message

        # Set a-prioris and parameter bounds:
        myfit.p_ini = pini
        myfit.bounds = bounds

        # Fit!!!!
        myfit.checkInputs()
        myfit.fit()

        # Write model visibilities into the MS:
        myfit.writeModel()

        # Read the gains (they are arranged as defined in the "allgains" models):
        FittedGains = []
        pf = 0
        for pi in range(len(ANTNAMES)):
            FittedGains.append([])
            if ANTNAMES[pi] != REFANT:
                FittedGains[-1] = myfit.result['Parameters'][3*pf:3*pf+3]
                pf += 1
            else:
                FittedGains[-1] = [0., 0., 0.]  # REFANT has null gains.

        # Print comparison between true gains and best-fit gains:
        message += '\n\n Least-Squares GFF: \nAntenna | Quantity     |  True       |  Global FF. |   Difference  \n\n'
        for ant in range(len(ANTS)):
            an = ANTNAMES[ant]

            # True gains:
            ph = (GAINS[an][2] - GAINS[REFANT][2])*r2d
            de = (GAINS[an][0] - GAINS[REFANT][0])*1.e9
            rt = (GAINS[an][1] - GAINS[REFANT][1])*1.e3

            # Best-fit gains:
            fph = FittedGains[ant][2]*r2d
            fde = FittedGains[ant][0]
            frt = FittedGains[ant][1]*1.e3

            # Solve phase ambiguity w.r.t. true phase gains:
            dph = ph-fph
            if dph > 180.:
                dph -= 360.
                fph += 360.
            if dph < -180.:
                dph += 360.
                fph -= 360.

            # Print differences:
            message += '  %s    |  Delay (ns)  |  % 9.3f  |  % 9.3f  |  % 9.3f \n' % (an, de, fde, de-fde)
            message += '  %s    |  Rate (mHz)  |  % 9.3f  |  % 9.3f  |  % 9.3f \n' % (an, rt, frt, rt-frt)
            message += '  %s    |  Phase (deg) |  % 9.3f  |  % 9.3f  |  % 9.3f \n\n' % (an, ph, fph, dph)

        tac = time.time()
        fitTime += tac-tic

        message += '\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n' % (tac-tic)
        print message

        # Divide data by model visibilities:
        print '\n NOW, CALIBRATING SCAN'

        ms.open('%s.dispersive.noisy.ms' % imname, nomodify=False)
        ms.selectinit(datadescid=0)
        ms.select({'scan_number': DOSCAN})

        ALLDATS = ms.getdata(['data', 'corrected_data', 'model_data'])
        ALLDATS['corrected_data'][:] = ALLDATS['data']/ALLDATS['model_data']
        ms.putdata(ALLDATS)
        ms.close()

        # Plot the calibrated (i.e., data/model) visibilities:
        print '\n PLOTTING...'

        tb.open('%s.dispersive.noisy.ms' % imname)
        SCAN = tb.getcol('SCAN_NUMBER') == DOSCAN
        ANT1 = tb.getcol('ANTENNA1')[SCAN]
        ANT2 = tb.getcol('ANTENNA2')[SCAN]
        DATA = tb.getcol('DATA')[0, :, SCAN]
        CORRDATA = tb.getcol('CORRECTED_DATA')[0, :, SCAN]

        fig = pl.figure()
        fig.subplots_adjust(wspace=0.01, hspace=0.01)
        NPLOT = len(ANTNAMES)-1
        k = 0
        for i in range(len(ANTNAMES)):
            if ri != i:
                BASI = np.logical_or((ANT1 == ri)*(ANT2 == i), (ANT2 == ri)*(ANT1 == i))
                AverData = np.average(DATA[BASI, :], axis=0)
                sub = fig.add_subplot(2, NPLOT, k+1)
                sub.plot(np.angle(AverData)*r2d, '.r')
                pl.setp(sub.get_xticklabels(), 'visible', False)
                if k > 0:
                    pl.setp(sub.get_yticklabels(), 'visible', False)
                else:
                    sub.set_ylabel('Phase (deg.)')
                pl.ylim((-180., 180.))
                sub.set_title(ANTNAMES[i])
                sub = fig.add_subplot(2, NPLOT, NPLOT + k+1)
                pl.setp(sub.get_xticklabels(), 'visible', False)
                if k > 0:
                    pl.setp(sub.get_yticklabels(), 'visible', False)
                else:
                    sub.set_ylabel('Amp (Norm.)')
                pl.ylim((0., 1.5))
                sub.plot(np.abs(AverData), '-b')
                k += 1
        pl.savefig('TEST8.ORIGINAL_FRINGES_%i.png' % DOSCAN)

        fig = pl.figure()
        fig.subplots_adjust(wspace=0.01, hspace=0.01)
        NPLOT = len(ANTNAMES)-1
        k = 0
        for i in range(len(ANTNAMES)):
            if ri != i:
                BASI = np.logical_or((ANT1 == ri)*(ANT2 == i), (ANT2 == ri)*(ANT1 == i))
                AverData = np.average(CORRDATA[BASI, :], axis=0)
                sub = fig.add_subplot(2, NPLOT, k+1)
                sub.plot(np.angle(AverData)*r2d, '.r')
                pl.setp(sub.get_xticklabels(), 'visible', False)
                if k > 0:
                    pl.setp(sub.get_yticklabels(), 'visible', False)
                else:
                    sub.set_ylabel('Phase (deg.)')
                pl.ylim((-180., 180.))
                sub.set_title(ANTNAMES[i])
                sub = fig.add_subplot(2, NPLOT, NPLOT + k+1)
                pl.setp(sub.get_xticklabels(), 'visible', False)
                if k > 0:
                    pl.setp(sub.get_yticklabels(), 'visible', False)
                else:
                    sub.set_ylabel('Amp (Norm.)')
                pl.ylim((0., 1.5))
                sub.plot(np.abs(AverData), '-b')
                k += 1

        pl.savefig('TEST8.NON_DISPERSIVE_GFF.SCAN-%i.png' % DOSCAN)

        # Release memory
        #  del myfit

        # Write results:
        resf = open('test8.dat', 'w')
        print >> resf, '\n\n\nTEST 8: GLOBAL FRINGE FITTING\n'
        print >> resf, message

        #  print >> resf, '  TOTAL FRINGING TIME:  %.2f SECONDS '%fitTime
        #  resf.close()

        # NOW FIT WITH THE DISPERSIVE-DELAY MODEL!!!!!
        print '\n\n\n   DISPERSIVE CASE \n\n\n'

        #  import time

        # Read antenna names and frequency info (already done, but doesn't harm):

        #  ANTS = GAINS.keys()
        #  tb.open(imname+'.dispersive.noisy.ms/ANTENNA')
        #  ANTNAMES = tb.getcol('NAME')
        #  tb.close()

        # Read frequency info:
        tb.open(imname+'.dispersive.noisy.ms/SPECTRAL_WINDOW')
        FREQS = tb.getcol('CHAN_FREQ')
        tb.close()

        # Square of average frequency:
        NuAvgSq = np.average(FREQS)**2.

        # FIRST, PERFORM AN ORDINARY GFF, TO ESTIMATE THE EFFECTIVE GROUP DELAY:
        # This fit is *exactly* the same as the one above.

        #  allgains = {}
        #  ai = 0
        #  ri = -1
        #  for ant in range(len(ANTNAMES)):
        #    pc = 3*ai
        #    if ANTNAMES[ant] != REFANT:
        #      allgains[ant] = '2.*3.1416*(p[%i]*(nu - nu0)*(1.e-9) + p[%i]*t) + p[%i]'%(pc,pc+1,pc+2)
        #      ai += 1
        #    else:
        #      ri = ant
        #  Npar = len(allgains.keys())*3

        #  print 'GLOBAL FRINGE FITTING (NON-DISPERSIVE MODEL)\n'

        #  fitTime = 0

        #  for DOSCAN in DOSCANs:

        #    tic = time.time()

        #    message += '\n##############\n   SCAN   %i\n##############\n\n'%DOSCAN

        #    myfit = uvm.uvmultifit(vis='%s.dispersive.noisy.ms'%imname, spw='0', column='data',scans = [DOSCAN],
        #                 stokes = 'RR', model = ['delta'], var=['0,0,1'], write='model',
        #                 p_ini=[0.0 for i in range(Npar)], OneFitPerChannel=False, finetune=True,
        #                 phase_gains = allgains,LMtune=[1.e-3,10.,1.e-4,5,1.e-3])

        #    QGains = myfit.QuinnFF(0,ri,0,1)
        #    nu = qa.convertfreq(NU)['value']
        #    r2d = 180./np.pi
        #    twopi = 2.*np.pi
        #    pini = []
        #    bounds = []
        #    Dtau = QGains[3]*1.e9 ; DRate = QGains[4]*1.e3

        #    if QGains != -1:

        #      for ant in range(len(ANTS)):
        #        an = ANTNAMES[ant]

        #        fph = QGains[2][ant]*r2d
        #        fde = QGains[0][ant]*1.e9
        #        frt = QGains[1][ant]*1.e3

        #        if an != REFANT:
        #          pini += [QGains[0][ant]*1.e9, QGains[1][ant], QGains[2][ant]]
        #          bounds += [[pini[-3]-Dtau,pini[-3]+Dtau],[pini[-2]-DRate,pini[-2]+DRate],[-twopi,twopi]]

        #    else:
        #
        #      for ant in range(len(ANTS)):
        #        an = ANTNAMES[ant]
        #        ph = (GAINS[an][2] - GAINS[REFANT][2])
        #        de = (GAINS[an][0] - GAINS[REFANT][0])
        #        rt = (GAINS[an][1] - GAINS[REFANT][1])
        #        pini += [de*1.e9, rt, fph]

        #    myfit.p_ini = pini
        #    myfit.bounds = bounds
        #    myfit.checkInputs()
        #    myfit.fit()

        #    FittedGains = []
        #    pf = 0
        #    for pi in range(len(ANTNAMES)):
        #      FittedGains.append([])
        #
        #      if ANTNAMES[pi] != REFANT:
        #        FittedGains[-1] = myfit.result['Parameters'][3*pf:3*pf+3]
        #        pf += 1
        #      else:
        #        FittedGains[-1] = [0.,0.,0.]

        # NOW, WE FIX THE "EFFECTIVE DELAY" AND SOLVE FOR THE DISPERSIVE ONE:
        print 'GLOBAL FRINGE FITTING (DISPERSIVE MODEL!!!)\n'
        allgains = {}
        ai = 0
        ri = -1
        for ant in range(len(ANTNAMES)):
            pc = 3*ai
            if ANTNAMES[ant] != REFANT:
                allgains[ant] = '2.*3.1416*(p[%i]*(nu - nu0)*(1.e-9) + %.16e*t  + %.16e*(p[%i]/%.16e*(nu-nu0)  + p[%i]/nu - p[%i]/nu0)) + p[%i]' % \
                    (pc+1, FittedGains[ant][1], TECFAC, pc, NuAvgSq, pc, pc, pc+2)
                pini[pc] = 0.0
                pini[pc+1] = FittedGains[ant][0]
                bounds[pc] = [None, None]
                bounds[pc+1] = [FittedGains[ant][0] -
                                2*Dtau, FittedGains[ant][0]+2*Dtau]
                ai += 1
            else:
                ri = ant

        # three parameters, since we FIX THE RATE
        Npar = len(allgains.keys())*3

        myfit.phase_gains = allgains
        myfit.p_ini = pini
        myfit.bounds = bounds
        myfit.checkInputs()
        myfit.initModel()
        myfit.fit()
        myfit.writeModel()

        # Read the disperive-model gains back:
        FittedGainsD = []
        pf = 0
        for pi in range(len(ANTNAMES)):
            FittedGainsD.append([])
            if ANTNAMES[pi] != REFANT:
                FittedGainsD[-1] = myfit.result['Parameters'][3*pf:3*pf+3]
                pf += 1
            else:
                FittedGainsD[-1] = [0., 0., 0.]

        # Print the true-fitted comparison:
        message += '\n\n Least-Squares GFF: \nAntenna | Quantity     |  True       |  Global FF. |   Difference  \n\n'
        for ant in range(len(ANTS)):
            an = ANTNAMES[ant]

            # True gains:
            ph = (GAINS[an][2] - GAINS[REFANT][2])*r2d
            de = (GAINS[an][0] - GAINS[REFANT][0])*1.e9
            rt = (GAINS[an][1] - GAINS[REFANT][1])*1.e3
            disp = (GAINS[an][3] - GAINS[REFANT][3])/TECFAC

            # Fitted gains:
            # + effect of Delay and TEC at 1.4GHz!!
            fph = FittedGainsD[ant][2]*r2d
            # i.e., effective + dispersive at NuAvg.
            fde = FittedGainsD[ant][1] + TECFAC * \
                FittedGainsD[ant][0]/NuAvgSq*1.e9
            frt = FittedGains[ant][1]*1.e3
            # This is the TEC:
            fdp = FittedGainsD[ant][0]

            # Amiguities of the phase gains:
            dph = ph-fph
            if dph > 180.:
                dph -= 360.
                fph += 360.
            if dph < -180.:
                dph += 360.
                fph -= 360.

            message += '  %s    |  Delay (ns)  |  % 9.3f  |  % 9.3f  |  % 9.3f \n' % (an, de, fde, de-fde)
            message += '  %s    |  Rate (mHz)  |  % 9.3f  |  % 9.3f  |  % 9.3f \n' % (an, rt, frt, rt-frt)

            # For a meaningful comparison of the phase gains, we still have to correct the phase by the effect of the
            # (nu-nu0) and TEC/nu0**2. at 1.4GHz!!!
            #      message += '  %s    |  Phase (deg) |  % 9.3f  |  % 9.3f  |  % 9.3f \n'%(an,ph,fph,dph)
            message += '  %s    |  TEC (TECU)  |  % 9.3f  |  % 9.3f  |  % 9.3f \n\n' % (an, -disp, -fdp, -disp+fdp)

        tac = time.time()
        fitTime += tac-tic

        message += '\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n' % (tac-tic)
        print message

        # Divide data by model (to get the calibrated visibilities):
        print '\n NOW, CALIBRATING SCAN'

        ms.open('%s.dispersive.noisy.ms' % imname, nomodify=False)
        ms.selectinit(datadescid=0)
        ms.select({'scan_number': DOSCAN})

        ALLDATS = ms.getdata(['data', 'corrected_data', 'model_data'])
        ALLDATS['corrected_data'][:] = ALLDATS['data']/ALLDATS['model_data']
        ms.putdata(ALLDATS)
        ms.close()

        print '\n PLOTTING...'
        tb.open('%s.dispersive.noisy.ms' % imname)
        SCAN = tb.getcol('SCAN_NUMBER') == DOSCAN
        ANT1 = tb.getcol('ANTENNA1')[SCAN]
        ANT2 = tb.getcol('ANTENNA2')[SCAN]
        DATA = tb.getcol('DATA')[0, :, SCAN]
        CORRDATA = tb.getcol('CORRECTED_DATA')[0, :, SCAN]

        fig = pl.figure()
        fig.subplots_adjust(wspace=0.01, hspace=0.01)
        NPLOT = len(ANTNAMES)-1
        k = 0
        for i in range(len(ANTNAMES)):
            if ri != i:
                BASI = np.logical_or((ANT1 == ri)*(ANT2 == i), (ANT2 == ri)*(ANT1 == i))
                AverData = np.average(CORRDATA[BASI, :], axis=0)
                sub = fig.add_subplot(2, NPLOT, k+1)
                sub.plot(np.angle(AverData)*r2d, '.r')
                pl.setp(sub.get_xticklabels(), 'visible', False)
                if k > 0:
                    pl.setp(sub.get_yticklabels(), 'visible', False)
                else:
                    sub.set_ylabel('Phase (deg.)')
                pl.ylim((-180., 180.))
                sub.set_title(ANTNAMES[i])
                sub = fig.add_subplot(2, NPLOT, NPLOT + k+1)
                pl.setp(sub.get_xticklabels(), 'visible', False)
                if k > 0:
                    pl.setp(sub.get_yticklabels(), 'visible', False)
                else:
                    sub.set_ylabel('Amp (Norm.)')
                pl.ylim((0., 1.5))
                sub.plot(np.abs(AverData), '-b')
                k += 1
        pl.savefig('TEST8.DISPERSIVE_GFF.SCAN-%i.png' % DOSCAN)

        # Release memory
        del myfit

        # Write results:
        #  resf = open('test8_disp.dat','w')
        print >> resf, '\n\n\nTEST 8: GLOBAL FRINGE FITTING\n'
        print >> resf, message

    print >> resf, '  TOTAL FRINGING TIME:  %.2f SECONDS ' % fitTime
    resf.close()
