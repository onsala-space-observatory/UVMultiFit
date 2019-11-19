import numpy as np
import scipy as sp
import pylab as pl
import os
from simutil import *

# What to do:
# DoSimObs = True
# DoFit = True
# casaexe = 'casa --nologger'

# Array to use:
arrayconfig = 'pdbi-a.cfg'

# Gain corruption factor:
Gfac = 1.5

# Center Freq:
Nu = '150.0GHz'

# Image size (pixels)
Npix = 1000

# Pixel size:
cell = '0.01arcsec'
cellN = 0.01

# Number of frequency channels:
Nnu = 100

# Channel width:
ChW = '0.01GHz'

# Total integration time:
Ttot = '3600s'

# For test 5:
imname = 'Gaussian'
# type   flux    major       minor       PA  freq   spec.ind.    Direction
s = ['gaussian', 3.5, '1.5arcsec', '1.0arcsec',
     '45.0deg', Nu, 0.0, "J2000 10h00m00.0s +10d00m02.0s"]

config = '.'.join(arrayconfig.split('.')[:-1])

vis = '%s/%s.%s.ms' % (imname, imname, arrayconfig.split('.')[0])

if DoSimObs:
    # TEST 5: ANTENNA GAIN FIT FOR PdB SIMULATION, COUPLED TO SOURCE PARAMETERS

    # FOR THE LATEST CASA VERSION (5.4):
    # SIMOBSERVE IS COMPLETELY BROKEN IF WE WANT TO SIMULATE ARRAYS OTHER THAN NRAO's ONES.
    # WE HAVE TO SWTICH TO THE SM TOOL.

    os.system('rm -rf %s*' % imname)
    cl.done()

    # for si in s:
    cl.addcomponent(dir=s[7], flux=s[1], fluxunit='Jy', freq=s[5], shape=s[0],
                    majoraxis=s[2], minoraxis=s[3], positionangle=s[4],
                    spectrumtype='spectral index', index=s[6])

    ia.fromshape(imname+'.model', [Npix, Npix, 1, Nnu], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad = qa.convert(qa.quantity(cell), "rad")['value']
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    rdi, ddi = s[7].split()[1:]
    cs.setreferencevalue([qa.convert(rdi, 'rad')['value'],
                          qa.convert(ddi, 'rad')['value']], type="direction")
    cs.setreferencevalue(s[5], 'spectral')
    cs.setreferencepixel(0, 'spectral')

    cs.setincrement(ChW, 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    ia.close()

    os.system('mkdir %s' % imname)
    sm.open(vis)

    # Setting the observatory and the observation:
    util = simutil('')
    antlist = os.getenv("CASAPATH").split(' ')[0] + "/data/alma/simmos/" + arrayconfig
    stnx, stny, stnz, stnd, antnames, arrname, arrpos = util.readantenna(antlist)

    # CASA dictatorship imposes limits to observatory names:
    ARR = me.observatory('VLBA')

    H0 = 0.0
    refdate = '2017/01/01/00:00:00'
    usehourangle = True
    mount = 'alt-az'
    integ = '2s'

    sm.setconfig(telescopename='VLBA', x=stnx, y=stny, z=stnz,
                 dishdiameter=stnd.tolist(),
                 mount=mount, antname=antnames, padname=antnames,
                 coordsystem='global', referencelocation=ARR)

    sm.setspwindow(spwname='spsim', freq=Nu,
                   deltafreq=ChW,
                   freqresolution=ChW,
                   nchannels=Nnu, refcode="BARY",
                   stokes='RR')

    sm.setfield(sourcename='SIM5', sourcedirection=s[7],
                calcode="TARGET", distance='0m')

    mereftime = me.epoch('TAI', refdate)

    sm.settimes(integrationtime=integ, usehourangle=True,
                referencetime=mereftime)
    sm.observemany(sourcenames=['SIM5'], spwname='spsim', starttimes=[
                   '0s'], stoptimes=[Ttot], project='polsimulate')
    sm.close()

    clearcal(vis=vis, addmodel=True)
    ft(vis=vis, usescratch=True, spw='0', field='SIM5', model=imname+'.model')

    ms.open(vis, nomodify=False)
    ms.selectinit(datadescid=0)
    CDAT = ms.getdata(['data', 'model_data'])
    CDAT['data'][:] = CDAT['model_data']
    ms.putdata(CDAT)
    ms.close()

    os.system('cp -rf %s %s.noisy' % (vis, vis))

    sm.openfromms('%s.noisy' % vis)
    sm.setdata(fieldid=[0], spwid=0)
    sm.setnoise(trx=50.,
                tau=0.05, tatmos=260., tground=290., tcmb=2.725,
                mode="tsys-manual", senscoeff=-1)
    sm.corrupt()
    sm.done()

    sm.close()

    # GAIN CORRUPTION:
    # The gain of antenna #2 will increase from 1 to Gfac during the
    # first half of the experiment. Then, it will be kept constant = Gfac

    # Get time range of observations:
    ms.open('%s.noisy' % vis)
    trange = ms.range(['time'])['time']/86400.
    ms.close()

    # time at half of the experiment:
    tmed = np.average(trange)

    # Change gain of antenna 2:
    ms.open('%s.noisy' % vis, nomodify=False)

    ms.selectinit(datadescid=0)
    data = ms.getdata(['data', 'antenna1', 'antenna2', 'time'])
    def GainFac(t):
        Gfunc = np.ones(np.shape(t))
        mask = t < tmed
        Gfunc[mask] = 1.0 + (t[mask]-trange[0])/(tmed-trange[0])*(Gfac - 1.)
        Gfunc[np.logical_not(mask)] = Gfac
        return Gfunc

    Gfunc = GainFac(data['time']/86400.)

    mask = data['antenna2'] == 2
    data['data'][:, :, mask] *= Gfunc[mask]
    mask = data['antenna1'] == 2
    data['data'][:, :, mask] *= Gfunc[mask]
    ms.putdata(data)
    ms.close()

if DoFit:
    print("---------------------------------------------")
    print("TEST 5")
    print("---------------------------------------------")
    # Get time range of observations:
    ms.open('%s.noisy' % vis)
    trange = ms.range(['time'])['time']/86400.
    ms.close()

    tmed = np.average(trange)

    tempfile = open('STEP5_FIT.py', 'w')

    # Source parameters:
    Gvar = '0.,0.,p[0],p[1],p[2],p[3]'

    ma = float(s[2].split('a')[0])
    mi = float(s[3].split('a')[0])
    pa = float(s[4].split('d')[0])

    # Last 2 params are those of the gain model function
    # (we will use a piece-wise linear function).
    # Both gains are set below their true value
    # (the solution is unstable if one is too high and
    # the other is too low).
    Pars = [s[1], ma, mi/ma, pa, 1., Gfac]
    BiasFac = [0.7, 1.3, 0.8, 0.75, 0.8, 0.8]

    Trues = Pars

    string = 'LMtune = [1.e-5,10.,1.e-5,10]'
    print >> tempfile, string

    string = 'S = [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]' % tuple(Pars)
    print >> tempfile, string

    string = "visname = '%s'" % ('%s.noisy' % vis)
    print >> tempfile, string

    # Define the gain function for antenna #2:
    gainFunc = 'pieceWise(p[4],p[5],%.8e,%.8e)' % (0., (tmed-trange[0])*86400.)
    string = "amp_gains = {2:'%s'}" % gainFunc
    print >> tempfile, string

    if s[0] == 'disk':
        s[0] = 'disc'
    if s[0] == 'gaussian':
        s[0] = 'Gaussian'

    string = "modelshape = ['%s']" % s[0]
    print >> tempfile, string

    NuF = float(Nu.split('G')[0])*1.e9
    string = "modvars = ['%s']" % Gvar
    print >> tempfile, string

    string = 'pini = [%.6e, %.6e, %.6e, %.6e, %.6e, %.6e]' % tuple(
        [Pars[pi]*BiasFac[pi] for pi in range(len(Pars))])
    print >> tempfile, string

    string = 'parbound = [[0.,None],[0.,None],[0.,None],[0.,90.],[0.75,1.25],[1.1,1.9]]'
    # ll = []
    # for pi in range(len(Pars)):
    #   ll += [Pars[pi]*(1. - BiasFac[pi]),Pars[pi]*(1. + BiasFac[pi])]
    # string = 'parbound = [[%.6e,%.6e],[%.6e,%.6e],[%.6e,%.6e],[%.6e,%.6e],[%.6e,%.6e],[%.6e,%.6e]]'%tuple(ll)
    print >> tempfile, string

    string = "STK = 'RR'"
    print >> tempfile, string

    lines = open('test5.py')
    for l in lines.readlines():
        print >> tempfile, l[:-1]
    lines.close()
    tempfile.close()

    pl.ioff()
    Cfile = 'TEST5.CLEAN'
    Rfile = 'TEST5.RESIDUALS'
    clearcal('%s.noisy' % vis)

    os.system('rm -rf %s.*' % Cfile)
    os.system('rm -rf %s.*' % Rfile)

    clean('%s.noisy' % vis,
          imagename=Cfile, cell=cell,
          imsize=Npix, niter=0)

    ia.open('%s.image' % Cfile)
    resdat = ia.getchunk()[:, :, 0, 0]
    ia.close()
    fig = pl.figure()
    sub = fig.add_subplot(111)
    IMsize = cellN*Npix/2.
    ims = sub.imshow(np.transpose(resdat), origin='lower',
                     extent=[IMsize, -IMsize, -IMsize, IMsize])
    sub.set_xlabel('RA offset (mas)')
    sub.set_ylabel('Dec offset (mas)')
    cb = pl.colorbar(ims)
    cb.set_label('Jy/beam')
    pl.savefig('%s.png' % Cfile)

    impeak = np.max(resdat)

    os.system('%s -c STEP5_FIT.py' % casaexe)

    os.system('rm -rf %s.*' % Rfile)
    clean('%s.noisy' % vis,
          imagename=Rfile, cell=cell,
          imsize=Npix, niter=0)

    impeak *= 0.01

    ia.open('%s.image' % Rfile)
    resdat = ia.getchunk()[:, :, 0, 0]
    ia.close()
    fig = pl.figure()
    sub = fig.add_subplot(111)
    IMsize = cellN*Npix/2.
    ims = sub.imshow(np.transpose(resdat), vmax=impeak,
                     origin='lower', extent=[IMsize, -IMsize, -IMsize, IMsize])
    sub.set_xlabel('RA offset (mas)')
    sub.set_ylabel('Dec offset (mas)')
    cb = pl.colorbar(ims)
    cb.set_label('Jy/beam')
    pl.savefig('%s.png' % Rfile)
