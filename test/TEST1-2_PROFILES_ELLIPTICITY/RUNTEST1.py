import os
import numpy as np
import pylab as pl

DoSimObs = False
DoFit = True
casaexe = 'casa --nologger'

arrayconfig = 'alma.out10.cfg'

# Center Freq:
Nu = '50.0GHz'
# Image size (pixels)
Npix = 1000
# Pixel size:
cell = '0.01arcsec'
cellN = 0.01
# Number of frequency channels:
Nnu = 100
# Channel width:
ChW = '0.01GHz'

# For test 1:
imname = 'Disc'
# type   flux    major       minor       PA  freq   spec.ind.    Direction
si = ['disk', 1.0, '0.2arcsec', '0.1arcsec', '60.0deg', Nu, 1.0, "J2000 10h00m00.0s -30d00m00.0s"]

# For test 2:
Diameter = 1.0
DiamDelta = 0.2
Sigma = 0.25

config = '.'.join(arrayconfig.split('.')[:-1])
if DoSimObs:
    # TEST 1: ELLIPTICITY AND CONTINUUM FIT:
    print('Generating %s' % imname)
    cl.done()
    cl.addcomponent(dir=si[7], flux=si[1], fluxunit='Jy', freq=si[5], shape=si[0],
                    majoraxis=si[2], minoraxis=si[3], positionangle=si[4],
                    spectrumtype='spectral index', index=si[6])

    ia.fromshape(imname+'.model', [Npix, Npix, 1, Nnu], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad = qa.convert(qa.quantity(cell), "rad")['value']
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    rdi, ddi = si[7].split()[1:]
    cs.setreferencevalue([qa.convert(rdi, 'rad')['value'],
                          qa.convert(ddi, 'rad')['value']], type="direction")
    cs.setreferencevalue(si[5], 'spectral')
    cs.setincrement(ChW, 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    ia.close()

    simobserve(project=imname, skymodel=imname+'.model',
               totaltime='300s', antennalist=arrayconfig)

    # TEST 2: RADIAL PROFILES
    IMsize = Npix*cellN
    NuF = float(Nu.split('G')[0])*1.e9

    LAS = 1./(IMsize/3600.*np.pi/180.)*3.e8/NuF

    X = np.linspace(-IMsize/2., IMsize/2., Npix)
    XX = np.power(np.outer(X, np.ones(Npix)), 2.)
    RSq = XX + np.transpose(XX)
    R = np.sqrt(RSq)

if DoFit:
    print("---------------------------------------------")
    print("TEST 1")
    print("---------------------------------------------")
    # TEST 1:
    # Note: the flux density of the disc at 50GHz is 0.2435 Jy, and not 1Jy. For some reason,
    # the cl CASA tool doesn't put the right flux density of the disc component (you can
    # check it with plotms and/or adding up the flux density from the pixels of the image model).
    os.system('rm -rf *.png')
    with open('STEP1_FIT.py', 'w') as tempfile:
        print("from NordicARC import uvmultifit as uvm", file=tempfile)

        size = float(si[2].split('a')[0])
        minor = float(si[3].split('a')[0])
        PA = float(si[4].split('d')[0])
        string = "S = [%.3e, %.3e, %.3e, %.3e, %.3e]" % (0.2435, si[6], size, minor/size, PA)
        print(string, file=tempfile)

        string = "visname = '%s'" % ("%s/%s.%s.noisy.ms" % (imname, imname, config))
        print(string, file=tempfile)

        if si[0] == 'disk':
            si[0] = 'disc'
        string = "modelshape = '%s'" % (si[0])
        print(string, file=tempfile)

        NuF = float(Nu.split('G')[0])*1.e9
        string = "modvars = '0,0, p[0]*(nu/%.4e)**p[1], p[2], p[3], p[4]'" % NuF
        print(string, file=tempfile)

        string = "pini = [%.3e, %.3e, %.3e, %.3e, %.3e]" % (0.8, 0., size*1.2, minor/size*0.8, 45.)
        print(string, file=tempfile)

        string = "parbound = [[0., None], [-2.,2.], [0.,None], [0.1,0.9], [0.,180.]]"
        print(string, file=tempfile)

        lines = open('test1.py')
        for l in lines.readlines():
            print(l[:-1], file=tempfile)
        lines.close()

    pl.ioff()
    Cfile = 'TEST1.CLEAN'
    Rfile = 'TEST1.RESIDUALS'
    clearcal('%s/%s.%s.noisy.ms' % (imname, imname, config))

    os.system('rm -rf %s.*' % Cfile)
    os.system('rm -rf %s.*' % Rfile)

    tclean('%s/%s.%s.noisy.ms' % (imname, imname, config),
           imagename=Cfile, cell=cell,
           imsize=2*Npix, niter=0)

    ia.open('%s.image' % Cfile)
    resdat = ia.getchunk()[:, :, 0, 0]
    ia.close()
    fig = pl.figure()
    sub = fig.add_subplot(111)
    IMsize = cellN*Npix/2.
    ims = sub.imshow(np.transpose(resdat), origin='lower',
                     extent=[IMsize, -IMsize, -IMsize, IMsize])
    sub.set_xlabel('RA offset (arcsec)')
    sub.set_ylabel('Dec offset (arcsec)')
    cb = pl.colorbar(ims)
    cb.set_label('Jy/beam')
    pl.savefig('%s.png' % Cfile)
    impeak = np.max(resdat)

    # os.system('%s -c STEP1_FIT.py' % casaexe)
    exec(open("STEP1_FIT.py").read())

    os.system('rm -rf %s.*' % Rfile)
    tclean('%s/%s.%s.noisy.ms' % (imname, imname, config),
           imagename=Rfile, cell=cell,
           imsize=2*Npix, niter=0)

    impeak *= 0.01

    ia.open('%s.image' % Rfile)
    resdat = ia.getchunk()[:, :, 0, 0]
    ia.close()
    fig = pl.figure()
    sub = fig.add_subplot(111)
    IMsize = cellN*Npix/2.
    ims = sub.imshow(np.transpose(resdat), origin='lower',
                     vmax=impeak, extent=[IMsize, -IMsize, -IMsize, IMsize])
    sub.set_xlabel('RA offset (arcsec)')
    sub.set_ylabel('Dec offset (arcsec)')
    cb = pl.colorbar(ims)
    cb.set_label('Jy/beam')
    pl.savefig('%s.png' % Rfile)
