import numpy as np
import scipy as sp
import pylab as pl
import os

# What to do:
# DoSimObs = False
# DoFit = True
# casaexe='casa --nologger'

# Array to use:
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
si = ['disk', 1.0, '0.2arcsec', '0.1arcsec', '60.0deg',
      Nu, 1.0, "J2000 10h00m00.0s -30d00m00.0s"]

# For test 2:
Diameter = 1.0
DiamDelta = 0.2
Sigma = 0.25
imname2 = ['GaussianRing', 'Gaussian', 'ring', 'sphere', 'bubble', 'expo', 'power-2']

config = '.'.join(arrayconfig.split('.')[:-1])
if DoSimObs:
    # TEST 1: ELLIPTICITY AND CONTINUUM FIT:
    print 'Generating %s' % imname
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

    for j, modnam in enumerate(imname2):
        print 'Generating %s' % modnam
        os.system('rm -rf %s.model' % modnam)
        os.system('cp -r %s.model  %s.model' % (imname, modnam))

        ia.open(modnam+'.model')
        data = ia.getchunk()
        for i in range(Nnu):
            SouRad = (Diameter+DiamDelta*i/Nnu)/2.
            mask = R >= SouRad
            if j == 0:
                data[:, :, 0, i] = np.exp(-np.power(R-SouRad, 2.)/(2.*Sigma**2.))
                data[:, :, 0, i] /= np.sum(data[:, :, 0, i])
            elif j == 1:
                data[:, :, 0, i] = np.exp(-RSq/(2.*(2.*SouRad/2.35482)**2.))
                data[:, :, 0, i] /= np.sum(data[:, :, 0, i])
            elif j == 2:
                data[:, :, 0, i] = np.exp(-np.power(R-SouRad, 2.)/(2.*(2.*cellN)**2.))
                data[:, :, 0, i] /= np.sum(data[:, :, 0, i])
            elif j == 3:
                data[:, :, 0, i] = SouRad**2. - RSq
                data[mask, 0, i] = 0.0
                data[:, :, 0, i] = np.sqrt(data[:, :, 0, i])
                data[:, :, 0, i] /= np.sum(data[:, :, 0, i])
            elif j == 4:
                data[:, :, 0, i] = SouRad**2. - RSq
                data[mask, 0, i] = 0.0
                data[np.logical_not(mask), 0, i] = 1. / np.sqrt(data[np.logical_not(mask), 0, i])
                data[:, :, 0, i] /= np.sum(data[:, :, 0, i])
            elif j == 5:
                data[:, :, 0, i] = np.exp(-R/SouRad*np.log(2.))
                data[:, :, 0, i] /= np.sum(data[:, :, 0, i])
            elif j == 6:
                data[:, :, 0, i] = 1/(RSq + SouRad**2.)
                # Normalize to flux-density from 0 to R:
                data[:, :, 0, i] /= 2.*np.pi * np.sum(data[np.logical_not(mask), 0, i])

        ia.putchunk(data)
        ia.close()

        simobserve(project=modnam, skymodel=modnam+'.model',
                   totaltime='300s', antennalist=arrayconfig)

        if j == 6:  # Remove window effect in the power-2 radial profile:
            flagdata('%s/%s.%s.ms' % (modnam, modnam, config),
                     mode='manual', uvrange='<%.1fm' % LAS)

            flagdata('%s/%s.%s.noisy.ms' % (modnam, modnam, config),
                     mode='manual', uvrange='<%.1fm' % LAS)

if DoFit:
    print("---------------------------------------------")
    print("TEST 1")
    print("---------------------------------------------")
    # TEST 1:
    # Note: the flux density of the disc at 50GHz is 0.2435 Jy, and not 1Jy. For some reason,
    # the cl CASA tool doesn't put the right flux density of the disc component (you can
    # check it with plotms and/or adding up the flux density from the pixels of the image model).
    os.system('rm -rf *.png')
    tempfile = open('STEP1_FIT.py', 'w')

    size = float(si[2].split('a')[0])
    minor = float(si[3].split('a')[0])
    PA = float(si[4].split('d')[0])
    string = "S = [%.3e, %.3e, %.3e, %.3e, %.3e]" % (0.2435, si[6], size, minor/size, PA)
    print >> tempfile, string

    string = "visname = '%s'" % ("%s/%s.%s.noisy.ms" % (imname, imname, config))
    print >> tempfile, string

    if si[0] == 'disk':
        si[0] = 'disc'
    string = "modelshape = '%s'" % (si[0])
    print >> tempfile, string

    NuF = float(Nu.split('G')[0])*1.e9
    string = "modvars = '0,0, p[0]*(nu/%.4e)**p[1], p[2], p[3], p[4]'" % NuF
    print >> tempfile, string

    string = "pini = [%.3e, %.3e, %.3e, %.3e, %.3e]" % (0.8, 0., size*1.2, minor/size*0.8, 45.)
    print >> tempfile, string

    string = "parbound = [[0., None], [-2.,2.], [0.,None], [0.1,0.9], [0.,180.]]"
    print >> tempfile, string

    lines = open('test1.py')
    for l in lines.readlines():
        print >> tempfile, l[:-1]
    lines.close()
    tempfile.close()

    pl.ioff()
    Cfile = 'TEST1.CLEAN'
    Rfile = 'TEST1.RESIDUALS'
    clearcal('%s/%s.%s.noisy.ms' % (imname, imname, config))

    os.system('rm -rf %s.*' % Cfile)
    os.system('rm -rf %s.*' % Rfile)

    cleancmd = "clean('%s/%s.%s.noisy.ms', imagename='%s', cell='%s', imsize=%d, niter=0)" % \
               (imname, imname, config, Cfile, cell, 2*Npix)
    print cleancmd
    clean('%s/%s.%s.noisy.ms' % (imname, imname, config),
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

    os.system('%s -c STEP1_FIT.py' % casaexe)

    os.system('rm -rf %s.*' % Rfile)
    cleancmd = "clean('%s/%s.%s.noisy.ms', imagename='%s', cell='%s', imsize=%d, niter=0)" % \
               (imname, imname, config, Rfile, cell, 2*Npix)
    print cleancmd
    clean('%s/%s.%s.noisy.ms' % (imname, imname, config),
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

    # TEST 2:
    for j, modnam in enumerate(imname2):
        print("---------------------------------------------")
        print("TEST2: %s" % (modnam))
        print("---------------------------------------------")
        tempfile = open('STEP2_FIT_%s.py' % modnam, 'w')
        if modnam == 'GaussianRing':
            var = 'shapevar = \'p[3],p[4],p[0],p[1],1.0,0.0,p[2]\''
            pi = 'pi = [0.8, %.3e,%.3e,0.,0.]' % (Diameter*0.8, Sigma*1.2)
            bb = 'parbound = [[0.,None],[0.,%.3e],[%.3e,%.3e],[-1.,1.],[-1.,1.]]' % (
                Diameter*1.5, Sigma*0.8, Sigma*1.5)
        else:
            var = 'shapevar = \'p[2],p[3],p[0],p[1],1.0,0.0\''
            pi = 'pi = [0.8, %.3e,0.,0.]' % (Diameter*1.2)
            bb = 'parbound = [[0.,None],[0.,%.3e],[-1.,1.],[-1.,1.]]' % (
                Diameter*1.5)

        string = 'Diameter = %.4e' % Diameter
        print >> tempfile, string
        string = 'DiamDelta = %.4e' % DiamDelta
        print >> tempfile, string
        string = 'Sigma = %.4e' % Sigma
        print >> tempfile, string

        string = 'visname = \'%s/%s.%s.noisy.ms\'' % (modnam, modnam, config)
        print >> tempfile, string
        string = 'modshape = \'%s\'' % modnam
        print >> tempfile, string
        print >> tempfile, var
        print >> tempfile, pi
        print >> tempfile, bb

        lines = open('test2.py')
        for l in lines.readlines():
            print >> tempfile, l[:-1]
        lines.close()
        tempfile.close()

        Cfile = 'TEST2.%s.CLEAN' % modnam
        Rfile = 'TEST2.%s.RESIDUALS' % modnam
        os.system('rm -rf %s.*' % Cfile)
        os.system('rm -rf %s.*' % Rfile)

        clearcal('%s/%s.%s.noisy.ms' % (modnam, modnam, config))

        cleancmd = "clean('%s/%s.%s.noisy.ms', imagename='%s', cell='%s', imsize=%d, niter=0)" % \
                   (modnam, modnam, config, Cfile, cell, 2*Npix)
        print cleancmd
        clean('%s/%s.%s.noisy.ms' % (modnam, modnam, config),
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
        os.system('%s -c STEP2_FIT_%s.py' % (casaexe, modnam))
        os.system('rm -rf %s.*' % Rfile)
        cleancmd = "clean('%s/%s.%s.noisy.ms', imagename='%s', cell='%s', imsize=%d, niter=0)" % \
                   (modnam, modnam, config, Cfile, cell, 2*Npix)
        print cleancmd
        clean('%s/%s.%s.noisy.ms' % (modnam, modnam, config),
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

    os.system('rm -rf test2.dat temp.py')
    os.system('cat ' + ' '.join(['test2.%s.dat' % im for im in imname2]) + ' > test2.dat')
    os.system('rm ' + ' '.join(['test2.%s.dat' % im for im in imname2]))
