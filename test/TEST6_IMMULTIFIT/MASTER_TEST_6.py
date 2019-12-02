import os
import numpy as np
import scipy as sp
import pylab as pl

# TEST 6: TWO OFFSET GAUSSIANS
#         FIT DONE WITH GRIDDED VISIBILITIES

#pylint: disable=undefined-variable

# What to do:
# DoSimObs = False
# DoFit = True
# casaexe = 'casa --nologger'

# Array to use:
arrayconfig = 'alma.out10.cfg'

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

# For test 3:
imname = 'Offset_Gaussians'
# type   flux    major       minor       PA  freq   spec.ind.    Direction
s = [['gaussian', 1.0, '0.2arcsec', '0.2arcsec', '0.0deg', Nu, 0.0, "J2000 10h00m00.0s -30d00m00.0s"],
     ['gaussian', 0.5, '0.4arcsec', '0.4arcsec', '0.0deg', Nu, 0.0, "J2000 10h00m00.0s -30d00m01.5s"]]

doff = -1.5  # Dec offset of Gaussian (arcsec.)
config = '.'.join(arrayconfig.split('.')[:-1])

if DoSimObs:
    print 'Generating %s' % imname
    cl.done()

    for si in s:
        cl.addcomponent(dir=si[7], flux=si[1], fluxunit='Jy', freq=si[5], shape=si[0],
                        majoraxis=si[2], minoraxis=si[3], positionangle=si[4],
                        spectrumtype='spectral index', index=si[6])

    ia.fromshape(imname+'.model', [Npix, Npix, 1, Nnu], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad = qa.convert(qa.quantity(cell), "rad")['value']
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    rdi, ddi = s[0][7].split()[1:]
    cs.setreferencevalue([qa.convert(rdi, 'rad')['value'],
                          qa.convert(ddi, 'rad')['value']], type="direction")
    cs.setreferencevalue(s[0][5], 'spectral')
    cs.setincrement(ChW, 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    ia.close()

    simobserve(project=imname, skymodel=imname+'.model',
               totaltime='300s', antennalist=arrayconfig)

if DoFit:
    print("---------------------------------------------")
    print("TEST 6")
    print("---------------------------------------------")
    pl.ioff()
    Cfile = 'TEST6.CLEAN'
    clearcal('%s/%s.%s.noisy.ms' % (imname, imname, config))
    os.system('rm -rf %s.*' % Cfile)

    clean('%s/%s.%s.noisy.ms' % (imname, imname, config),
          imagename=Cfile, cell=cell, weighting='briggs', robust=0.5,
          psfmode='hogbom', ftmachine='ft', imagermode='', imsize=Npix, niter=0)

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
    tempfile = open('STEP6_FIT.py', 'w')

    # p[0] -> flux of first Gaussian
    # p[1] -> Major axis of first Gaussian
    # p[2] -> Pos. Offset of second Gaussian
    # p[3] -> Flux density of second Gaussian
    # p[4] -> Size of second Gaussian
    Gvar = '0.,p[2],p[3],p[4],1.,0.'
    D1var = '0.,0.,p[0],p[1],1.,0.'

    shapes = [str(si[0]) for si in s]
    fluxes = [float(si[1]) for si in s]
    sizes = [float(si[2].split('a')[0]) for si in s]
    minors = [float(si[3].split('a')[0]) for si in s]
    PAs = [float(si[4].split('d')[0]) for si in s]

    Pars = [fluxes[0], sizes[0], doff, fluxes[1], sizes[1]]
    BiasFac = [0.7, 1.3, 0.75, 1.2, 1.3]
    Trues = Pars

    string = 'S = [%.3e, %.3e, %.3e, %.3e, %.3e]' % tuple(Pars)
    print >> tempfile, string

    string = "immname = '%s'" % ('%s.residual' % Cfile)
    print >> tempfile, string

    string = "psfname = '%s'" % ('%s.psf' % Cfile)
    print >> tempfile, string

    for i in range(len(shapes)):
        if shapes[i] == 'disk':
            shapes[i] = 'disc'
        if shapes[i] == 'gaussian':
            shapes[i] = 'Gaussian'

    string = 'modelshape = [%s]' % (','.join(["'" + shi +"'" for shi in shapes]))
    print >> tempfile, string

    NuF = float(Nu.split('G')[0])*1.e9
    string = "modvars = ['%s','%s']" % (D1var, Gvar)
    print >> tempfile, string

    string = "pini = [%.3e, %.3e, %.3e, %.3e, %.3e]" % tuple(
        [Pars[pi]*BiasFac[pi] for pi in range(len(Pars))])
    print >> tempfile, string

    string = "parbound = [[0.,None], [0.,None], [-15,15], [0.,None], [0.,None]]"
    print >> tempfile, string

    lines = open('test6.py')
    for l in lines.readlines():
        print >> tempfile, l[:-1]
    lines.close()
    tempfile.close()

    os.system('%s -c STEP6_FIT.py' % casaexe)
