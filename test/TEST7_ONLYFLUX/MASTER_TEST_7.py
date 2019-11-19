import numpy as np
import scipy as sp
import pylab as pl
import os

# TEST 7: SEVERAL SOURCES. ONLY FIT THE FLUX DENSITY.

# What to do:
# DoSimObs = False
# DoFit = True
# casaexe='casa --nologger'

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
imname = 'Only_Flux'
#      type       flux   major         minor         PA     freq  spec.ind.       Direction
s = [['disk',     1.0,  '0.5arcsec',  '0.5arcsec',  '0.0deg', Nu, 0.0, "J2000 10h00m00.0s   -30d00m00.0s"],
     ['disk',     0.75, '0.25arcsec', '0.25arcsec', '0.0deg', Nu, 0.0, "J2000 10h00m00.154s -30d00m00.0s"],
     ['gaussian', 0.25, '0.4arcsec',  '0.4arcsec',  '0.0deg', Nu, 0.0, "J2000 10h00m00.0s   -30d00m02.0s"],
     ['gaussian', 2.5,  '0.1arcsec',  '0.1arcsec',  '0.0deg', Nu, 0.0, "J2000 09h59m59.692s -30d00m00.0s"],
     ['gaussian', 1.5,  '0.4arcsec',  '0.4arcsec',  '0.0deg', Nu, 0.0, "J2000 09h59m59.846s -29d59m58.0s"],
     ['disk',     0.15, '0.85arcsec', '0.85arcsec', '0.0deg', Nu, 0.0, "J2000 09h59m59.692s -30d00m02.0s"]
     ]

OFFSETS = [[0., 0.], [2., 0.], [0., -2.], [-4., 0.], [-2., 2], [-4., -2.]]
BiasFac = [0.6, 1.6, 1.8, 0.2, 0.8, 1.2]

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
    print("TEST 7")
    print("---------------------------------------------")
    tempfile = open('STEP7_FIT.py', 'w')

    # p[0] -> Flux of first disc
    # p[1] -> Flux of second disc
    # p[2] -> Flux of Gaussian
    shapes = [str(si[0]) for si in s]
    fluxes = [float(si[1]) for si in s]
    sizes = [float(si[2].split('a')[0]) for si in s]
    minors = [float(si[3].split('a')[0]) for si in s]
    PAs = [float(si[4].split('d')[0]) for si in s]

    Pvars = []

    for ci, fli in enumerate(fluxes):
        Pvars.append('%.3e, %.3e, p[%i],%.3e,1.,0.' % (
            OFFSETS[ci][0], OFFSETS[ci][1], ci, sizes[ci]))

    Pars = fluxes  # [fluxes[0],fluxes[1],fluxes[2]]

# For some reason, the disc flux in CASA is wrong by a factor ~0.2477 (??)
# You can check this, e.g. on the model image with the viewer.
    for ci, comp in enumerate(shapes):
        if comp == 'disk':
            Pars[ci] *= 0.2477

    Trues = Pars

    string = ('S = [%.3e' + ',%.3e'*(len(Pars)-1) + ']') % tuple(Pars)
    print >> tempfile, string

    string = "visname = '%s'" % ('%s/%s.%s.noisy.ms' % (imname, imname, config))
    print >> tempfile, string

    for i in range(len(shapes)):
        if shapes[i] == 'disk':
            shapes[i] = 'disc'
        if shapes[i] == 'gaussian':
            shapes[i] = 'Gaussian'

    string = "modelshape = [%s]" % (','.join(["'" + shi +"'" for shi in shapes]))
    print >> tempfile, string

    NuF = float(Nu.split('G')[0])*1.e9
    string = ("modvars = ['%s'" + ",'%s'" * (len(Pars)-1) + ']') % tuple(Pvars)
    print >> tempfile, string

    string = ("pini = [%.3e" + ",%.3e"*(len(Pars)-1) + "]") % \
                      tuple([Pars[pi]*BiasFac[pi] for pi in range(len(Pars))])
    print >> tempfile, string

    string = 'parbound = None'
    print >> tempfile, string

    lines = open('test7.py')
    for l in lines.readlines():
        print >> tempfile, l[:-1]
    lines.close()
    tempfile.close()

    pl.ioff()
    Cfile = 'TEST7.CLEAN'
    Rfile = 'TEST7.RESIDUALS'
    clearcal('%s/%s.%s.noisy.ms' % (imname, imname, config))
    os.system('rm -rf %s.*' % Cfile)
    os.system('rm -rf %s.*' % Rfile)

    tclean('%s/%s.%s.noisy.ms' % (imname, imname, config),
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
    sub.set_xlabel('RA offset (arcsec)')
    sub.set_ylabel('Dec offset (arcsec)')
    cb = pl.colorbar(ims)
    cb.set_label('Jy/beam')
    pl.savefig('%s.png' % Cfile)

    impeak = np.max(resdat)

    os.system('%s -c STEP7_FIT.py' % casaexe)

    os.system('rm -rf %s.*' % Rfile)
    tclean('%s/%s.%s.noisy.ms' % (imname, imname, config),
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
    sub.set_xlabel('RA offset (arcsec)')
    sub.set_ylabel('Dec offset (arcsec)')
    cb = pl.colorbar(ims)
    cb.set_label('Jy/beam')
    pl.savefig('%s.png' % Rfile)
