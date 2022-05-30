import os
import numpy as np
import pylab as pl

# TEST 3: DISC WITH A HOLE PLUS AN OFFSET GAUSSIAN

#pylint: disable=undefined-variable

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
imname = 'Discs_Gaussian'
# type   flux    major       minor       PA  freq   spec.ind.    Direction
s = [['disk',   1.0, '0.5arcsec', '0.5arcsec', '0.0deg', Nu, 0.0, "J2000 10h00m00.0s -30d00m00.0s"],
     ['disk', -0.25, '0.25arcsec', '0.25arcsec', '0.0deg', Nu, 0.0, "J2000 10h00m00.0s -30d00m00.0s"],
     ['gaussian', 0.5, '0.4arcsec', '0.4arcsec', '0.0deg', Nu, 0.0, "J2000 10h00m00.0s -30d00m02.0s"]]

doff = -2.  # Dec offset of Gaussian (arcsec.)
config = '.'.join(arrayconfig.split('.')[:-1])
if DoSimObs:
    # TEST 1: MULTI-COMPONENT CORRELATED MODEL.
    print('Generating %s' % imname)
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
    print("TEST 3")
    print("---------------------------------------------")

    # p[0] -> flux of outer disc
    # p[1] -> Major axis of outer disc
    # p[2] -> Ratio of outer-to-inner disc
    # p[3] -> Pos. Offset of Gaussian
    # p[4] -> Flux density of Gaussian
    # p[5] -> Size of Gaussian

    Gvar = '0.,p[3],p[4],p[5],1.,0.'
    D1var = '0.,0.,p[0],p[1],1.,0.'
    D2var = '0.,0.,-p[0]*p[2]**2.,p[1]*p[2],1.,0.'

    shapes = [str(si[0]) for si in s]
    fluxes = [float(si[1]) for si in s]
    sizes = [float(si[2].split('a')[0]) for si in s]
    minors = [float(si[3].split('a')[0]) for si in s]
    PAs = [float(si[4].split('d')[0]) for si in s]

    Pars = [fluxes[0], sizes[0], sizes[1]/sizes[0], doff, fluxes[2], sizes[2]]
    BiasFac = [0.7, 1.3, 0.8, 0.75, 1.2, 1.3]
    Trues = Pars

    # For some reason, the disc flux in CASA is wrong by a factor ~0.2477 (??)
    # You can check this, e.g. on the model image with the viewer.
    Trues[0] *= 0.2477

    string = 'S = [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]' % tuple(Pars)
    with open('STEP3_FIT.py', 'w') as tempfile:
        print(string, file=tempfile)

        string = "visname = '%s'" % ('%s/%s.%s.noisy.ms' % (imname, imname, config))
        print(string, file=tempfile)

        for idx, shape in enumerate(shapes):
            if shape == 'disk':
                shapes[idx] = 'disc'
            if shape == 'gaussian':
                shapes[idx] = 'Gaussian'

        string = 'modelshape = [%s]' % (','.join(["'" + shi + "'" for shi in shapes]))
        print(string, file=tempfile)

        NuF = float(Nu.split('G')[0])*1.e9
        string = "modvars = ['%s','%s','%s']" % (D1var, D2var, Gvar)
        print(string, file=tempfile)

        string = 'pini = [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]' % tuple(
            [Pars[pi]*BiasFac[pi] for pi in range(len(Pars))])
        print(string, file=tempfile)

        string = 'parbound = [[0.,None], [0.,None], [0.,None], [-15,15], [0.,None], [0.,None]]'
        print(string, file=tempfile)

        lines = open('test3.py')
        for l in lines.readlines():
            print(l[:-1], file=tempfile)
        lines.close()

    pl.ioff()
    Cfile = 'TEST3.CLEAN'
    Rfile = 'TEST3.RESIDUALS'
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

    # os.system('%s -c STEP3_FIT.py' % casaexe)
    exec(open("STEP3_FIT.py").read())

    # ms.open('%s/%s.%s.noisy.ms'%(imname,imname,config),nomodify=False)
    # ms.selectinit(datadescid=0)
    # data = ms.getdata(['corrected_data','model_data'])
    # data['corrected_data'][:] = data['model_data']
    # ms.putdata(data)
    # ms.close()
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
