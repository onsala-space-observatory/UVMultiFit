import numpy as np
import scipy as sp
import pylab as pl
import os
import sys

# What to do:
# DoSimObs = False
# DoFit = True
# casaexe = 'casa --nologger'

imname = 'WideField'
arrayconfig = 'vla.d.cfg'

phref = "J2000 10h00m00.08s 30d00m02.0s"
Nu = '1.5GHz'
Npix = 4000
cell = '0.9arcsec'
cellN = 0.9  # in arcsec
Nnu = 10
ChW = '0.01GHz'

# Will make a grid of NxN, with equidistant sources separated "incr." degrees.
N = 4
incr = 0.20  # in degrees

# The flux densities of the components:
# Fi = [1.0, 1.5, 1.2, 0.8, 0.9, 2.0, 1.4, 1.0, 1.5, 1.2, 0.8, 0.9, 2.0, 1.4, 0.75, 0.8]
Fi = np.ones(N*N)

config = '.'.join(arrayconfig.split('.')[:-1])
visname = '%s/%s.%s.noisy.ms' % (imname, imname, config)

if DoSimObs:
    ia.fromshape(imname+'.model', [Npix, Npix, 1, Nnu], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad = qa.convert(qa.quantity(cell), "rad")['value']
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    rdi, ddi = phref.split()[1:]
    cs.setreferencevalue([qa.convert(rdi, 'rad')['value'],
                          qa.convert(ddi, 'rad')['value']], type="direction")
    cs.setreferencevalue(Nu, 'spectral')
    cs.setincrement(ChW, 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.close()

    ia.open(imname+'.model')
    data = ia.getchunk()
    for i in range(N):
        for j in range(N):
            ipix = int(Npix/2 + (-float(N-1)/2. + i)*incr*3600./cellN)
            jpix = int(Npix/2 + (-float(N-1)/2. + j)*incr*3600./cellN)
            print i, j, ipix, jpix
            data[ipix, jpix, :, :] = Fi[i+j]

    ia.putchunk(data)
    ia.close()

    simobserve(project=imname, skymodel=imname+'.model',
               indirection='phref',
               totaltime='1200s', antennalist=arrayconfig)

    print "Done"

if DoFit:
    print("---------------------------------------------")
    print("TEST 4")
    print("---------------------------------------------")
    pl.ioff()
    Cfile = 'TEST4.CLEAN'
    Rfile = 'TEST4.RESIDUALS'

    tempfile = open('STEP4_FIT.py', 'w')

    field = ''

    if False:
        clearcal(vis=visname, addmodel=True)

        os.system('rm -rf %s.*' % Cfile)
        clean(vis=visname, field=field, imagename=Cfile, cell=cell, imsize=Npix,
              interactive=False, niter=0, imagermode='mosaic', pbcor=True, phasecenter=phref)

        ia.open('%s.image' % Cfile)
        resdat = ia.getchunk()[:, :, 0, 0]
        ia.close()

        ia.open('%s.flux' % Cfile)
        flxdat = ia.getchunk()[:, :, 0, 0]
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
        #  RAs = [int((-float(N-1)/2. + i)*incr*3600./cellN)*cellN for i in range(N)]
        #  Decs = [int((-float(N-1)/2. + i)*incr*3600./cellN)*cellN for i in range(N)]

    model = ['delta' for i in range(N*N)]
    RAs = [incr*3600.*((N-1.)/2.-i) for i in range(N)]
    Decs = [incr*3600.*(-(N-1.)/2.+i) for i in range(N)]

    modvars = []
    pini = []
    S = []
    i = 0
    for ra in RAs:
        for dec in Decs:
            modvars.append('%8.3f+p[%i], %8.3f+p[%i], p[%i]' % (ra, 3*i, dec, 3*i+1, 3*i+2))
            # modvars.append('%8.3f, %8.3f, p[%i]'%(ra,dec,i))

            pini += [0., 0., 0.8]
            S += [0., 0., 1.]
            i += 1

    print >> tempfile, "field = '%s'" % field
    print >> tempfile, "msout = '%s'" % visname
    print >> tempfile, "model = [%s]" % (','.join(["'%s'" % s for s in model]))
    print >> tempfile, "modvars = [%s]" % (','.join(["'%s'" % s for s in modvars]))
    print >> tempfile, "pini = [%s]" % (','.join(map(str, pini)))
    print >> tempfile, "phref = '%s'" % phref
    print >> tempfile, "S = [%s]" % (','.join(map(str, S)))
    print >> tempfile, "bounds = [%s]" % (','.join(['[-5.,5], [-5.,5.], [0.,None]' for pi in range(N*N)]))

    iff = open('test4.py')
    for line in iff.readlines():
        print >> tempfile, line[:-1]
    iff.close()
    tempfile.close()

    os.system('%s -c STEP4_FIT.py' % casaexe)
    if True:
        os.system('rm -rf %s.*' % Rfile)

        # ANOTHER WAY OF COMPUTING THE RESIDUALS, AND KEEP THE MODEL
        # IN THE MODEL COLUMN:
        #  ms.open(visname,nomodify=False)
        #  ms.selectinit(datadescid=0)
        #  data = ms.getdata(['model_data','data','corrected_data'])
        #  data['corrected_data'][:] = data['data']-data['model_data']
        #  ms.putdata(data)
        #  ms.close()

        clean(vis=visname, imagename=Rfile, cell=cell, imsize=Npix, field=field,
              interactive=False, niter=0, imagermode='mosaic', pbcor=True, phasecenter=phref)

        ia.open('%s.image' % Rfile)
        resdat = ia.getchunk()[:, :, 0, 0]
        ia.close()
        ia.open('%s.flux' % Rfile)
        flxdat = ia.getchunk()[:, :, 0, 0]
        ia.close()

        try:
            impeak *= 0.01
        except:
            impeak = np.max(resdat)

        fig = pl.figure()
        sub = fig.add_subplot(111)
        IMsize = cellN*Npix/2.
        ims = sub.imshow(np.transpose(resdat/flxdat), vmax=impeak,
                         origin='lower', extent=[IMsize, -IMsize, -IMsize, IMsize])
        sub.set_xlabel('RA offset (arcsec)')
        sub.set_ylabel('Dec offset (arcsec)')
        cb = pl.colorbar(ims)
        cb.set_label('Jy/beam')
        pl.savefig('%s.png' % Rfile)
