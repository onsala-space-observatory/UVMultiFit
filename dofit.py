import os
import numpy as np

from casatools import image

from NordicARC import uvmultifit as uvm

import matplotlib.pyplot as plt

vis = 'Disc/Disc.alma.out10.noisy.ms'
model = ['disc']
Nu = '50.0GHz'
NuF = float(Nu.split('G')[0])*1.e9
modvars = "0,0, p[0]*(nu/%.4e)**p[1], p[2], p[3], p[4]" % NuF
modbounds = [[0.0, None], [-2.0, 2.0], [0.0, None], [0.1, 0.9], [0.0, 180.0]]

si = ['disc', 1.0, '0.2arcsec', '0.1arcsec', '60.0deg',
      Nu, 1.0, "J2000 10h00m00.0s -30d00m00.0s"]
size = float(si[2].split('a')[0])
minor = float(si[3].split('a')[0])
initial = [0.8, 0.0, size*1.2, minor/size*0.8, 45.0]

if False:
    from casatasks import clearcal
    from casatasks import tclean

    Cfile = 'TEST1.CLEAN'
    Rfile = 'TEST1.RESIDUALS'
    clearcal(vis)

    os.system('rm -rf %s.*' % Cfile)
    os.system('rm -rf %s.*' % Rfile)

    Npix = 1000
    cell = '0.01arcsec'
    cellN = float(cell.replace("arcsec", ""))
    tclean(vis, imagename=Cfile, cell=cell, imsize=2*Npix, niter=0)

    ia = image()
    ia.open('%s.image' % Cfile)
    resdat = ia.getchunk()[:, :, 0, 0]
    ia.close()

    fig = plt.figure()
    sub = fig.add_subplot(111)
    IMsize = cellN*Npix/2.
    ims = sub.imshow(np.transpose(resdat), origin='lower',
                     extent=[IMsize, -IMsize, -IMsize, IMsize])
    sub.set_xlabel('RA offset (arcsec)')
    sub.set_ylabel('Dec offset (arcsec)')
    cb = plt.colorbar(ims)
    cb.set_label('Jy/beam')
    # plt.savefig('%s.png' % Cfile)
    plt.show()
    impeak = np.max(resdat)

pbeam = False
r = uvm.uvmultifit(vis=vis, spw='0',
                   model=model, OneFitPerChannel=False,
                   var=modvars, write='residuals',
                   p_ini=initial, pbeam=pbeam,
                   bounds=modbounds)

S = np.array([2.435e-01, 1.000e+00, 2.000e-01, 5.000e-01, 6.000e+01])

print(f"disc flux at 50GHz (Jy): {r['Parameters'][0]:7.3f} +/- {r['Uncertainties'][0]:.3f}, true: {S[0]:7.3f}")
print(f"disc spectral index:     {r['Parameters'][1]:7.3f} +/- {r['Uncertainties'][1]:.3f}, true: {S[1]:7.3f}")
print(f"disc size (as):          {r['Parameters'][2]:7.3f} +/- {r['Uncertainties'][2]:.3f}, true: {S[2]:7.3f}")
print(f"disc axis ratio:         {r['Parameters'][3]:7.3f} +/- {r['Uncertainties'][3]:.3f}, true: {S[3]:7.3f}")
print(f"disc pos.angle (deg):    {r['Parameters'][4]:7.3f} +/- {r['Uncertainties'][4]:.3f}, true: {S[4]:7.3f}")
