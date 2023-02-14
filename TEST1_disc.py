import os
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

from casatools import componentlist
from casatools import image
from casatools import quanta
from casatasks import simobserve

from simulation import Simulation

from NordicARC import uvmultifit as uvm

cl = componentlist()
ia = image()
qa = quanta()

sim = Simulation(array_config="alma.out10.cfg",
                 center_freq=50.0,
                 pixels=1000,
                 cell_size=0.01,
                 freq_channels=100,
                 channel_width=0.01,
                 image_name="Disc",
                 total_time=300)

Nu = f"{sim.center_freq}GHz"
cell = f"{sim.cell_size}arcsec"
ChW = f"{sim.channel_width}GHz"
totaltime = f"{sim.total_time}s"
imname = sim.image_name

si = {"type": "disk", "flux": 1.0,
      "major": "0.2arcsec",  "minor": "0.1arcsec", "PA": "60.0deg",
      "freq": Nu, "spectral_index": 1.0,
      "direction": "J2000 10h00m00.0s -30d00m00.0s"}

vis = "{0}/{0}.alma.out10.noisy.ms".format(sim.image_name)

if not Path(vis).exists():
    config = '.'.join(sim.array_config.split('.')[:-1])

    print('Generating %s' % sim.image_name)
    cl.done()
    cl.addcomponent(dir=si["direction"], flux=si["flux"], fluxunit='Jy', freq=si["freq"],
                    shape=si["type"],
                    majoraxis=si["major"], minoraxis=si["minor"], positionangle=si["PA"],
                    spectrumtype='spectral index', index=si["spectral_index"])

    skymodel = sim.image_name+'.model'
    ia.fromshape(skymodel, [sim.pixels, sim.pixels, 1, sim.freq_channels], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad = qa.convert(qa.quantity(cell), "rad")['value']
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    rdi, ddi = si["direction"].split()[1:]
    cs.setreferencevalue([qa.convert(rdi, 'rad')['value'],
                          qa.convert(ddi, 'rad')['value']], type="direction")
    cs.setreferencevalue(si["freq"], 'spectral')
    cs.setincrement(ChW, 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    ia.close()

    simobserve(project=sim.image_name, skymodel=skymodel,
               totaltime=totaltime, antennalist=sim.array_config)

model = ['disc']
NuF = float(Nu.split('G')[0])*1.e9
modvars = f"0,0, p[0]*(nu/{sim.center_freq*1.0e9:.4e})**p[1], p[2], p[3], p[4]"
modbounds = [[0.0, None], [-2.0, 2.0], [0.0, None], [0.1, 0.9], [0.0, 180.0]]

size = float(si["major"].split('a')[0])
minor = float(si["minor"].split('a')[0])
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

S = np.array([1.0000, 1.000e+00, 2.000e-01, 5.000e-01, 6.000e+01])

print(f"disc flux at 50GHz (Jy): {r['Parameters'][0]:7.3f} +/- {r['Uncertainties'][0]:.3f}, true: {S[0]:7.3f}")
print(f"disc spectral index:     {r['Parameters'][1]:7.3f} +/- {r['Uncertainties'][1]:.3f}, true: {S[1]:7.3f}")
print(f"disc size (as):          {r['Parameters'][2]:7.3f} +/- {r['Uncertainties'][2]:.3f}, true: {S[2]:7.3f}")
print(f"disc axis ratio:         {r['Parameters'][3]:7.3f} +/- {r['Uncertainties'][3]:.3f}, true: {S[3]:7.3f}")
print(f"disc pos.angle (deg):    {r['Parameters'][4]:7.3f} +/- {r['Uncertainties'][4]:.3f}, true: {S[4]:7.3f}")
