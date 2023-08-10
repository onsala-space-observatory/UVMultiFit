import sys
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

sim = Simulation(array_config="alma.out10.cfg", center_freq=50.0, pixels=1000, cell_size=0.01,
                 freq_channels=100, channel_width=0.01, image_name="GaussianRing", total_time=300)

image_size = sim.pixels * sim.cell_size
Nu = f"{sim.center_freq}GHz"
cell = f"{sim.cell_size}arcsec"
ChW = f"{sim.channel_width}GHz"
totaltime = f"{sim.total_time}s"
imname = sim.image_name
diameter = 1.0
diameter_delta = 0.2
sigma = 0.25

si = {"type": "disk", "flux": 1.0,
      "major": "0.2arcsec",  "minor": "0.1arcsec", "PA": "60.0deg",
      "freq": Nu, "spectral_index": 1.0,
      "direction": "J2000 10h00m00.0s -30d00m00.0s"}

skymodel = sim.image_name+'.model'
vis = "{0}/{0}.alma.out10.noisy.ms".format(sim.image_name)
print(vis)

if not Path(skymodel).exists():
    if Path("Disc/Disc.alma.out10.noisy.ms").exists():
        os.system('cp -r %s.model  %s.model' % ("Disc", sim.image_name))
        ia.open(skymodel)

        X = np.linspace(-image_size/2.0, image_size/2.0, sim.pixels)
        XX = np.power(np.outer(X, np.ones(sim.pixels)), 2.0)
        RSq = XX + np.transpose(XX)
        R = np.sqrt(RSq)
        data = ia.getchunk()
        for i in range(sim.freq_channels):
            source_radius = (diameter + diameter_delta*i/sim.freq_channels)/2.0
            mask = R >= source_radius
            data[:, :, 0, i] = np.exp(-np.power(R-source_radius, 2.0)/(2.0*sigma**2.0))
            data[:, :, 0, i] /= np.sum(data[:, :, 0, i])

        ia.putchunk(data)
        ia.close()
    else:
        print("Please run the Disc model first! (python3 disc_profile.py)")
        sys.exit(1)

if not Path(vis).exists():
    config = '.'.join(sim.array_config.split('.')[:-1])

    print('Generating %s' % sim.image_name)

    simobserve(project=sim.image_name, skymodel=skymodel,
               totaltime=totaltime, antennalist=sim.array_config)

model = [sim.image_name]
NuF = float(Nu.split('G')[0])*1.e9
modvars = "p[3],p[4],p[0],p[1],1.0,0.0,p[2]"
modbounds = [[0.0, None], [0.0, diameter*1.5], [sigma*0.8, sigma*1.2], [-1.0, 1.0], [-1.0, 1.0]]

size = float(si["major"].split('a')[0])
minor = float(si["minor"].split('a')[0])
initial = [0.8, diameter*0.8, sigma*1.2, 0.0, 0.0]

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
    ims = sub.imshow(np.transpose(resdat), origin='lower',
                     extent=[image_size, -image_size, -image_size, image_size])
    sub.set_xlabel('RA offset (arcsec)')
    sub.set_ylabel('Dec offset (arcsec)')
    cb = plt.colorbar(ims)
    cb.set_label('Jy/beam')
    # plt.savefig('%s.png' % Cfile)
    plt.show()
    impeak = np.max(resdat)

pbeam = False
r = uvm.uvmultifit(vis=vis, spw='0',
                   model=model, OneFitPerChannel=True,
                   var=modvars, write='residuals',
                   p_ini=initial, pbeam=pbeam,
                   bounds=modbounds)

print(f"RADIAL PROFILES. SHAPE: {sim.image_name}")
maxdev = 100.*np.max(np.abs(r['Parameters'][0][:, 0] - 1.0)/1.0)
print(f" Model {sim.image_name}. Maximum Flux deviation: {maxdev:.4f} %")
ModD = diameter + np.linspace(0, diameter_delta, np.shape(r['Frequency'][0])[0])
maxdev = 100.*np.max(np.abs(r['Parameters'][0][:, 1] - ModD)/ModD)
print(f" Model {sim.image_name}. Maximum Size deviation: {maxdev:.4f} %")

default = "modelfit.dat"
if Path(default).exists():
    os.rename(default, os.path.basename(__file__).replace(".py", ".out"))
