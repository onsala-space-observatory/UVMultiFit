import sys
import os
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
# from icecream import ic

from casatools import componentlist
from casatools import image
from casatools import quanta
from casatasks import simobserve

from simulation import Simulation

from NordicARC import uvmultifit as uvm

cl = componentlist()
ia = image()
qa = quanta()

sim = Simulation(array_config="vla.d.cfg", center_freq=1.5, pixels=4000, cell_size=0.9,
                 freq_channels=10, channel_width=0.01, image_name="WideField", total_time=1200)

image_size = sim.pixels * sim.cell_size
Nu = f"{sim.center_freq}GHz"

# Gain corruption factor:
Gfac = 1.5

cell = f"{sim.cell_size}arcsec"
ChW = f"{sim.channel_width}GHz"

totaltime = f"{sim.total_time}s"
diameter = 1.0
diameter_delta = 0.2
sigma = 0.25

si = {"type": "delta", "flux": 1.0,
      "major": "0.5arcsec", "minor": "0.5arcsec", "PA": "0.0deg",
      "freq": Nu, "spectral_index": 0.0,
      "direction": "J2000 10h00m00.0s 30d00m02.0s"}

imname = sim.image_name
skymodel = sim.image_name+'.model'
vis = "{0}/{0}.{1}.noisy.ms".format(sim.image_name, os.path.splitext(sim.array_config)[0])

N = 4
incr = 0.20  # in degrees
Fi = np.ones(N*N)
if not Path(vis).exists():
    print('Generating %s' % sim.image_name)

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
    ia.close()

    ia.open(skymodel)
    data = ia.getchunk()
    for i in range(N):
        for j in range(N):
            ipix = int(sim.pixels/2 + (-float(N-1)/2. + i)*incr*3600./sim.cell_size)
            jpix = int(sim.pixels/2 + (-float(N-1)/2. + j)*incr*3600./sim.cell_size)
            print(f"pixel[{i},{j}] = {ipix}, {jpix}")
            data[ipix, jpix, :, :] = Fi[i+j]

    ia.putchunk(data)
    ia.close()

    simobserve(project=imname, skymodel=skymodel, indirection="phref",
               totaltime=totaltime, antennalist=sim.array_config)
    print("Done")

model = ['delta' for i in range(N*N)]
RAs = [incr*3600.*((N-1.)/2.-i) for i in range(N)]
Decs = [incr*3600.*(-(N-1.)/2.+i) for i in range(N)]
phref = si["direction"]

modvars = []
initial = []
bounds = []
S = []
i = 0
for ra in RAs:
    for dec in Decs:
        modvars.append('%8.3f+p[%i], %8.3f+p[%i], p[%i]' % (ra, 3*i, dec, 3*i+1, 3*i+2))
        initial += [0., 0., 0.8]
        bounds += [[-5.,5], [-5.,5.], [0.,None]]
        S += [0., 0., 1.]
        i += 1

# ic(modvars)
# ic(initial)
# ic(phref)
# ic(bounds)

NuF = float(Nu.split('G')[0])*1.e9

field = ''
r = uvm.uvmultifit(vis=vis, model=model, var=modvars, field=field,
                   p_ini=initial, NCPU=4, column='data', write='residuals',
                   OneFitPerChannel=False, phase_center=phref, pbeam=True,
                   dish_diameter=25.0, ldfac=1.0, LMtune=[1.e-3, 10., 1.e-5, 10], bounds=bounds)

print("TEST 4: WIDEFIELD\n")
for pi in range(len(initial)):
    print("Parameter %i: %.4f +/- %.4f | True value %.2f" % (
        pi, r['Parameters'][pi], r['Uncertainties'][pi], S[pi]))

default = "modelfit.dat"
if Path(default).exists():
    os.rename(default, os.path.basename(__file__).replace(".py", ".out"))
