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

sim = Simulation(array_config="alma.out10.cfg", center_freq=150.0, pixels=1000, cell_size=0.01,
                 freq_channels=100, channel_width=0.01, image_name="Discs_Gaussian", total_time=300)

image_size = sim.pixels * sim.cell_size
Nu = f"{sim.center_freq}GHz"

cell = f"{sim.cell_size}arcsec"
ChW = f"{sim.channel_width}GHz"
totaltime = f"{sim.total_time}s"
diameter = 1.0
diameter_delta = 0.2
sigma = 0.25

s = [{"type": "disk",
      "flux": 1.0,
      "major": "0.5arcsec",
      "minor": "0.5arcsec",
      "PA": "0.0deg",
      "freq": Nu,
      "spectral_index": 0.0,
      "direction": "J2000 10h00m00.0s -30d00m00.0s"},
     {"type": "disk",
      "flux": -0.25,
      "major": "0.25arcsec",
      "minor": "0.25arcsec",
      "PA": "0.0deg",
      "freq": Nu,
      "spectral_index": 0.0,
      "direction": "J2000 10h00m00.0s -30d00m00.0s"},
     {"type": "gaussian",
      "flux": 0.5,
      "major": "0.4arcsec",
      "minor": "0.4arcsec",
      "PA": "0.0deg",
      "freq": Nu,
      "spectral_index": 0.0,
      "direction": "J2000 10h00m00.0s -30d00m02.0s"}]

imname = sim.image_name
skymodel = sim.image_name+'.model'
vis = "{0}/{0}.alma.out10.noisy.ms".format(sim.image_name)

if not Path(vis).exists():
    cl.done()
    config = '.'.join(sim.array_config.split('.')[:-1])

    print('Generating %s' % sim.image_name)

    for si in s:
        cl.addcomponent(dir=si["direction"], flux=si["flux"], fluxunit="Jy", freq=si["freq"], shape=si["type"],
                        majoraxis=si["major"], minoraxis=si["minor"], positionangle=si["PA"],
                        spectrumtype="spectral index", index=si["spectral_index"])

    ia.fromshape(skymodel, [sim.pixels, sim.pixels, 1, sim.freq_channels], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad = qa.convert(qa.quantity(cell), "rad")['value']
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    rdi, ddi = s[0]["direction"].split()[1:]
    cs.setreferencevalue([qa.convert(rdi, 'rad')['value'],
                          qa.convert(ddi, 'rad')['value']], type="direction")
    cs.setreferencevalue(s[0]["freq"], 'spectral')
    cs.setincrement(ChW, 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    ia.close()

    simobserve(project=imname, skymodel=skymodel,
               totaltime='300s', antennalist=sim.array_config)

doff = -2.  # Dec offset of Gaussian (arcsec.)
model = [sim.image_name]

Gvar = "0.,p[3],p[4],p[5],1.,0."
D1var = "0.,0.,p[0],p[1],1.,0."
D2var = "0.,0.,-p[0]*p[2]**2.,p[1]*p[2],1.,0."

shapes = [str(si["type"]) for si in s]
fluxes = [float(si["flux"]) for si in s]
sizes = [float(si["major"].split('a')[0]) for si in s]
minors = [float(si["minor"].split('a')[0]) for si in s]
PAs = [float(si["PA"].split('d')[0]) for si in s]

Pars = [fluxes[0], sizes[0], sizes[1]/sizes[0], doff, fluxes[2], sizes[2]]
BiasFac = [0.7, 1.3, 0.8, 0.75, 1.2, 1.3]

# For some reason, the disc flux in CASA is wrong by a factor ~0.2477 (??)
# You can check this, e.g. on the model image with the viewer.
Trues = Pars
Trues[0] *= 0.2477

S = Pars
# S = [fluxes[0], sizes[0], sizes[1]/sizes[0], doff, fluxes[2], sizes[2]]
# S = [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e] % tuple(Pars)

# visname = '%s' % ('%s/%s.%s.noisy.ms' % (imname, imname, config))
for idx, shape in enumerate(shapes):
    if shape == 'disk':
        shapes[idx] = 'disc'
    if shape == 'gaussian':
        shapes[idx] = 'Gaussian'

# modelshape = [%s]' % (','.join(["'" + shi + "'" for shi in shapes]))


NuF = float(Nu.split('G')[0])*1.e9
modvars = [D1var, D2var, Gvar]
initial = [Pars[pi]*BiasFac[pi] for pi in range(len(Pars))]
modbounds = [[0.0, None], [0.0, None], [0.0, None], [-15, 15], [0.0, None], [0.0, None]]

pbeam = False
r = uvm.uvmultifit(vis=vis, spw='0',
                   model=shapes, OneFitPerChannel=False,
                   var=modvars, write='residuals',
                   p_ini=initial, pbeam=pbeam,
                   bounds=modbounds)

print("TEST 3: MULTI-COMPONENT WITH CORRELATED VARIABLES\n")
for pi in range(len(initial)):
    print("Parameter %i: %.4f +/- %.4f | True value %.2f" % (
        pi, r['Parameters'][pi], r['Uncertainties'][pi], S[pi]))

default = "modelfit.dat"
if Path(default).exists():
    os.rename(default, os.path.basename(__file__).replace(".py", ".out"))
