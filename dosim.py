from casatools import componentlist
from casatools import image
from casatools import quanta
from casatasks import simobserve

from simulation import Simulation

cl = componentlist()
ia = image()
qa = quanta()

sim = Simulation(array_config="alma.out10.cfg", center_freq=50.0, pixels=1000, cell_size=0.01,
                 freq_channels=100, channel_width=0.01, image_name="Disc", total_time=300)

Nu = f"{sim.center_freq}GHz"
cell = f"{sim.cell_size}arcsec"
ChW = f"{sim.channel_width}GHz"
totaltime = f"{sim.total_time}s"
imname = sim.image_name

si = {"type": "disk", "flux": 1.0,
      "major": "0.2arcsec",  "minor": "0.1arcsec", "PA": "60.0deg",
      "freq": Nu, "spectral_index": 1.0,
      "direction": "J2000 10h00m00.0s -30d00m00.0s"}

# si = ['disk', 1.0, '0.2arcsec', '0.1arcsec', '60.0deg', Nu, 1.0, "J2000 10h00m00.0s -30d00m00.0s"]

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
