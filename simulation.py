from dataclasses import dataclass

@dataclass
class Simulation:
    """Class for describing a UVMultiFit simulation."""
    array_config: str
    center_freq: float
    pixels: int
    cell_size: float
    freq_channels: int
    channel_width: float
    image_name: str
    total_time: int

if __name__ == "__main__":
    sim = Simulation("alma.out10.cfg", center_freq=50.0, pixels=1000, cell_size=0.01,
                     freq_channels=100, channel_width=0.01, image_name="Disc", total_time=300)
    arrayconfig = sim.array_config
    Nu = f"{sim.center_freq}GHz"
    Npix = sim.pixels
    cellN = sim.cell_size
    cell = f"{cellN}arcsec"
    Nnu = sim.freq_channels
    ChW = f"{sim.channel_width}GHz"
    imname = sim.image_name
    totaltime = f"{sim.total_time}s"
    si = ['disk', 1.0, '0.2arcsec', '0.1arcsec', '60.0deg',
          Nu, 1.0, "J2000 10h00m00.0s -30d00m00.0s"]
    modvars = f"0,0, p[0]*(nu/{sim.center_freq*1.0e9:.4e})**p[1], p[2], p[3], p[4]"

    config = '.'.join(arrayconfig.split('.')[:-1])
    print(si)
    print(config)
    print(modvars)
    print('Generating %s' % imname)
