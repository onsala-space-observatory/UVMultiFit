from typing import List

from .modeler import Modeler

class MeasurementSet():

    def __init__(self, vis: str, spw: str = '0', field: int = 0, scans: List = [],
                 corrected: bool = False, uniform: bool = False, uvtaper: float = 0.0,
                 chanwidth: int = 1, timewidth: int = 1, stokes: str = 'I',
                 MJDrange: List[float] = [-1.0, -1.0], ldfac: float = 1.22,
                 phase_center: str = '', pbeam: bool = False, wgt_power: float = 1.0,
                 dish_diameter: float = 0.0) -> None:

        self.vis = vis
        self.pbeam = pbeam
        self.dish_diameter = dish_diameter

        self.spwlist = []                        # type: List
        self.antnames = []                       # type: List
        self.userDiameter = []                   # type: List
        self.tArr = []                           # type: List
        self.averfreqs = []                      # type: List
        self.refpos = [0.0, 0.0]                 # type: List[float]
        self.Nspw = 0                            # type: int
        ...

    def check_measurementset(self) -> bool:
        ...

    def read_data(self, takeModel: bool = False) -> bool:
        ...

    def write_model(self, column: str, model: Modeler) -> bool:
        ...
