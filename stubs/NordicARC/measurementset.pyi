from typing import List

class MeasurementSet():

    def __init__(self, vis: str, spw: str = '0', field: int = 0, scans: List = [],
                 corrected: bool = False, uniform: bool = False, uvtaper: float = 0.0,
                 chanwidth: int = 1, timewidth: int = 1, stokes: str = 'I',
                 MJDrange: List[float] = [-1.0, -1.0], ldfac: float = 1.22,
                 phase_center: str = '', pbeam: bool = False, wgt_power: float = 1.0,
                 dish_diameter: float = 0.0) -> None:
        ...
