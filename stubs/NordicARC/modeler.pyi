from typing import List, Dict

class Modeler():
    def __init__(self, model: List[str] = ['delta'], var: List[str] = ['p[0], p[1], p[2]'],
                 p_ini: List[float] = [0.0, 0.0, 1.0], bounds: List = None,
                 OneFitPerChannel: bool = False, only_flux: bool = False,
                 fixed: List = [], fixedvar: List = [], scalefix: str = '1.0',
                 phase_gains: Dict = {}, amp_gains: Dict = {},
                 method: str = "levenberg", HankelOrder: int = 80,
                 LMtune: List[float] = [1.e-3, 10., 1.e-5, 200, 1.e-3],
                 SMPtune: List[float] = [1.e-4, 1.e-1, 200],
                 NCPU: int = 4, proper_motion: float = 0.0) -> None:
        ...
