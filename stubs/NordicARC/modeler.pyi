from typing import List, Dict

from .measurementset import MeasurementSet

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
        self.OneFitPerChannel = OneFitPerChannel
        self.model = model
        self.var = var
        self.fixed = fixed
        self.takeModel: bool = 'model_column' in self.fixed
        self.fixed = fixed
        self.fixedvar = fixedvar
        self.scalefix = scalefix
        self.p_ini = p_ini
        self.bounds = bounds
        self.bounds = bounds
        self.p_ini = p_ini

        ...

    def init_data(self, nspw: int, ms: MeasurementSet) -> bool:
        ...

    def init_model(self, nspw: int, tArr: List, averfreqs: List, refposnspw: List[float]) -> bool:
        ...

    def fit(self, ms: MeasurementSet, write_model: int, cov_return: bool, redo_fixed: bool = True, reset_flags: bool = False):
        ...
