import re
import time
import logging

from typing import List, Dict

import numpy as np
from scipy import special        # type: ignore

import uvmultimodel as uvmod     # type: ignore

from .utils import get_list_of_strings, check_proper_motion
from .simplex import _mod_simplex


class Model:

    def __init__(self, name, index, nparams):
        self.name = name
        self.index = index
        self.nparams = nparams

    def __repr__(self):
        txt = "Model("
        txt += f"name = {self.name}, "
        txt += f"index = {self.index}, "
        txt += f"nparams = {self.nparams})"
        return txt


class Modeler():
    """ Class to deal with model equations and fitting procedures.

    It convert strings representing models and parameters into compiled equations, to be used in a ChiSq
    visibility fitting. It also interacts with the C++ extension of UVMultiFit."""

    LIGHTSPEED = 2.99792458e+8
    FOURIERFACTOR = (2.0 * np.pi) * np.pi / 180.0 / 3600.0
    logger = logging.getLogger("modeler")
    isNumerical = ['GaussianRing']

    implemented_models = [Model('delta',        index=0, nparams=3),
                          Model('Gaussian',     index=1, nparams=6),
                          Model('disc',         index=2, nparams=6),
                          Model('ring',         index=3, nparams=6),
                          Model('sphere',       index=4, nparams=6),
                          Model('bubble',       index=5, nparams=6),
                          Model('expo',         index=6, nparams=6),
                          Model('power-2',      index=7, nparams=6),
                          Model('power-3',      index=8, nparams=6),
                          Model('GaussianRing', index=9, nparams=7)]

    def __init__(self, model: List[str] = ['delta'], var: List[str] = ['p[0], p[1], p[2]'],
                 p_ini: List[float] = [0.0, 0.0, 1.0], bounds: List = None,
                 OneFitPerChannel: bool = False, only_flux: bool = False,
                 fixed: List = [], fixedvar: List = [], scalefix: str = '1.0',
                 phase_gains: Dict = {}, amp_gains: Dict = {},
                 method: str = "levenberg", HankelOrder: int = 80,
                 LMtune: List[float] = [1.e-3, 10., 1.e-5, 200, 1.e-3],
                 SMPtune: List[float] = [1.e-4, 1.e-1, 200],
                 NCPU: int = 4, proper_motion: float = 0.0) -> None:
        """ Just the constructor of the 'modeler' class."""
        self.logger.debug("modeler::__init__")
        self.model = get_list_of_strings(model)

        self.propRA = np.zeros(len(model), dtype=np.float64)
        self.propDec = np.zeros(len(model), dtype=np.float64)
        self.fixed = fixed
        if isinstance(fixedvar, list) and all(isinstance(elem, (int, float)) for elem in fixedvar):
            self.fixedvar = fixedvar
        else:
            self.logger.error("fixed variables must be a list of numerical values (ints or floats)")
            self.fixedvar = []

        self.var = get_list_of_strings(var)  # variables of the model components.
        self.scalefix = scalefix
        self._hankel_order = HankelOrder if self.is_numerical_model(self.model) else 0
        self.p_ini = p_ini
        self.bounds = self.check_bounds(bounds, p_ini)
        self._spectral_mode = OneFitPerChannel

        # Lists of compiled functions (one per model component).
        # These will take p and nu and return the variables of the model:

        self.varfunc = [0.0 for component in model]
        self.proper_motion = check_proper_motion(proper_motion, self.model)
        self._NCPU = NCPU
        self.t0 = 0.0
        self.t1 = 1.e12
        self._only_flux = False
        if isinstance(only_flux, bool):
            self._only_flux = only_flux
        else:
            self.logger.error("'only_flux' must be boolean")

        self.Nants = 0
        self.initiated = False
        self.addfixed = False
        self.expka = 2. * np.log(2.)
        self.pow2ka = 1. / (np.pi * np.log(2.))
        self.pow3ka = np.sqrt(2.**(2. / 3.) - 1.)
        self.FouFac = (2. * np.pi) * np.pi / 180. / 3600.
        self.LtSpeed = 2.99792458e+8
        self.deg2rad = np.pi / 180.
        self.calls = 0   # Total number of calls during the fit
        self.minnum = np.sqrt(np.finfo(np.float64).eps)
        # Tells us if the fixedmodel array exists and should be used.
        self.removeFixed = False
        self.varfixed = []     # type: List
        self.dpar = []         # type: List
        self.method = method
        self.Hessian = np.zeros((1, 1), dtype=np.float64)
        self.Gradient = []         # type: List
        self.iscancoords = []      # type: List
        self.varbuffer = []        # type: List
        self.isGain = []           # type: List
        self.fittablebool = []     # type: List
        self.fittable = []         # type: List
        self.wgtcorr = []          # type: List
        self.dt = []               # type: List
        self.dtArr = []            # type: List
        self.dtIdx = []            # type: List
        self.par2 = []             # type: List
        # Arrays of data and pointers (to share access with C library):
        self.data = []             # type: List
        self.wgt = []              # type: List
        self.uv = []               # type: List
        self.offset = []           # type: List
        self.output = []           # type: List
        self.freqs = []            # type: List
        self.fixedmodel = []       # type: List
        self.ants = []             # type: List
        # Model indices (to let the C++ library know which model is which component):
        self.imod = []             # type: List
        self.ifixmod = []          # type: List
        # spw and channel to fit (-1 means fit to the continuum):
        self.currspw = 0
        self.currchan = 0
        # C++ library:
        self.Ccompmodel = uvmod.modelcomp
        # Parameters computed in unbound space:
        self.par2 = []
        # Levenberg-Marquardt parameters:
        self.LMtune = LMtune
        self.SMPtune = SMPtune
        # Gains:
        self.GainBuffer = []       # type: List
        #   self.phaseGainBuffer = []
        self.parDependence = []    # type: List
        self.strucvar = []         # type: List
        KGaus = np.sqrt(1. / (4. * np.log(2.)))
        # Some useful functions:
        self.compiledScaleFixed = lambda p, nu: 1.0 + 0.0
        self.LorentzLine = lambda nu, nu0, P, G: P * 0.25 * G * G / (np.power(nu - nu0, 2.) + (0.5 * G)**2.)
        self.GaussLine = lambda nu, nu0, P, G: P * np.exp(-np.power((nu - nu0) / (G * KGaus), 2.))

        self.pieceWise = lambda t, p0, p1, t0, t1: np.clip(p0 + (p1 - p0) * (t - t0) / (t1 - t0), p0, p1)

        self.wgtEquation = lambda D, Kf: -D * Kf

        # List of currently-supported model components:
        self.phase_gains = phase_gains
        self.amp_gains = amp_gains

        mixed_phase, self.phase_gainsNu, self.phase_gainsT, self.phase_gainsNuT = self.check_gains(self.phase_gains)
        mixed_amp, self.amp_gainsNu, self.amp_gainsT, self.amp_gainsNuT = self.check_gains(self.amp_gains)

        self.gainFunction = [self.phase_gainsNuT, self.amp_gainsNuT, self.phase_gainsNu,
                             self.amp_gainsNu, self.phase_gainsT, self.amp_gainsT]
        self.check_model_consistency()
        if mixed_phase != mixed_amp:
            self.logger.error("inconsistent 'amp_gains' and 'phase_gains'")
        self.useGains = list(self.phase_gainsNuT.keys()) + list(self.amp_gainsNuT.keys()) \
            + list(self.phase_gainsNu.keys()) + list(self.amp_gainsNu.keys()) \
            + list(self.phase_gainsT.keys()) + list(self.amp_gainsT.keys())

        self.isMixed = mixed_phase

    def dump(self):
        temp = vars(self)
        print(f"{self.__class__.__name__}(")
        for item in sorted(temp.keys()):
            print(f"  {item}: {temp[item]}")
        print(")")

    def __repr__(self):
        txt = "Modeler("
        txt += f"model = {self.model}, "
        txt += f"var = {self.var}, "
        txt += f"fixed = {self.fixed}, "
        txt += f"fixedvar = {self.fixedvar}, "
        txt += f"scalefix = {self.scalefix}, "
        txt += f"NCPU = {self._NCPU}, "
        txt += f"flux_only = {self._only_flux}, "
        txt += f"hankel_order = {self._hankel_order}, "
        txt += f"phase_gains = {self.phase_gains}, "
        txt += f"phase_gains = {self.phase_gains})"
        return txt

    @classmethod
    def is_numerical_model(cls, model):
        """Test if the model is numeric.

        Checks whether any of the numerical models (currently only "GaussianRing") is used by the actual model.

        Args:
            model (str): the actual model

        Returns:
            bool: True if model contains a numerical model.
        """
        cls.logger.debug("Modeler::is_numerical_model")
        for mods in cls.isNumerical:
            if mods in model:
                return True
        return False

    @classmethod
    def check_gains(cls, gains):
        cls.logger.debug("Modeler::check_gains")
        if not isinstance(gains, dict):
            cls.logger.error("gains must be a dictionary!")
        mixed_gains = any(isinstance(key, int) for key in gains.keys())
        gainsNu = {}
        gainsT = {}
        gainsNuT = {}
        for key in gains.keys():
            if key in ['nuG', 'tG']:
                if mixed_gains:
                    cls.logger.error("you cannot define split gains and mixed gains in the same fit")
                if key == 'nuG':
                    gainsNu = gains['nuG']
                elif key == 'tG':
                    gainsT = gains['tG']
            elif not isinstance(key, int):
                cls.logger.error("the keys of 'phase_gains' must be integers or 'nuG/tG'!")
            else:
                gainsNuT[key] = gains[key]
        return (mixed_gains, gainsNu, gainsT, gainsNuT)

    @property
    def NCPU(self):
        return self._NCPU

    @NCPU.setter
    def NCPU(self, no_of_cpus):
        self._NCPU = no_of_cpus

    @property
    def flux_only(self):
        return self._only_flux

    @flux_only.setter
    def flux_only(self, flag):
        self._only_flux = flag
        if flag and len(self.p_ini) != len(self.model):
            self.logger.error(f"only_flux=True, but number of parameters ({len(self.p_ini)}) "
                              f"!= number of model components ({len(self.model)}).")

    @property
    def hankel_order(self):
        return self._hankel_order

    @hankel_order.setter
    def hankel_order(self, order):
        self._hankel_order = order

    @property
    def spectral_mode(self):
        return self._spectral_mode

    @spectral_mode.setter
    def spectral_mode(self, flag):
        self._spectral_mode = flag

    @classmethod
    def get_parameter_indices(cls, par):
        cls.logger.debug("Modeler::get_parameter_indices")
        indices = []
        paridx = zip([m.start() + 1 for m in re.finditer(r'\[', par)],
                     [m.start() for m in re.finditer(r'\]', par)])
        for m, n in paridx:
            indices.append(int(par[m:n]))
        return indices

    @classmethod
    def check_bounds(cls, bounds, p_ini):
        """Set (and check) the bounds in the fit"""
        cls.logger.debug("Modeler::check_bounds")
        if not bounds:
            bounds = [[None, None] for p in p_ini]
        elif len(bounds) != len(p_ini):
            cls.logger.error(f"length of 'bounds' ({len(bounds)}) "
                             f"is not equal to number of inital values ({len(p_ini)})")
            bounds = [[None, None] for p in p_ini]

        for b, bound in enumerate(bounds):
            if not bound:
                bounds[b] = [None, None]
                continue
            if len(bounds[b]) != 2:
                cls.logger.error("bounds should come in pairs of two")
                bounds[b] = [None, None]
                continue
            if bounds[b][0] and p_ini[b] <= bounds[b][0]:
                cls.logger.error(f"initial value ({p_ini[b]:.2e}) of parameter {b} "
                                 f"is smaller (or equal) than its lower boundary ({bounds[b][0]:.2e})")
                bounds[b][0] = None
                continue
            if bounds[b][1] and p_ini[b] >= bounds[b][1]:
                cls.logger.error(f"initial value ({p_ini[b]:.2e}) of parameter {b} "
                                 f"is larger (or equal) than its upper boundary ({bounds[b][1]:.2e})")
                bounds[b][1] = None
        return bounds

    def get_model_index(self, component):
        model_names = [m.name for m in self.implemented_models]
        if component not in model_names:
            self.logger.error(f"model component '{component}' is not known!")
            return -1
        return model_names.index(component)

    def is_valid_model(self, component, nparams):
        """Check that model 'component' is among the implemented models and takes 'nparams' parameters"""

        model_index = self.get_model_index(component)
        if model_index == -1:
            self.logger.error(f"model component '{component}' is not known!")
            return False

        if self.implemented_models[model_index].nparams != nparams:
            self.logger.error(f"wrong number of variables ({nparams}) in '{component}' model")
            return False

        return True

    def check_model_consistency(self):
        """Get the number of parameters and check the model consistency"""
        self.logger.debug("Modeler::check_model_consistency")

        indices = []
        # Typical function that can be used.
        lf = ['GaussLine', 'GaussLine', 'LorentzLine', 'LorentzLine', 'power', 'maximum', 'minimum']
        for i, component in enumerate(self.model):
            vic = self.var[i].count
            checkpars = [p.strip() for p in self.var[i].split(',')]
            nfuncs = sum(list(map(vic, lf)))  # 2*(self.var[i].count('GaussLine')+self.var[i].count('LorentzLine'))
            indices += self.get_parameter_indices(self.var[i])

            if not self.is_valid_model(component, len(checkpars) - nfuncs):
                return False

        # Scalefix must be a string representing a function:
        if not isinstance(self.scalefix, str):
            self.logger.error("'scalefix' should be a string!")
            return False

        indices += self.get_parameter_indices(self.scalefix)
        maxpar = max(indices)
        if len(self.p_ini) != maxpar + 1:
            self.logger.error(f"'p_ini' is of length {len(self.p_ini)}, but {maxpar + 1} parameters are used")

        self.takeModel = 'model_column' in self.fixed
        if self.takeModel:
            self.logger.info("MODEL COLUMN will be taken as fixed model.")
            self.fixed = ['model_column']
            self.fixedvar = []
        else:
            for i, component in enumerate(self.fixed):
                checkpars = self.fixedvar[i].split(',')
                if not self.is_valid_model(component, len(checkpars)):
                    return False
        return indices

    def compile_all_models(self):
        """ Compile all models (fixed, variable, and fixed-scale. Not to be called directly by the user. """
        self.logger.debug("Modeler::compile_all_models")
        ok = self._compile_model()
        ok = ok and self._compile_fixed_model()
        ok = ok and self._compile_scale_factor()
        ok = ok and self._compile_gains()
        return ok

    def execute_code(self, code):
        self.logger.debug("Modeler::execute_code")
        try:
            self.logger.info(f"executing: '{code}'")
            exec(code, locals())
            return True
        except Exception:
            return False

    def _compile_mixed_gain_function(self, gainFunction, which):
        self.logger.debug("Modeler::_compile_mixed_gain_function")
        n = ["phase", "amp"].index(which)
        for ni in gainFunction[n].keys():
            tempstr = gainFunction[n][ni].replace(
                "pieceWise(", "self.pieceWise(t, ").replace(
                    "t", "t[:,__import__('numpy').newaxis]").replace(
                        "nu0", "%.12f" % self.freqs[0][0])
            modstr = "self." + which + "AntsFunc[" + str(ni) + "] = lambda t, nu, p: " + tempstr
            parpos = [x.start() for x in re.finditer(r"p\[", tempstr)]
            for p0 in parpos:
                p1 = tempstr[p0:].find("]") + p0
                self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

            if not self.execute_code(modstr):
                self.logger.error(f"syntax error in {which} gain of antenna {ni}")
                return False
        return True

    def _compile_unmixed_gain_function(self, gainFunction, which):
        self.logger.debug("Modeler::_compile_unmixed_gain_function")
        n = ["phase", "amp"].index(which) + 2
        for ni in gainFunction[n].keys():
            tempstr = gainFunction[n][ni].replace("nu0", "%.12f" % self.freqs[0][0])
            if "t" in tempstr:
                self.logger.error("a frequency-dependent gain cannot depend on time")
                return False

            modstr = "self." + which + "AntsFuncNu[" + str(ni) + "] = lambda nu, p: " + tempstr
            parpos = [x.start() for x in re.finditer(r"p\[", tempstr)]
            for p0 in parpos:
                p1 = tempstr[p0:].find(']') + p0
                self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

            if not self.execute_code(modstr):
                self.logger.error(f"syntax error in {which} gain of antenna {ni}")
                return False
        return True

    def _compile_timedependent_gain_function(self, gainFunction, which):
        self.logger.debug("Modeler::_compile_timedependent_gain_function")
        n = ["phase", "amp"].index(which) + 4
        for ni in gainFunction[n].keys():
            tempstr = gainFunction[n][ni].replace(
                "pieceWise(", "self.pieceWise(t,").replace(
                    "nu0", "%.12f" % self.freqs[0][0])

            if "nu" in tempstr:
                self.logger.error("a time-dependent gain cannot depend on frequency")
                return False

            modstr = "self." + which + "AntsFuncT[" + str(ni) + "] = lambda t, p: " + tempstr
            parpos = [x.start() for x in re.finditer(r"p\[", tempstr)]
            for p0 in parpos:
                p1 = tempstr[p0:].find("]") + p0
                self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

            if not self.execute_code(modstr):
                self.logger.error("syntax error in phase gain of antenna {ni}")
                return False
        return True

    def _compile_gains(self):
        """ Compile the functions related to the antenna gains."""
        self.logger.debug("Modeler::_compile_gains")
        if self.isMixed:
            self.phaseAntsFunc = [lambda t, nu, p: 0.0 for i in range(self.Nants)]
            self.ampAntsFunc = [lambda t, nu, p: 1.0 for i in range(self.Nants)]
        else:
            self.phaseAntsFuncNu = [lambda nu, p: 0.0 for i in range(self.Nants)]
            self.ampAntsFuncNu = [lambda nu, p: 1.0 for i in range(self.Nants)]
            self.phaseAntsFuncT = [lambda t, p: 0.0 for i in range(self.Nants)]
            self.ampAntsFuncT = [lambda t, p: 1.0 for i in range(self.Nants)]

        self.parDependence = [[0] for i in range(self.Nants)]

        if self.isMixed:
            self._compile_mixed_gain_function(self.gainFunction, "phase")
            self._compile_mixed_gain_function(self.gainFunction, "amp")
        else:
            self._compile_unmixed_gain_function(self.gainFunction, "phase")
            self._compile_unmixed_gain_function(self.gainFunction, "amp")
            self._compile_timedependent_gain_function(self.gainFunction, "phase")
            self._compile_timedependent_gain_function(self.gainFunction, "amp")

        for ni in range(self.Nants):
            self.parDependence[ni] = np.unique(self.parDependence[ni]).astype(np.int32)

    def _compile_model(self):
        """ Compiles the variable model, according to the contents of the 'model' and 'parameters' lists."""
        self.logger.debug("Modeler::_compile_model")
        # Define variable model:
        self.varfunc = [0.0 for component in self.model]

        for ii, component in enumerate(self.model):
            tempstr = self.var[ii].replace(
                'LorentzLine(', 'self.LorentLine(nu, ').replace(
                    'GaussLine(', 'self.GaussLine(nu, ').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            if self._only_flux:
                try:
                    tempstr2 = tempstr.split(',')
                    if len(tempstr2) > 3:
                        self.strucvar.append(list(map(float, tempstr2[:2] + tempstr2[3:])))
                    else:
                        self.strucvar.append(list(map(float, tempstr2[:2])))
                except Exception:
                    self.logger.error("if only_flux=True, all variables but the flux must be constants")
                    return False

            modstr = 'self.varfunc[' + str(ii) + '] = lambda p, nu: [' + tempstr + ']'
            if not self.execute_code(modstr):
                self.logger.error(f"syntax error in component number {ii} of the variable model")
                return False

            model_index = self.get_model_index(component)
            if model_index == -1:
                self.logger.error(f"component '{component}' is unknown")
                return False
            self.imod[ii] = model_index
        return True

    def _compile_fixed_model(self):
        """ Compiles the fixed model, according to the contents of the 'fixed' and 'fixedpars' lists."""

        self.logger.debug("Modeler::_compile_fixed_model")
        if len(self.fixed) > 0 and 'model_column' in self.fixed:
            return True

        # Define fixed model:
        self.fixedvarfunc = [0.0 for component in self.fixed]

        self.ifixmod = np.zeros(len(self.fixed), dtype=np.int32)  # []

        for ii, component in enumerate(self.fixed):
            tempstr = self.fixedvar[ii].replace(
                'LorentzLine(', 'self.LorentzLine(nu, ').replace(
                    'GaussLine(', 'self.GaussLine(nu, ').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            modstr = 'self.fixedvarfunc[' + str(ii) + '] = lambda p, nu: [' + tempstr + ']'
            if not self.execute_code(modstr):
                self.logger.error(f"syntax error in component number {ii} of the fixed model")
                return False

            model_index = self.get_model_index(component)
            if model_index == -1:
                self.logger.error(f"component '{component}' is unknown")
                return False
            self.ifixmod[ii] = model_index
        return True

    def _compile_scale_factor(self):
        """ Compiles the scaling factor for the fixed model """

        self.logger.debug("Modeler::_compile_scale_factor")
        tempstr = self.scalefix.replace(
            'LorentzLine(', 'self.LorentzLine(nu, ').replace(
                'GaussLine(', 'self.GaussLine(nu, ').replace(
                    'nu0', '%.12f' % self.freqs[0][0])

        scalefixedstr = 'self.compiledScaleFixed = lambda p, nu: ' + tempstr + ' + 0.0'
        if not self.execute_code(scalefixedstr):
            self.logger.error("syntax error in the flux-scale equation")
            return False
        return True

    ############################################
    #
    #  GET UNBOUND PARAMETER SPACE FROM BOUND PARAMETERS
    #
    def _getPar2(self, mode=0):
        """ Function to change fitting parameters to/from the unbound space from/to the bound space.

        It also comptues the gradient for the Jacobian matrix. Unbound values are in index 0, gradient is
        in index 1, and bound values are in index 2. The equations for the changes of variables are taken
        from the MINUIT package."""

        self.logger.debug("Modeler::_getPar2")
        if self.bounds is None:
            self.par2[1, :] = 1.0
            if mode == 0:
                self.par2[0, :] = self.par2[2, :]
            else:
                self.par2[2, :] = self.par2[0, :]
            return
        if mode == 0:
            par2 = self.par2[2, :]
            for i, p in enumerate(par2):
                if all(b is None for b in self.bounds[i]):
                    self.par2[:, i] = [p, 1.0, p]
                elif self.bounds[i][0] is None and self.bounds[i][1] is not None:
                    self.par2[0, i] = np.sqrt(np.power(self.bounds[i][1] - p + 1, 2.) - 1.)
                elif self.bounds[i][0] is not None and self.bounds[i][1] is None:
                    self.par2[0, i] = np.sqrt(np.power(p - self.bounds[i][0] + 1, 2.) - 1.)
                else:
                    Kbound = (self.bounds[i][1] - self.bounds[i][0]) / 2.
                    self.par2[0, i] = np.arcsin((p - self.bounds[i][0]) / Kbound - 1.0)

        else:
            par2 = self.par2[0, :]
            for i, p in enumerate(par2):
                if all(b is None for b in self.bounds[i]):
                    self.par2[2, i] = p
                elif self.bounds[i][0] is None and self.bounds[i][1] is not None:
                    self.par2[2, i] = self.bounds[i][1] + 1. - np.sqrt(p**2. + 1.)
                elif self.bounds[i][0] is not None and self.bounds[i][1] is None:
                    self.par2[2, i] = self.bounds[i][0] - 1. + np.sqrt(p**2. + 1.)
                else:
                    Kbound = (self.bounds[i][1] - self.bounds[i][0]) / 2.
                    self.par2[2, i] = self.bounds[i][0] + Kbound * (np.sin(p) + 1.)

        for i, p in enumerate(par2):
            if all(b is None for b in self.bounds[i]):
                self.par2[1, i] = 1.0
            elif self.bounds[i][0] is None and self.bounds[i][1] is not None:
                self.par2[1, i] = -self.par2[0, i] / np.sqrt(self.par2[0, i]**2. + 1.)
            elif self.bounds[i][0] is not None and self.bounds[i][1] is None:
                self.par2[1, i] = self.par2[0, i] / np.sqrt(self.par2[0, i]**2. + 1.)
            else:
                Kbound = (self.bounds[i][1] - self.bounds[i][0]) / 2.
                self.par2[1, i] = Kbound * np.cos(self.par2[0, i])

    ############################################
    #
    #  LEVENBERG-MARQUARDT LEAST-SQUARE MINIMIZATION
    #
    def LMMin(self, pars):
        """ Implementation of the Levenberg-Marquardt algorithm.
            Not to be called directly by the user. """

        self.logger.debug("Modeler::LMMin")
        NITER = int(self.LMtune[3] * len(pars))
        self.calls = 0

        Lambda = float(self.LMtune[0])
        kfac = float(self.LMtune[1])
        functol = float(self.LMtune[2])
        if len(self.LMtune) > 4:
            partol = float(self.LMtune[4])
        else:
            partol = functol
        Chi2 = 0
        self.par2[2, :] = pars
        self._getPar2()

        Hessian2 = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        Gradient2 = np.zeros(np.shape(self.Gradient), dtype=np.float64)
        HessianDiag = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupHess = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupGrad = np.zeros(np.shape(self.Gradient), dtype=np.float64)
        Inverse = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupP = np.copy(self.par2[0, :])

        self.Hessian[:, :] = 0.0
        self.Gradient[:] = 0.0

        if self._only_flux and self.currchan == -1:
            nnu = max([len(self.freqs[sp]) for sp in range(len(self.freqs))])
            self.logger.info("computing structure-only parameters")

            for midx, p in enumerate(pars):
                if len(self.strucvar[midx]) > 3:
                    tempvar = self.strucvar[midx][:2] + [p] + self.strucvar[midx][2:]
                else:
                    tempvar = self.strucvar[midx] + [p]
                for i, tv in enumerate(tempvar):
                    for j in range(len(pars) + 1):
                        self.varbuffer[j][midx, i, :nnu] = tv

        CurrChi = self.residuals(self.par2[2, :], -1, dof=False)
        Hessian2[:, :] = self.Hessian * (self.par2[1, :])[np.newaxis, :] * (self.par2[1, :])[:, np.newaxis]
        Gradient2[:] = self.Gradient * self.par2[1, :]
        backupHess[:, :] = Hessian2
        backupGrad[:] = Gradient2

        controlIter = 0
        for i in range(NITER):
            controlIter += 1
            for n in range(len(pars)):
                HessianDiag[n, n] = Hessian2[n, n]
            try:
                goodsol = True
                Inverse[:] = np.linalg.pinv(Hessian2 + Lambda * HessianDiag)
                Dpar = np.dot(Inverse, Gradient2)
                DirDer = sum([Hessian2[n, n] * Dpar[n] * Dpar[n] for n in range(len(pars))])
                DirDer2 = np.sum(Gradient2 * Dpar)
                TheorImpr = DirDer - 2. * DirDer2
            except Exception:
                goodsol = False
                Dpar = 0.0
                TheorImpr = -10.0
            self.par2[0, :] += Dpar
            self._getPar2(mode=1)
            if goodsol:
                self.Hessian[:, :] = 0.0
                self.Gradient[:] = 0.0
                Chi2 = self.residuals(self.par2[2, :], -1, dof=False)
                Hessian2[:, :] = self.Hessian * (self.par2[1, :])[np.newaxis, :] * (self.par2[1, :])[:, np.newaxis]
                Gradient2[:] = self.Gradient * self.par2[1, :]
                RealImpr = Chi2 - CurrChi
            else:
                RealImpr = 1.0

            if TheorImpr != 0.0:
                Ratio = RealImpr / TheorImpr
            else:
                Ratio = 0.0

            if Ratio < 0.25:
                if RealImpr < 0.0:
                    temp = np.sqrt(kfac)
                else:
                    temp = kfac

            elif Ratio > 0.75:
                temp = 1. / np.sqrt(kfac)

            elif not goodsol:
                temp = kfac

            else:
                temp = 1.0

            Lambda *= temp
            if Chi2 == 0.0 or CurrChi == 0.0:
                break

            relchi = Chi2 / CurrChi
            if relchi < 1:
                relchi = 1. / relchi
            todivide = np.copy(backupP)
            todivide[todivide == 0.0] = 1.0
            totest = [np.abs(tt) for tt in self.par2[0, :] / todivide if tt != 0.0]  # totest[totest==0.0] = 1.0
            relpar = max([{True: 1. / pb, False: pb}[pb < 1] for pb in totest])
            if relchi - 1.0 < functol or relpar - 1.0 < partol:
                self.par2[0, :] = backupP
                Hessian2[:, :] = backupHess
                Gradient2[:] = backupGrad
                self._getPar2(mode=1)
                break

            if CurrChi < Chi2:
                self.par2[0, :] = backupP
                Hessian2[:, :] = backupHess
                Gradient2[:] = backupGrad
                self._getPar2(mode=1)
            else:
                CurrChi = Chi2
                backupHess[:, :] = Hessian2
                backupGrad[:] = Gradient2
                backupP[:] = self.par2[0, :]

        if controlIter == NITER:
            self.logger.warning("Reached maximum number of iterations! The algorithm may not have converged!")
            # self.logger.warning("Please, check if the parameter values are meaningful.")

        self._getPar2(mode=1)

        try:
            return True, [self.par2[2, :], np.linalg.pinv(self.Hessian), Chi2]
        except np.linalg.LinAlgError:
            return False, []

    ############################################
    #
    #  NON-ALGEBRAIC MODELS
    #
    # Compute elements of Taylor expansion of the source's Hankel transform:
    def _grid_model(self, imod, tempvar):
        """ Compute elements of Taylor expansion of the source's Hankel transform."""

        self.logger.debug("Modeler::_grid_model")
        n = self._hankel_order - 1

        if imod == 'GaussianRing':   # Gaussian Ring

            a = 1.0
            k = 1. / (2. * np.power(tempvar[6], 2.)) * tempvar[3] * tempvar[3] / 4.
            m = 2 * n + 1

            # Recurrence relation:
            merf = (1. - special.erf(-np.sqrt(k) * a))

            m0 = np.sqrt(np.pi / k) / 2. * merf
            m1 = np.exp(-k * a * a) / 2. / k + np.sqrt(np.pi / k) * a / 2. * merf

            tempvar.append(np.ones(np.shape(m1)))  # term for n=0, normalized

            if m in [0, 1]:
                return tempvar  # Only order n=0.

            res = a * m1 + 1. / (2 * k) * m0  # order m=2
            resaux = np.copy(res)
            res2 = np.copy(m1)
            for mi in range(3, m + 1):
                if isinstance(res, np.float64):
                    res = a * res + (mi - 1) / (2 * k) * res2
                    res2 = resaux
                    resaux = res
                else:
                    res[:] = a * res + (mi - 1) / (2 * k) * res2
                    res2[:] = resaux
                    resaux[:] = res
                if np.mod(mi + 1, 2) == 0:  # (i.e., orders n=1, 2, 3...)
                    tempvar.append(res / m1 * np.power(-1., (mi - 1) / 2)
                                   / np.power(np.math.factorial((mi - 1) / 2), 2.))
            return tempvar
        raise ValueError("Model %i was not correctly interpreted!" % imod)

    def _compute_antenna_gain_corrections(self, pars, spwrange, nui):
        self.logger.debug("Modeler::_compute_antenna_gain_corrections")
        # Compute antenna-gain corrections:
        for spw in spwrange:
            Nchan = len(self.freqs[spw])
            Nint = len(self.dtArr[spw])
            for i in self.useGains:
                for j, pardep in self.parDependence[i]:
                    ptemp = list(pars)   # [pi for pi in pars]
                    if j > 0:
                        j2 = pardep[j] - 1
                        ptemp[j2] += self.dpar[j2]
                    if nui == -1:
                        if self.isMixed:
                            self.GainBuffer[spw][i][j][:] = self.ampAntsFunc[i](self.dtArr[spw],
                                                                                self.freqs[spw], ptemp) * \
                                np.exp(1.j * self.phaseAntsFunc[i](self.dtArr[spw], self.freqs[spw], ptemp))
                        else:
                            self.GainBuffer[spw][i][j][:Nchan] = self.ampAntsFuncNu[i](self.freqs[spw], ptemp) * \
                                np.exp(1.j * self.phaseAntsFuncNu[i](self.freqs[spw], ptemp))
                            self.GainBuffer[spw][i][j][Nchan:Nchan + Nint] = self.ampAntsFuncT[i](
                                self.dtArr[spw], ptemp) * np.exp(1.j * self.phaseAntsFuncT[i](self.dtArr[spw],
                                                                                              ptemp))

                    else:
                        if self.isMixed:
                            self.GainBuffer[spw][i][j][:, nui] = \
                                np.squeeze(self.ampAntsFunc[i](self.dtArr[spw], self.freqs[spw][nui], ptemp)) * \
                                np.exp(1.j * np.squeeze(self.phaseAntsFunc[i](self.dtArr[spw],
                                                                              self.freqs[spw][nui], ptemp)))
                        else:
                            self.GainBuffer[spw][i][j][nui] = self.ampAntsFuncNu[i](self.freqs[spw][nui], ptemp) * \
                                np.exp(1.j * self.phaseAntsFuncNu[i](self.freqs[spw][nui], ptemp))
                            self.GainBuffer[spw][i][j][Nchan:Nchan + Nint] = \
                                self.ampAntsFuncT[i](self.dtArr[spw], ptemp) * \
                                np.exp(1.j * self.phaseAntsFuncT[i](self.dtArr[spw], ptemp))

    def _compute_fixed_model(self, pars, spwrange, nui):
        self.logger.debug("Modeler::_compute_fixed_model")
        self.removeFixed = True
        # isfixed = True
        modbackup = self.imod[0]
        for spw in spwrange:
            self.output[spw][:] = 0.0
            for midx, mi in enumerate(self.ifixmod):
                tempvar = self.fixedvarfunc[midx](pars, self.freqs[spw])
                self.imod[0] = mi
                if self.imod[0] in self.isNumerical:
                    tempvar = self._grid_model(self.imod[0], tempvar)
                for i, tv in enumerate(tempvar):
                    nnu = len(self.freqs[spw])
                    self.varbuffer[0][0, i, :nnu] = tv
                self.Ccompmodel(spw, nui, 0)

        self.imod[0] = modbackup
        return 0

    def residuals(self, pars, mode=-1, dof=True):
        """ Compute the residuals, fixed model, and/or covariance matrix and Chi square.

        This method has a wide and flexible usage, depending on the value of 'mode'

        Parameters
        ----------
        p : `list`
          List of parameter values at which to compute the residuals.
        mode : `int`
          Controls what is computed and returned:
          mode ==  0. Compute the fixed model. Fill-in the output array with it.
          mode == -1. Compute the Hessian matrix and the Error vector. Return the Chi2.
          mode == -2. Only return the Chi2.
          mode == -3. Add the variable model to the output array.
          mode == -4. Write the residuals to the output array.
          mode == -5. Save the calibrated data to the output array.

        The so-called 'output array' is the data that will be saved into the measurement set(s)
        then the "writeModel" method of the parent UVMultiFit instance is called."""

        self.logger.debug("Modeler::residuals")
        # f = open("casa6", "a")

        #  varsize = self.maxNvar + self._hankel_order
        if mode in [0, -3]:
            self.calls = 0
        else:
            self.calls += 1

        # DELTA FOR DERIVATIVES:
        for j, p in enumerate(pars):
            #  self.dpar[j] = self.minnum
            if np.abs(p) < 1.0:
                self.dpar[j] = self.minnum
            else:
                self.dpar[j] = np.abs(p) * self.minnum

        nui, spwrange = self.fit_range(self.currchan, self.currspw, self.freqs)
        if len(self.useGains) > 0 and mode != 0:
            self._compute_antenna_gain_corrections(pars, spwrange, nui)

        if mode == 0:  # Just compute the fixed model and return
            self._compute_fixed_model(pars, spwrange, nui)
            return 0

        currmod = self.model
        currvar = self.varfunc
        # currimod = self.imod
        # isfixed = False

        ChiSq = 0.0
        ndata = 0
        for spw in spwrange:
            _, nnu = np.shape(self.output[spw])
            if nui == -1:
                freqs = self.freqs[spw]
                idx = slice(0, nnu)
            else:
                freqs = self.freqs[spw][nui]
                idx = 0
            scalefx = self.compiledScaleFixed(pars, freqs)
            self.varfixed[0][:nnu] = scalefx
            for j in range(len(pars)):
                ptemp = list(pars)   # [pi for pi in pars]
                ptemp[j] += self.dpar[j]
                # Variables of current component
                scalefx = self.compiledScaleFixed(ptemp, self.freqs[spw])
                self.varfixed[j + 1][idx] = scalefx

            for midx, modi in enumerate(currmod):
                ptemp = [float(pi) for pi in pars]

                if nui == -1 and self._only_flux:  # Continuum
                    currflux = pars[midx]
                    # nstrucpars = len(self.strucvar[midx])
                    self.varbuffer[0][midx, 2, :nnu] = currflux
                else:
                    # Variables of current component
                    tempvar = currvar[midx](pars, freqs)
                    if modi in self.isNumerical:
                        tempvar = self._grid_model(modi, tempvar)

                    for i, tv in enumerate(tempvar):
                        self.varbuffer[0][midx, i, idx] = tv

                if nui == -1 and self._only_flux:
                    ptemp[midx] += self.dpar[midx]
                    self.varbuffer[midx + 1][midx, 2, :nnu] = ptemp[midx]
                else:
                    for j in range(len(pars)):
                        ptemp = list(pars)   # [pi for pi in pars]
                        ptemp[j] += self.dpar[j]
                        # Variables of current component
                        tempvar = currvar[midx](ptemp, freqs)
                        if modi in self.isNumerical:
                            tempvar = self._grid_model(modi, tempvar)
                        for i, tv in enumerate(tempvar):
                            self.varbuffer[j + 1][midx, i, idx] = tv

            ChiSq += self.Ccompmodel(spw, nui, mode)

            if mode == -2:
                if nui == -1:
                    # Only add the unflagged data to compute the DoF
                    ndata += np.sum(self.wgt[spw] > 0.0)
                else:
                    # Only add the unflagged data to compute the DoF
                    ndata += np.sum(self.wgt[spw][:, nui] > 0.0)

        if mode in [-2, -1]:
            if nui < 0:
                self.logger.info(f"iteration #{self.calls}: achieved ChiSq: {ChiSq:.8e}")

        if ChiSq <= 0.0:
            raise ValueError("Invalid Chi Square!"
                             + "Maybe the current fitted value of flux (and/or size) is negative!"
                             + "Please, set BOUNDS to the fit!")

        # print("ChiSq: %.7f" % (ChiSq,), file=f)
        # f.close()
        if mode in [-1, -2, -3]:
            if dof:
                self.calls = 0
                return [ChiSq, ndata]
            return ChiSq
        return None

    ############################################
    #
    #  COMPUTE CHI SQUARE FOR A MODEL REALIZATION (NO DERIVATIVES)
    #
    def chi_square(self, p, bounds=None, p_ini=[]):
        """ Just a wrapper of the 'residuals' function, usually called by simplex."""

        self.logger.debug("Modeler::chi_square")
        inside = True
        if bounds is not None:
            for i, bound in enumerate(bounds):
                if isinstance(bound, list):
                    vmin = bound[0] is not None and p[i] <= bound[0]
                    vmax = bound[1] is not None and p[i] >= bound[1]
                    if vmin or vmax:
                        inside = False

        # If trying to explore outside the boundaries, set the function to be
        # equal to that at p_ini (i.e., likely to be the worst from all the
        # candidate function minima):
        if inside:
            p_comp = list(p)
        else:
            p_comp = list(p_ini)

        return self.residuals(p_comp, mode=-2, dof=False)

    def init_data(self, nspw, ms):
        """ Initiates the data pointers of the 'modeler' instance.

        The 'modeler' instance stores pointers to the data and metadata, the compiled model (and
        fixedmodel), the parameter values, and all the methods to compute residuals, Chi Squared, etc."""

        self.logger.debug("Modeler::init_data")
        # Reset pointers for the modeler:
        # self.deleteData()

        # Set number of spectral windows:
        gooduvm = uvmod.setNspw(int(nspw))

        if gooduvm != nspw:
            self.logger.error("Error in the C++ extension!")
            return False

        # Maximum number of frequency channels (i.e., the maximum from all the selected spws):
        self.maxnfreq = 0

        # Fill in data pointers for each spw:
        for spidx in range(nspw):
            self.data.append([])
            self.dt.append([])
            self.dtArr.append([])
            self.dtIdx.append([])
            self.wgt.append([])
            self.wgtcorr.append([])
            self.uv.append([])
            self.offset.append([])
            self.output.append([])
            self.fixedmodel.append([])
            self.ants.append([])

            self.fittable.append([])
            self.fittablebool.append([])
            self.isGain.append([])

            self.iscancoords = self.iscancoords

            self.freqs.append(np.require(ms.averfreqs[spidx], requirements=['C', 'A']))
            self.maxnfreq = max(self.maxnfreq, len(ms.averfreqs[spidx]))

            # Only have to multiply by freq, to convert these into lambda units:
            ulambda = np.require(self.FOURIERFACTOR * ms.u[spidx] / self.LIGHTSPEED, requirements=['C', 'A'])
            vlambda = np.require(self.FOURIERFACTOR * ms.v[spidx] / self.LIGHTSPEED, requirements=['C', 'A'])
            wlambda = np.require(self.FOURIERFACTOR * ms.w[spidx] / self.LIGHTSPEED, requirements=['C', 'A'])

            # Data, uv coordinates, and weights:

            # Data are saved in two float arrays (i.e., for real and imag):
            self.data[-1] = ms.averdata[spidx]  # [0], ms.averdata[spidx][1]]
            self.wgt[-1] = ms.averweights[spidx]
            self.uv[-1] = list([ulambda, vlambda, wlambda])
            self.offset[-1] = list([ms.RAshift[spidx], ms.Decshift[spidx], ms.Stretch[spidx]])
            self.ants[-1] = list([ms.ant1[spidx], ms.ant2[spidx]])
            PBFactor = -np.sqrt(ms.KfacWgt[ms.ant1[spidx]] * ms.KfacWgt[ms.ant2[spidx]])

            # Set number of antennas and whether each one has a fittable gain:
            self.Nants = self.Nants
            self.isGain[-1] = np.require(np.zeros(len(ms.t[spidx]), dtype=np.int8), requirements=['C', 'A'])
            for i in self.useGains:
                mask0 = self.ants[-1][0] == i
                mask1 = self.ants[-1][1] == i
                self.isGain[-1][np.logical_or(mask0, mask1)] = True

            # Release memory:
            del ulambda, vlambda, wlambda

            # Time spent on observation:
            self.dt[-1] = np.require(ms.t[spidx] - self.t0, requirements=['C', 'A'])
            self.dtArr[-1] = np.require(ms.tArr[spidx] - self.t0, requirements=['C', 'A'])
            self.dtIdx[-1] = np.require(ms.tIdx[spidx], requirements=['C', 'A'])

            # Array to save the residuals (or model, or any output from the C++ library):
            self.output[-1] = np.require(np.zeros(np.shape(ms.averdata[spidx]),
                                                  dtype=np.complex128), requirements=['C', 'A'])
            if self.takeModel:
                try:
                    self.output[-1][:] = ms.avermod[spidx]
                except Exception:
                    self.logger.error("You already used the model column! Should read it again!")

            ########
            # Array of booleans, to determine if a datum enters the fit:
            self.fittable[-1] = self.wgt[-1][:, 0] > 0.0
            # self.wgtcorr[-1] = np.require(np.copy(np.tile(PBFactor,(len(self.model), 1))),
            #                                       requirements=['C', 'A'])
            self.wgtcorr[-1] = np.require(np.copy(PBFactor), requirements=['C', 'A'])

            del PBFactor

            self.fittablebool[-1] = np.require(np.copy(self.fittable[-1]).astype(bool),
                                               requirements=['C', 'A'])

            gooduvm = uvmod.setData(spidx, self.uv[-1][0], self.uv[-1][1], self.uv[-1][2],
                                    self.wgt[-1], self.data[-1],
                                    self.output[-1], self.freqs[-1],
                                    self.fittable[-1], self.wgtcorr[-1], self.dt[-1],
                                    self.dtArr[-1], self.dtIdx[-1],
                                    self.offset[-1][0], self.offset[-1][1], self.offset[-1][2],
                                    self.ants[-1][0], self.ants[-1][1],
                                    self.isGain[-1], self.Nants)

            # if not gooduvm:
            #     self.logger.error("Error in the C++ extension!")
            #     return False
            #
            # try:
            #     for spidx in range(nspw - 1, -1, -1):
            #         del self.avermod[spidx]
            # except Exception:
            #     pass
            # gc.collect()
        # self.logger.debug("leaving init_data")
        return True

    def init_model(self, nspw, tArr, averfreqs, refpos):
        """ Allocates memory for the modeler data, which will be used by the C++ extension.

        Also compiles the model variables. It is a good idea to run this method everytime that
        the user reads new data to be fitted (i.e., after a call to readData())

        It is indeed MANDATORY to run this function if the model to be fitted has changed."""

        self.logger.debug("Modeler::init_model")
        if not self.initiated:
            # self.clearPointers(1)
            # self.deleteModel()
            # else:
            self.logger.info("UVMultiFit compiled model(s) do not seem to exist yet.")

        # Set number of threads:
        ncpu = int(self.NCPU)
        gooduvm = uvmod.setNCPU(ncpu)
        if gooduvm != ncpu:
            self.logger.error("Error in the C++ extension!")
            return False

        # Fill-in proper motions (in as/day)
        for i in range(len(self.model)):
            self.propRA[i] = self.proper_motion[i][0] / 365.
            self.propDec[i] = self.proper_motion[i][1] / 365.

        #####
        # Array to save the variables of the model and the 'scalefix' value,
        # all of them as a function of frequency. The C++ library will read the variables from here
        # at each iteration and for each model component:
        self.initiated = True
        self.logger.debug("model initiated")
        self.logger.debug(f"length of model: {str(len(self.model))}")
        N = len(self.p_ini)
        maxNvar = max(m.nparams for m in self.implemented_models)
        self.varbuffer = [np.zeros((len(self.model), maxNvar + self._hankel_order, self.maxnfreq))
                          for i in range(N + 1)]
        self.varfixed = [np.zeros(self.maxnfreq) for i in range(N + 1)]
        self.dpar = np.zeros(N, dtype=np.float64)
        self.par2 = np.zeros((3, N), dtype=np.float64)

        self.Hessian = np.zeros((N, N), dtype=np.float64)
        self.Gradient = np.zeros(N, dtype=np.float64)
        self.imod = np.zeros(len(self.model), dtype=np.int32)

        # Compile equations that set the model variables and gains:
        self.logger.info("going to compile models")
        self.compile_all_models()

        self.GainBuffer = [[[] for AI in range(self.Nants)] for spidx in range(nspw)]

        for spidx in range(nspw):
            for AI in range(self.Nants):
                if self.isMixed:
                    self.GainBuffer[spidx][AI] = [np.ones(
                        (len(tArr[spidx]),
                         len(averfreqs[spidx])),
                        dtype=np.complex128) for i in range(len(self.parDependence[AI]))]
                else:
                    self.GainBuffer[spidx][AI] = [np.ones(
                        len(tArr[spidx]) + len(averfreqs[spidx]),
                        dtype=np.complex128) for i in range(len(self.parDependence[AI]))]

        compFixed = len(self.fixed) > 0

        # Fill in all the information and data pointers for the modeler:
        self.logger.info("going to run setModel")
        gooduvm = uvmod.setModel(self.imod, self.Hessian, self.Gradient,
                                 self.varbuffer, self.varfixed, self.dpar,
                                 self.propRA, self.propDec, refpos,
                                 self.parDependence, self.GainBuffer, compFixed, self.isMixed)

        # Delta to estimate numerical derivatives:
        # self.minnum = self.minDp

        if not gooduvm:
            self.logger.error("Error in the C++ extension!")
            return False

        goodinit = uvmod.setWork()
        if goodinit != 10:
            self.logger.error("Memory allocation error!")

        # self.bounds = self.bounds
        # self.LMtune = [float(lmi) for lmi in self.LMtune]

        return True

    @classmethod
    def fit_range(cls, currchan, currspw, freqs):
        cls.logger.debug("Modeler::fit_range")
        nui = -1 if currchan < 0 else int(currchan)
        spwrange = list(range(len(freqs))) if currspw < 0 else [int(currspw)]
        return (nui, spwrange)

    def _perform_fit(self, write_model):
        self.logger.debug("Modeler::_perform_fit")
        if self.method == 'simplex':
            fitsimp = _mod_simplex(self.chi_square, self.p_ini,
                                   args=(self.bounds, self.p_ini),
                                   relxtol=self.SMPtune[0], relstep=self.SMPtune[1],
                                   maxiter=self.SMPtune[2] * len(self.p_ini))
            fit = [fitsimp[0], np.zeros((len(self.p_ini), len(self.p_ini))), fitsimp[1]]
        else:
            converged, fit = self.LMMin(self.p_ini)
            if converged:
                # Estimate the parameter uncertainties and save the model in the output array:
                if write_model == 1:
                    _ = self.residuals(fit[0], mode=-3)
                elif write_model == 2:
                    _ = self.residuals(fit[0], mode=-4)
                elif write_model == 3:
                    _ = self.residuals(fit[0], mode=-5)
        return fit

    def fit(self, ms, write_model, cov_return, redo_fixed=True, reset_flags=False):
        """ Fits the data, using the models previously compiled with ``init_model()``.

        Parameters
        ----------

        **redo_fixed** : `bool`
          It is True by default. If False, the fixed model will **not** be recomputed throughout
          the fit.

          .. warning:: Setting ``redo_fixed=False`` can be dangerous if you are fitting in
             spectral-line mode and have channels with very different frequencies (i.e., a wide
             fractional bandwidth), since the UV coordinates will **not** be re-projected in that
             case.

        **reset_flags** : `bool`
          Default is False. This is used to clear out the status of *bad data* that may have been
          set by special routines of **UVMultiFit** (e.g., the Quinn Fringe Fitter). Default is
          to NOT reset flags (this option should be OK most of the time)."""

        self.logger.debug("Modeler::fit")
        tic = time.time()

        # self.bounds = self.bounds
        # self.LMtune = [float(lmi) for lmi in self.LMtune]

        npars = len(self.p_ini)
        nspwtot = ms.spwlist[-1][3] + len(ms.spwlist[-1][2])

        notfit = [[] for si in range(nspwtot)]

        # Select data according to time range:
        datatot = 0
        self.allflagged = False

        for si in range(nspwtot):
            if reset_flags:
                self.fittablebool[si][:] = True
            if ms.MJDrange[0] > 0.0 and ms.MJDrange[1] > 0.0:
                self.logger.info("selecting data by Modified Julian Date")
                timeWindow = np.logical_and(ms.t[si] >= ms.MJDrange[0],
                                            ms.t[si] <= ms.MJDrange[1])
                self.fittablebool[si][:] = np.logical_and(timeWindow, self.fittablebool[si])
                del timeWindow

            self.fittable[si][:] = self.fittablebool[si][:].astype(np.int8)

        #    else:
        #
        #      self.fittable[si][:] = 1
        #      self.fittablebool[si][:] = 1

        # Check if there is data available:
            unflagged = np.sum(self.wgt[si][self.fittablebool[si], :] != 0.0, axis=0)
            if self._spectral_mode:
                if np.sum(unflagged == 0.0) > 0:
                    ch = list(np.where(unflagged == 0.0))
                    self.logger.error("not enough data for this time range, channels: {}".format(ch))
                    self.allflagged = True
                    notfit[si] = list(np.where(unflagged == 0.0)[0])
            else:
                datatot += np.sum(unflagged)

        if datatot == 0 and not self._spectral_mode:
            self.logger.error("not enough data for this time range")
            self.allflagged = True
            return None

        #  for si in range(nspwtot):
        #    self.wgtcorr[si][:] = -self.KfacWgt
        self.logger.info("now fitting model")

        # Initialize model:
        # TODO test re-initialization of model by direct call of init_model
        # if reinit_model:
        #     if self.initiated:
        #         goodinit = self.init_model()
        #         if not goodinit:
        #             self.logger.error(f"bad model (re)initialization")
        #             return None

        for i in range(len(self.p_ini) + 1):
            self.varbuffer[i][:] = 0.0
            self.varfixed[i][:] = 0.0

        # FIT!!!
        ndata = 0.0

        ##################
        # CASE OF SPECTRAL-MODE FIT:
        if self._spectral_mode:
            self.logger.debug("spectral mode fit")
            fitparams = [[] for j in range(nspwtot)]
            fiterrors = [[] for j in range(nspwtot)]
            covariance = [[] for j in range(nspwtot)]
            ChiSq = [[] for j in range(nspwtot)]
            Nvis = [[] for j in range(nspwtot)]

            # Fit channel-wise for each spw:
            for si in range(nspwtot):
                rang = np.shape(self.wgt[si])[1]
                for nuidx in range(rang):
                    self.logger.info("fitting channel {} of {} in spw {}".format(nuidx+1, rang, si))
                    self.currspw = si
                    self.currchan = nuidx

                    # Compute fixed model:
                    if redo_fixed and len(self.fixed) > 0 and not self.takeModel:
                        nui, spwrange = self.fit_range(self.currchan, self.currspw, self.freqs)
                        self._compute_fixed_model(self.par2[2, :], spwrange, nui)

                    # Fit with simplex (if asked for):

                    if nuidx not in notfit[si]:
                        fit = self._perform_fit(write_model)
                        if not fit:
                            return None
                    else:
                        fit = [[0.0 for pi in self.p_ini], np.zeros((len(self.p_ini), len(self.p_ini))), 0]

                    fitparams[si].append([float(f) for f in fit[0]])

                    # Only add the unflagged data to compute the DoF
                    ndata = float(np.sum(self.wgt[si][:, nuidx] > 0.0))

                    if ndata > 0.0:
                        ChiSq[si].append(fit[2] / ndata)   # Reduced ChiSquared
                    else:
                        # There are 0 'really-free' parameters?! Watch out!
                        ChiSq[si].append(float(fit[2]))

                    Nvis[si].append(ndata)

                    fiterrors[si].append([np.sqrt(fit[1][i, i] * ChiSq[si][nuidx]) for i in range(npars)])
                    covariance[si].append(fit[1] * ChiSq[si][nuidx])

                self.fitpars = fitparams[si]

        ##################
        # CASE OF CONTINUUM-MODE FIT:
        else:
            self.logger.info("continuum mode fit")

            # This will tell the modeller to solve in continuum mode:
            self.currspw = -1
            self.currchan = -1

            # Compute fixed model:
            if redo_fixed and len(self.fixed) > 0 and not self.takeModel:
                self.logger.debug("Generating fixed model. May take some time")
                nui, spwrange = self.fit_range(self.currchan, self.currspw, self.freqs)
                self._compute_fixed_model(self.par2[2, :], spwrange, nui)

            fit = self._perform_fit(write_model)
            if not fit:
                return None

            fitparams = fit[0]
            for si in range(nspwtot):
                # Only add the unflagged data to compute the DoF
                ndata += float(np.sum(self.wgt[si] > 0.0))

            if fit[2] > 0.0:
                ChiSq = fit[2] / ndata  # Reduced chi squared.
            else:
                # There are 0 'really-free' parameters?! Watch out!!
                ChiSq = fit[2]
            Nvis = ndata

            fiterrors = [np.sqrt(fit[1][i, i] * ChiSq) for i in range(npars)]
            covariance = fit[1] * ChiSq

        self.logger.info("the reduced Chi Squared will be set to 1 by re-scaling the visibility weights.")
        # Free some memory:

        #####
        # Set the 'result' property:
        if not self._spectral_mode:
            fit_results = {'Frequency': ms.averfreqs[0][0], 'Parameters': np.array(fitparams),
                           'Uncertainties': np.array(fiterrors), 'Reduced Chi squared': ChiSq,
                           'Fit': fit, 'Degrees of Freedom': Nvis}
        else:
            Freq = []  # {}
            Par = []  # {}
            Err = []  # {}
            Chi2 = []  # {}
            for sp in range(nspwtot):
                Freq.append(ms.averfreqs[sp])
                Par.append(np.array(fitparams[sp]))
                Err.append(np.array(fiterrors[sp]))
                Chi2.append(np.array(ChiSq[sp]))
            fit_results = {'Frequency': Freq,
                           'Parameters': Par,
                           'Uncertainties': Err,
                           'Reduced Chi squared': Chi2,
                           'Fit': fit,
                           'Degrees of Freedom': Nvis}

        if cov_return:
            fit_results['covariance'] = covariance

        tac = time.time()
        self.logger.info(f"fit took {(tac - tic):.2f} seconds")
        #  uvmod.unsetWork()

        return fit_results

if __name__ == "__main__":
    mod = Modeler(['delta'])
    mod.logger.info("logged by modeler")
    print(mod)

    mod.NCPU = 8
    mod.flux_only = False
    mod.hankel_order = 40
    mod.dump()
