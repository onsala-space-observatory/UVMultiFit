import sys
import re
import logging
import numpy as np
from scipy import special

class modeler():
    """ Class to deal with model equations and fitting procedures.

    If convert strings representing models and parameters into compiled equations, to be used in a ChiSq
    visibility fitting. It also interacts with the C++ extension of UVMultiFit.

    This class should NOT be instantiated by the user (it is called from the UVMultiFit class."""

    logger = logging.getLogger("modeler")

    ############################################
    #
    #  FREE MEMORY
    #
    # def __del__(self):
    #     self.deleteData()
    #     self.deleteModel()

    ############################################
    #
    #  FREE MEMORY JUST FOR THE MODEL-RELATED DATA:
    #
    # def deleteModel(self):
    #     """ Free pointers to the model-related arrays and parameters."""
    #
    #     for mdi in range(len(self.varbuffer) - 1, -1, -1):  # [::-1]:
    #         del self.varbuffer[mdi]
    #     #    del self.varbuffer
    #
    #     for mdi in range(len(self.varfixed) - 1, -1, -1):  # [::-1]:
    #         del self.varfixed[mdi]
    #     #    del self.varfixed
    #
    #     #    for mdi in range(len(self.dpar)-1, -1, -1):  # -[::-1]:
    #     #      del self.dpar[mdi]
    #     del self.dpar
    #
    #     #    for mdi in range(len(self.par2)-1, -1, -1):  # [::-1]:
    #     #      del self.par2[mdi]
    #     del self.par2
    #
    #     del self.Hessian, self.Gradient, self.imod

    # def deleteData(self):
    #     """ Free pointers to the data-related arrays and gain buffers."""
    #
    #     for mdi in range(len(self.data) - 1, -1, -1):  # [::-1]:
    #         del self.data[mdi]
    #     #    del self.data
    #
    #     for mdi in range(len(self.wgt) - 1, -1, -1):  # [::-1]:
    #         del self.wgt[mdi]
    #     #    del self.wgt
    #
    #     for mdi in range(len(self.uv) - 1, -1, -1):  # [::-1]:
    #         try:
    #             del self.uv[mdi][2], self.uv[mdi][1], self.uv[mdi][0]
    #             del self.uv[mdi]
    #         except Exception:
    #             pass
    #     #    del self.uv
    #
    #     for mdi in range(len(self.offset) - 1, -1, -1):  # [::-1]:
    #         try:
    #             del self.offset[mdi][2], self.offset[mdi][1], self.offset[mdi][0]
    #             del self.offset[mdi]
    #         except Exception:
    #             pass
    #
    #     #    del self.offset
    #     for mdi in range(len(self.ants) - 1, -1, -1):  # [::-1]:
    #         try:
    #             del self.ants[mdi][1], self.ants[mdi][0]
    #             del self.ants[mdi]
    #         except Exception:
    #             pass
    #     #    del self.ants
    #     for mdspw in range(len(self.GainBuffer) - 1, -1, -1):  # [::-1]:
    #         NA = len(self.GainBuffer[mdspw])
    #         for a in range(NA - 1, -1, -1):
    #             NP = len(self.GainBuffer[mdspw][a])
    #             for mdp in range(NP - 1, -1, -1):
    #                 del self.GainBuffer[mdspw][a][mdp]
    #             del self.GainBuffer[mdspw][a]
    #         del self.GainBuffer[mdspw]
    #     #    del self.GainBuffer
    #
    #     for mdi in range(len(self.dt) - 1, -1, -1):  # [::-1]:
    #         del self.dt[mdi]
    #     #    del self.dt
    #
    #     for mdi in range(len(self.dtArr) - 1, -1, -1):  # [::-1]:
    #         del self.dtArr[mdi]
    #     #    del self.dtArr
    #
    #     for mdi in range(len(self.dtIdx) - 1, -1, -1):  # [::-1]:
    #         del self.dtIdx[mdi]
    #     #    del self.dtIdx
    #
    #     for mdi in range(len(self.output) - 1, -1, -1):  # [::-1]:
    #         del self.output[mdi]
    #     #    del self.output
    #
    #     for mdi in range(len(self.freqs) - 1, -1, -1):  # [::-1]:
    #         del self.freqs[mdi]
    #
    #     for mdi in range(len(self.fittable) - 1, -1, -1):  # self.fittable[::-1]:
    #         del self.fittable[mdi]
    #     #    del self.fittable
    #
    #     for mdi in range(len(self.wgtcorr) - 1, -1, -1):  # [::-1]:
    #         del self.wgtcorr[mdi]
    #     #    del self.wgtcorr
    #
    #     for mdi in range(len(self.fittablebool) - 1, -1, -1):  # [::-1]:
    #         del self.fittablebool[mdi]
    #     #    del self.fittablebool
    #
    #     for mdi in range(len(self.isGain) - 1, -1, -1):  # [::-1]:
    #         del self.isGain[mdi]
    #     #    del self.isGain
    #
    #     for mdi in range(len(self.iscancoords) - 1, -1, -1):  # [::-1]:
    #         Npar = len(self.iscancoords[mdi])
    #         for mdp in range(Npar - 1, -1, -1):
    #             del self.iscancoords[mdi][mdp]
    #         del self.iscancoords[mdi]
    #     #    del self.iscancoords

    ############################################
    #
    #  CREATE INSTANCE
    #
    def __init__(self):
        """ Just the constructor of the 'modeler' class."""
        logging.basicConfig(level=logging.INFO,
                            format='%(name)s - %(levelname)s - %(message)s')
        self.logger.debug("modeler::__init__")
        self.Nants = 0
        self.initiated = False
        self.addfixed = False
        self.expka = 2. * np.log(2.)
        self.pow2ka = 1. / (np.pi * np.log(2.))
        self.pow3ka = np.sqrt(2.**(2. / 3.) - 1.)
        self.FouFac = (2. * np.pi) * np.pi / 180. / 3600.
        self.LtSpeed = 2.99792458e+8
        self.deg2rad = np.pi / 180.
        self.failed = False
        self.calls = 0   # Total number of calls during the fit
        self.minnum = np.sqrt(np.finfo(np.float64).eps)
        # Tells us if the fixedmodel array exists and should be used.
        self.removeFixed = False
        self.varfixed = []
        self.dpar = []
        self.Hessian = np.zeros((1, 1), dtype=np.float64)
        self.Gradient = []
        self.iscancoords = []
        self.varbuffer = []
        self.isGain = []
        self.fittablebool = []
        self.fittable = []
        self.wgtcorr = []
        self.dt = []
        self.dtArr = []
        self.dtIdx = []
        self.imod = []
        self.par2 = []
        # Arrays of data and pointers (to share access with C library):
        self.data = []
        self.dt = []
        self.dtArr = []
        self.dtIdx = []
        self.wgt = []
        self.wgtcorr = []
        self.iscancoords = []
        self.uv = []
        self.offset = []
        self.output = []
        self.freqs = []
        self.fixedmodel = []
        self.fittable = []
        self.fittablebool = []
        self.ants = []
        # Buffer arrays to save the values of the variables:
        self.varbuffer = []
        # Model indices (to let the C++ library know which model is which component):
        self.imod = []
        self.ifixmod = []
        # spw and channel to fit (-1 means fit to the continuum):
        self.currspw = 0
        self.currchan = 0
        # C++ library:
        self.Ccompmodel = lambda spw, nui, mode: 0.0
        # Parameters computed in unbound space:
        self.par2 = []
        self.bounds = []
        # Levenberg-Marquardt parameters:
        self.LMtune = []
        # Gains:
        self.GainBuffer = []
        #   self.phaseGainBuffer = []
        self.parDependence = []
        self.isGain = []
        self.strucvar = []
        KGaus = np.sqrt(1. / (4. * np.log(2.)))
        # Some useful functions:
        self._compiledScaleFixed = lambda p, nu: 1.0 + 0.0
        self.LorentzLine = lambda nu, nu0, P, G: P * 0.25 * G * G / (np.power(nu - nu0, 2.) + (0.5 * G)**2.)
        self.GaussLine = lambda nu, nu0, P, G: P * np.exp(-np.power((nu - nu0) / (G * KGaus), 2.))

        self.pieceWise = lambda t, p0, p1, t0, t1: np.clip(p0 + (p1 - p0) * (t - t0) / (t1 - t0), p0, p1)

        self.wgtEquation = lambda D, Kf: -D * Kf
        self.KfacWgt = 1.0

        # List of currently-supported model components:
        self.allowedmod = ['delta', 'Gaussian', 'disc', 'ring', 'sphere',
                           'bubble', 'expo', 'power-2', 'power-3', 'GaussianRing']
        self.resultstring = ''

    def setup(self, model, parameters, fixedmodel, fixedparameters, scalefix, NCPU, only_flux, HankelOrder,
              isNumerical, useGains, gainFunction, isMixed):
        """ Setup the model components, fitting parameters, and gain functions. Compile the equations.

        Not to be called by the user."""

        self.logger.debug("modeler::setup")

        # self.minnum = np.finfo(np.float64).eps  # Step for Jacobian computation
        self.propRA = np.zeros(len(model), dtype=np.float64)
        self.propDec = np.zeros(len(model), dtype=np.float64)
        self.fixed = fixedmodel
        self.fixedvar = fixedparameters
        self.model = model
        self.var = parameters  # variables of the model components.
        self.scalefix = scalefix
        self.HankelOrder = HankelOrder
        self.isNumerical = isNumerical
        self.isMixed = isMixed

        # Lists of compiled functions (one per model component).
        # These will take p and nu and return the variables of the model:

        self.varfunc = [0.0 for component in model]
        self.fixedvarfunc = [0.0 for component in fixedmodel]

        self.NCPU = NCPU
        self.t0 = 0.0
        self.t1 = 1.e12
        self.only_flux = only_flux
        self.useGains = useGains

        self.gainFunction = gainFunction

    ############################################
    #
    #  COMPILE AT RUNTIME
    #
    # Compile the functions to make p, nu -> var:
    def compileAllModels(self):
        """ Compile all models (fixed, variable, and fixed-scale. Not to be called directly by the user. """
        self.logger.debug("modeler::compileAllModels")
        self.resultstring = ''
        self.failed = False
        self._compileModel()
        self._compileFixedModel()
        self._compileScaleFixed()
        self._compileGains()

    def executeCode(self, code):
        try:
            self.logger.info(code)
            exec(code, locals())
            return True
        except Exception:
            return False

    #  print self.gainFunction
    def _compileGains(self):
        """ Compile the functions related to the antenna gains."""
        self.logger.debug("modeler::_compileGains")
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
            for ni in self.gainFunction[0].keys():
                tempstr = self.gainFunction[0][ni].replace(
                    'pieceWise(', 'self.pieceWise(t, ').replace(
                        't', 't[:,__import__(\'numpy\').newaxis]').replace(
                            'nu0', '%.12f' % self.freqs[0][0])
                modstr = 'self.phaseAntsFunc[' + str(ni) + '] = lambda t, nu, p: ' + tempstr

                parpos = [x.start() for x in re.finditer(r'p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                if not self.executeCode(modstr):
                    self.failed = True
                    self.resultstring = 'Syntax error in phase gain of antenna %i' % (ni)
                    return

            for ni in self.gainFunction[1].keys():
                tempstr = self.gainFunction[1][ni].replace(
                    'pieceWise(', 'self.pieceWise(t,').replace(
                        't', 't[:,__import__(\'numpy\').newaxis]').replace(
                            'nu0', '%.12f' % self.freqs[0][0])

                modstr = 'self.ampAntsFunc[' + str(ni) + '] = lambda t, nu, p: ' + tempstr

                parpos = [x.start() for x in re.finditer(r'p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                if not self.executeCode(modstr):
                    self.failed = True
                    self.resultstring = 'Syntax error in amp gain of antenna %i' % (ni)
                    return

        else:

            for ni in self.gainFunction[2].keys():
                tempstr = self.gainFunction[2][ni].replace('nu0', '%.12f' % self.freqs[0][0])
                if "t" in tempstr:
                    self.failed = True
                    self.resultstring = 'A frequency-dependent gain cannot depend on time!'
                    return

                modstr = 'self.phaseAntsFuncNu[' + str(ni) + '] = lambda nu, p: ' + tempstr

                parpos = [x.start() for x in re.finditer(r'p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                if not self.executeCode(modstr):
                    self.failed = True
                    self.resultstring = 'Syntax error in phase gain of antenna %i' % (ni)
                    return

            for ni in self.gainFunction[3].keys():
                tempstr = self.gainFunction[3][ni].replace('nu0', '%.12f' % self.freqs[0][0])
                if "t" in tempstr:
                    self.failed = True
                    self.resultstring = 'A frequency-dependent gain cannot depend on time!'
                    return

                modstr = 'self.ampAntsFuncNu[' + str(ni) + '] = lambda nu, p: ' + tempstr

                parpos = [x.start() for x in re.finditer(r'p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                if not self.executeCode(modstr):
                    self.failed = True
                    self.resultstring = 'Syntax error in amp gain of antenna %i' % (ni)
                    return

            for ni in self.gainFunction[4].keys():
                tempstr = self.gainFunction[4][ni].replace(
                    'pieceWise(', 'self.pieceWise(t,').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

                if "nu" in tempstr:
                    self.failed = True
                    self.resultstring = 'A time-dependent gain cannot depend on frequency!'
                    return

                modstr = 'self.phaseAntsFuncT[' + str(ni) + '] = lambda t, p: ' + tempstr

                parpos = [x.start() for x in re.finditer(r'p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                if not self.executeCode(modstr):
                    self.failed = True
                    self.resultstring = 'Syntax error in phase gain of antenna %i' % (ni)
                    return

            for ni in self.gainFunction[5].keys():
                tempstr = self.gainFunction[5][ni].replace(
                    'pieceWise(', 'self.pieceWise(t,').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

                if "nu" in tempstr:
                    self.failed = True
                    self.resultstring = 'A time-dependent gain cannot depend on frequency!'
                    return

                modstr = 'self.ampAntsFuncT[' + str(ni) + '] = lambda t, p: ' + tempstr

                parpos = [x.start() for x in re.finditer(r'p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                if not self.executeCode(modstr):
                    self.failed = True
                    self.resultstring = 'Syntax error in amp gain of antenna %i' % (ni)
                    return

        for ni in range(self.Nants):
            self.parDependence[ni] = np.unique(self.parDependence[ni]).astype(np.int32)

    def _compileModel(self):
        """ Compiles the variable model, according to the contents of the 'model' and 'parameters' lists."""

        self.logger.debug("modeler::_compileModel")
        # Define variable model:
        self.varfunc = [0.0 for component in self.model]

        for ii, component in enumerate(self.model):
            # print(ii, component, self.var[ii])
            # print(self.freqs)
            tempstr = self.var[ii].replace(
                'LorentzLine(', 'self.LorentLine(nu, ').replace(
                    'GaussLine(', 'self.GaussLine(nu, ').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            if self.only_flux:
                try:
                    tempstr2 = tempstr.split(',')
                    if len(tempstr2) > 3:
                        self.strucvar.append(list(map(float, tempstr2[:2] + tempstr2[3:])))
                    else:
                        self.strucvar.append(list(map(float, tempstr2[:2])))
                except Exception:
                    print(tempstr.split(','))
                    self.resultstring = ' If only_flux=True, all variables but the flux must be constants! Aborting!'
                    self.failed = True

            modstr = 'self.varfunc[' + str(ii) + '] = lambda p, nu: [' + tempstr + ']'
            if not self.executeCode(modstr):
                self.failed = True
                self.resultstring = 'Syntax error in component number %i of the variable model' % (ii)
                return

            if component not in self.allowedmod:
                self.resultstring = " Component '" + component + "' is unknown. Aborting!"
                self.failed = True
                return
            self.imod[ii] = self.allowedmod.index(component)

    def _compileFixedModel(self):
        """ Compiles the fixed model, according to the contents of the 'fixed' and 'fixedpars' lists."""

        self.logger.debug("modeler::_compileFixedModels")
        if len(self.fixed) > 0 and 'model_column' in self.fixed:
            return

        # Define fixed model:
        self.fixedvarfunc = [0.0 for component in self.fixed]

        self.ifixmod = np.zeros(len(self.fixed), dtype=np.int32)  # []

        for ii, component in enumerate(self.fixed):
            tempstr = self.fixedvar[ii].replace(
                'LorentzLine(', 'self.LorentzLine(nu, ').replace(
                    'GaussLine(', 'self.GaussLine(nu, ').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            modstr = 'self.fixedvarfunc[' + str(ii) + '] = lambda p, nu: [' + tempstr + ']'
            if not self.executeCode(modstr):
                self.resultstring = 'Syntax error in component number %i of the fixed model' % (ii)
                self.failed = True
                return

            if component not in self.allowedmod:
                self.resultstring = " Component '" + component + "' is unknown. Aborting!"
                self.failed = True
                return
            self.ifixmod[ii] = self.allowedmod.index(component)

    def _compileScaleFixed(self):
        """ Compiles the scaling factor for the fixed model """

        self.logger.debug("modeler::_compileScaleFixed")
        tempstr = self.scalefix.replace(
            'LorentzLine(', 'self.LorentzLine(nu, ').replace(
                'GaussLine(', 'self.GaussLine(nu, ').replace(
                    'nu0', '%.12f' % self.freqs[0][0])

        scalefixedstr = 'self._compiledScaleFixed = lambda p, nu: ' + tempstr + ' + 0.0'
        if not self.executeCode(scalefixedstr):
            self.resultstring = 'Syntax error in the flux-scale equation'
            self.failed = True

    ############################################
    #
    #  GET UNBOUND PARAMETER SPACE FROM BOUND PARAMETERS
    #
    def getPar2(self, mode=0):
        """ Function to change fitting parameters to/from the unbound space from/to the bound space.

        It also comptues the gradient for the Jacobian matrix. Unbound values are in index 0, gradient is
        in index 1, and bound values are in index 2. The equations for the changes of variables are taken
        from the MINUIT package."""

        self.logger.debug("modeler::getPar2")
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
                if self.bounds[i] is None or (self.bounds[i][0] is None and self.bounds[i][1] is None):
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
                if self.bounds[i] is None or (self.bounds[i][0] is None and self.bounds[i][1] is None):
                    self.par2[2, i] = p
                elif self.bounds[i][0] is None and self.bounds[i][1] is not None:
                    self.par2[2, i] = self.bounds[i][1] + 1. - np.sqrt(p**2. + 1.)
                elif self.bounds[i][0] is not None and self.bounds[i][1] is None:
                    self.par2[2, i] = self.bounds[i][0] - 1. + np.sqrt(p**2. + 1.)
                else:
                    Kbound = (self.bounds[i][1] - self.bounds[i][0]) / 2.
                    self.par2[2, i] = self.bounds[i][0] + Kbound * (np.sin(p) + 1.)

        for i, p in enumerate(par2):
            if self.bounds[i] is None or (self.bounds[i][0] is None and self.bounds[i][1] is None):
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
        """ Implementation of the Levenberg-Marquardt algorithm. Not to be called directly by the user. """

        self.logger.debug("modeler::LMMin")
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
        self.getPar2()

        Hessian2 = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        Gradient2 = np.zeros(np.shape(self.Gradient), dtype=np.float64)
        HessianDiag = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupHess = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupGrad = np.zeros(np.shape(self.Gradient), dtype=np.float64)
        Inverse = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupP = np.copy(self.par2[0, :])

        self.Hessian[:, :] = 0.0
        self.Gradient[:] = 0.0

        if self.only_flux and self.currchan == -1:
            nnu = max([len(self.freqs[sp]) for sp in range(len(self.freqs))])
            self.logger.info("Computing structure-only parameters")

            # for midx in range(len(p)):
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
            self.getPar2(mode=1)
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
                self.getPar2(mode=1)
                break

            if CurrChi < Chi2:
                self.par2[0, :] = backupP
                Hessian2[:, :] = backupHess
                Gradient2[:] = backupGrad
                self.getPar2(mode=1)
            else:
                CurrChi = Chi2
                backupHess[:, :] = Hessian2
                backupGrad[:] = Gradient2
                backupP[:] = self.par2[0, :]

        if controlIter == NITER:
            sys.stdout.write("REACHED MAXIMUM NUMBER OF ITERATIONS!\n"
                             + "The algorithm may not have converged!\n"
                             + "Please, check if the parameter values are meaningful.\n")
            sys.stdout.flush()

        self.getPar2(mode=1)

        try:
            return [self.par2[2, :], np.linalg.pinv(self.Hessian), Chi2]
        except Exception:
            return False

    ############################################
    #
    #  NON-ALGEBRAIC MODELS
    #
    # Compute elements of Taylor expansion of the source's Hankel transform:
    def gridModel(self, imod, tempvar):
        """ Compute elements of Taylor expansion of the source's Hankel transform."""

        self.logger.debug("modeler::gridModel")
        n = self.HankelOrder - 1

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
        raise ValueError("Model %i was not correctly interpreted!\n" % imod)

    ############################################
    #
    #  COMPUTE RESIDUALS FOR A MODEL REALIZATION
    #
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

        self.logger.debug("uvmultifit residuals")
        #  varsize = self.maxNvar+self.HankelOrder
        if mode in [0, -3]:
            self.calls = 0
        else:
            self.calls += 1

        if self.currchan < 0:  # Continuum (i.e., all data modelled at once)
            nui = -1
        else:  # Spectral mode (i.e., model only the current channel of the current spw)
            nui = int(self.currchan)

        if self.currspw < 0:  # Continuum (i.e., all data modelled at once)
            spwrange = range(len(self.freqs))
        else:  # Spectral mode (i.e., model only the current channel of the current spw)
            spwrange = [int(self.currspw)]

        # DELTA FOR DERIVATIVES:
        for j, p in enumerate(pars):
            #  self.dpar[j] = self.minnum
            if np.abs(p) < 1.0:
                self.dpar[j] = self.minnum
            else:
                self.dpar[j] = np.abs(p) * self.minnum

        # Compute antenna-gain corrections:
        if len(self.useGains) > 0 and mode != 0:
            for spw in spwrange:
                Nchan = len(self.freqs[spw])
                Nint = len(self.dtArr[spw])
                for i in self.useGains:
                    for j in range(len(self.parDependence[i])):
                        ptemp = list(pars)   # [pi for pi in pars]
                        if j > 0:
                            j2 = self.parDependence[i][j] - 1
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

        if mode == 0:  # Just compute the fixed model and return
            self.removeFixed = True
            # isfixed = True
            modbackup = self.imod[0]
            for spw in spwrange:
                self.output[spw][:] = 0.0
                for midx, mi in enumerate(self.ifixmod):
                    tempvar = self.fixedvarfunc[midx](pars, self.freqs[spw])
                    self.imod[0] = mi
                    if self.imod[0] in self.isNumerical:
                        tempvar = self.gridModel(self.imod[0], tempvar)
                    for i, tv in enumerate(tempvar):
                        nnu = len(self.freqs[spw])
                        self.varbuffer[0][0, i, :nnu] = tv
                    self.Ccompmodel(spw, nui, 0)

            self.imod[0] = modbackup
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
                scalefx = self._compiledScaleFixed(pars, self.freqs[spw])
                self.varfixed[0][:nnu] = scalefx
                for j in range(len(pars)):
                    ptemp = list(pars)   # [pi for pi in pars]
                    ptemp[j] += self.dpar[j]
                    # Variables of current component
                    scalefx = self._compiledScaleFixed(ptemp, self.freqs[spw])
                    self.varfixed[j + 1][:nnu] = scalefx
            else:
                scalefx = self._compiledScaleFixed(pars, self.freqs[spw][nui])
                self.varfixed[0][0] = scalefx
                for j in range(len(pars)):
                    ptemp = list(pars)   # [pi for pi in pars]
                    ptemp[j] += self.dpar[j]
                    # Variables of current component
                    scalefx = self._compiledScaleFixed(ptemp, self.freqs[spw][nui])
                    self.varfixed[j + 1][0] = scalefx

            for midx, modi in enumerate(currmod):
                ptemp = [float(pi) for pi in pars]

                if nui == -1:  # Continuum

                    if self.only_flux:
                        currflux = pars[midx]
                        # nstrucpars = len(self.strucvar[midx])
                        self.varbuffer[0][midx, 2, :nnu] = currflux

                    else:
                        # Variables of current component
                        tempvar = currvar[midx](pars, self.freqs[spw])
                        if modi in self.isNumerical:
                            tempvar = self.gridModel(modi, tempvar)

                        for i, tv in enumerate(tempvar):
                            self.varbuffer[0][midx, i, :nnu] = tv

                    if self.only_flux:
                        ptemp[midx] += self.dpar[midx]
                        self.varbuffer[midx + 1][midx, 2, :nnu] = ptemp[midx]

                    else:

                        for j in range(len(pars)):
                            ptemp = list(pars)   # [pi for pi in pars]
                            ptemp[j] += self.dpar[j]
                        # Variables of current component
                            tempvar = currvar[midx](ptemp, self.freqs[spw])
                            if modi in self.isNumerical:
                                tempvar = self.gridModel(modi, tempvar)
                            for i, tv in enumerate(tempvar):
                                self.varbuffer[j + 1][midx, i, :nnu] = tv

                else:  # Spectral mode

                    # Variables of current component
                    tempvar = currvar[midx](pars, self.freqs[spw][nui])
                    if modi in self.isNumerical:
                        tempvar = self.gridModel(modi, tempvar)
                    for i, tv in enumerate(tempvar):
                        self.varbuffer[0][midx, i, 0] = tv

                    for j in range(len(pars)):
                        ptemp = list(pars)    # [pi for pi in pars]
                        ptemp[j] += self.dpar[j]  # self.minnum
                        # Variables of current component
                        tempvar = currvar[midx](ptemp, self.freqs[spw][nui])
                        if modi in self.isNumerical:
                            tempvar = self.gridModel(modi, tempvar)
                        for i, tv in enumerate(tempvar):
                            self.varbuffer[j + 1][midx, i, 0] = tv

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
                msg = "Iteration # %i: achieved ChiSq: %.8e " % (self.calls, ChiSq)
                self.logger.info(msg)

        if ChiSq <= 0.0:
            raise ValueError("Invalid Chi Square!"
                             + "Maybe the current fitted value of flux (and/or size) is negative!"
                             + "Please, set BOUNDS to the fit!")

        if mode in [-1, -2, -3]:
            if dof:
                self.calls = 0
                return [ChiSq, ndata]
            return ChiSq
        return None

    ############################################
    #
    #  COMPUTE QUI SQUARE FOR A MODEL REALIZATION (NO DERIVATIVES)
    #
    def ChiSquare(self, p, bounds=None, p_ini=[]):
        """ Just a wrapper of the 'residuals' function, usually called by simplex."""

        self.logger.debug("uvmultifit ChiSquare")
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
