import os
import shutil
import logging
import numpy as np
from casatools import image

from .uvmultifit import uvmultifit

ia = image()

class immultifit(uvmultifit):

    """Similar to uvmultifit, but the engine uses images created by the CASA task clean.

    All keywords are equal to those of uvmultifit, with the exception of these:

    :Parameters:
    ------------
    **psf** : `str`
      The image (or image cube) of the Point Spread Function. Usually, the name of this image
      is in the form ``blablabla.psf``
    **residual** : `str`
      The image (or image cube) of the clean residuals. Usually, the name is in the form ``blablabla.residual``.

      .. note:: You should clean with ``niter=0``, to ensure that no signal is deconvolved.

      .. note:: If you want primary-beam correction, you should use the ``blablabla.image``,
                instead of ``blablabla.residual``.

    **reinvert** : `bool`
      If immultifit has already been run and the user wants to refit the same data,
      the Fourier inversion of the images (or image cubes) can be avoided if reinvert=False. Default is True
    **dBcut** : `double`
      Mask the uv points with a signal lower than dBCut of the peak (in dB). Use this to minimize
      convolution-like artifacts in UV space. Default: -30dB => SNR = 1000.

    .. note:: The keywords ``vis``, ``spw``, ``timewidth``, ``chanwidth``, ``phase/amp_gain``, etc.
              are **not** used by immultifit.

    .. note:: Instead of writing models or residuals into a measurement set, ``immultifit``
              writes *images* with the same gridding as the original images used.
    """

    def __init__(self, reinvert=True, start=0, nchan=-1, psf='', residual='', dBcut=-30., **kwargs):
        """ Constructor."""

        logging.basicConfig(level=logging.DEBUG,
                            format='%(name)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger("immultifit")
        self.logger.debug("immultifit::__init__")
        self.Nspw = 1
        self.takeModel = False
        self.psf = str(psf)
        self.residual = str(residual)
        self.start = int(start)
        self.nchan = int(nchan)
        self.reinvert = bool(reinvert)
        # self.scipylm = True
        self.dBcut = dBcut
        uvmultifit.__init__(self, **kwargs)

    def writeModel(self):
        """ Redefined function to write fitting information into images.

        Instead of writing model (or residual) visibilities into measurement sets, the task make a new image
        with the post-fit residuals, with the same gridding as the image used for the fit."""

        self._printInfo("There is no measurement set (i.e., you are running 'immultifit').\n"
                        + "Will generate an image, instead.")
        if os.path.exists(self.residual + '.immultifit'):
            shutil.rmtree(self.residual + '.immultifit')
        os.system('cp -r %s %s' % (self.residual, self.residual + '.immultifit'))

        sq2 = np.sqrt(2.)

        # Set the correct total flux of model:
        q = (self.u[0]**2. + self.v[0]**2.)**0.5
        zerosp = q < 1.0
        for nui in range(self.start, self.nchan):
            self.mymodel.output[0][:, nui - self.start][zerosp] = np.max(self.mymodel.output[0][:, nui - self.start])

        ia.open(self.residual + '.immultifit')
        resim = ia.getchunk()
        imdims = np.shape(resim)
        # npix = float(imdims[0] * imdims[1])
        temparr = np.zeros((imdims[0], imdims[1]), dtype=np.complex128)
        for i in range(0, self.start):
            resim[:, :, :, i] = 0.0
        for i in range(self.nchan, imdims[-1]):
            resim[:, :, :, i] = 0.0

        for i in range(self.start, self.nchan):
            self._printInfo("Doing frequency channel %i of %i" % (i + 1 - self.start, self.nchan - self.start))
            for j in range(imdims[0]):
                for k in range(imdims[1]):
                    cpix = imdims[1] * j + k
                    temparr[j, k] = (self.mymodel.output[0][cpix, i - self.start]
                                     ) * self.averweights[0][cpix, i - self.start]
            # temparr[j, k] -= 1.j*(self.mymodel.output[0][1][cpix, i-self.start])*\
            #                       self.averweights[0][cpix, i-self.start]
            resim[:, :, self.stokes, i] -= np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(temparr))).real / sq2

        ia.putchunk(resim)
        ia.setbrightnessunit('Jy/beam')
        ia.close()

        self._printInfo("Will now save the unconvolved model image.")
        modname = '.'.join(self.residual.split('.')[:-1]) + '.fitmodel.immultifit'
        if os.path.exists(modname):
            shutil.rmtree(modname)
        os.system('cp -r %s %s' % (self.residual, modname))

        ia.open(modname)
        resim = ia.getchunk()
        imdims = np.shape(resim)
        # npix = float(imdims[0] * imdims[1])
        temparr = np.zeros((imdims[0], imdims[1]), dtype=np.complex128)
        for i in range(self.start, self.nchan):
            self._printInfo("Doing frequency channel %i of %i" % (i + 1 - self.start, self.nchan - self.start))
            for j in range(imdims[0]):
                for k in range(imdims[1]):
                    cpix = imdims[1] * j + k
                    temparr[j, k] = self.mymodel.output[0][cpix, i - self.start]
                    #    temparr[j, k] -= 1.j*(self.mymodel.output[0][1][cpix, i-self.start])
            resim[:, :, self.stokes, i] = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(temparr))).real / sq2

        ia.putchunk(resim)
        ia.setbrightnessunit('Jy/pixel')
        ia.close()

    ############################################
    #
    #  SANITY CHECKS AND PARAMETER SETTINGS
    #
    def checkInputs(self):
        """ Function re-definition for the immultifit class."""

        self._printInfo("In image mode, the whole image is taken.\n"
                        + "Will override the 'spw', 'field', 'MJDrange', and 'scan' parameters.")
        self._printInfo("In image mode, the channelization used is that of the image.\n"
                        + "Will override the 'chanwidth' and 'timewidth' parameters.")

        self.savemodel = True
        success = self._checkOrdinaryInputs()

        if not success:
            return False

        try:
            self.stokes = abs(int(self.stokes))
        except Exception:
            self._printInfo("In image mode, 'stokes' must be an integer\n(i.e., the Stokes column of the image).")
            polimag = ['I', 'Q', 'U', 'V']
            if self.stokes in polimag:
                stchan = polimag.index(self.stokes)
                self._printInfo("Since stokes is '%s', will try with image column %i" % (self.stokes, stchan))
                self.stokes = stchan
            else:
                self._printError("Please, set the 'stokes' keyword to an image pol. channel")
                return False

        try:
            ia.open(self.residual)
            imdims = np.array(ia.shape())
            imsum = ia.summary()
            ia.close()
        except Exception:
            self._printError("The residual image is not an image (or does not exist)!")
            return False

        try:
            ia.open(self.psf)
            imdims2 = np.array(ia.shape())
            ia.close()
        except Exception:
            self._printError("The PSF image is not an image (or does not exist)!")
            return False

        if max(np.abs(imdims - imdims2)) > 0:
            self._printError("The dimensions of the PSF image and the image of residuals do NOT coindice!")
            return False

        if imdims[2] > self.stokes:
            self._printInfo("Selecting stokes column %i" % self.stokes)
        else:
            self._printError("The images only have %i stokes columns, but you selected stokes=%i" %
                             (imdims[2], self.stokes))
            return False

        self.MJDrange = [0., 0.]
        self.spwlist = [[0, 0, [range(imdims[3])], 0]]

        RAd = np.where(imsum['axisnames'] == 'Right Ascension')[0][0]
        Decd = np.where(imsum['axisnames'] == 'Declination')[0][0]
        # np.copy(csys.torecord()['direction0']['crval'])
        self.refpos = np.array([imsum['refval'][RAd], imsum['refval'][Decd]])

        if (self.nchan < 0):
            self.nchan = imdims[3]

        self.nchan += self.start

        if self.start >= imdims[3]:
            self._printError("Starting channel (%i) must be smaller than number of channels (%i)" %
                             (self.start, imdims[3]))
        if self.nchan > imdims[3]:
            self._printError("Ending channel (%i) must be smaller than number of channels (%i)" %
                             (self.nchan, imdims[3]))

# Try to compile the equations for the variables:
#     if data_changed:
#       try:
#         del self.mymodel
#       except Exception:
#         pass

#       self.mymodel = modeler(self.model, self.var, self.fixed, self.fixedvar, self.scalefix,
#                              self.NCPU, self.only_flux, self.applyHankel, self.isNumerical,
#                              self.useGains, [self.phase_gains, self.amp_gains])
#     else:
#       self.mymodel.var = self.var
#       self.mymodel.model = self.model
#       self.mymodel.fixed = self.fixed
#       self.mymodel.fixedvar = self.fixedvar
#       self.mymodel.scalefix = self.scalefix
#       self.mymodel.calls = 0
#       self.mymodel.removeFixed=False

#     if self.mymodel.failed:
#       self._printError(self.mymodel.resultstring)
#       return False

        if self.proper_motion == 0.0:
            self.proper_motion = [[0., 0.] for i in self.model]
        elif not isinstance(self.proper_motion, list):
            self._printError("'proper_motion' must be a list!")

        self._setWgtEq()
        return True

    ############################################
    #
    #  NO DDE WEIGHTING
    #
    def _setWgtEq(self):
        self.mymodel.KfacWgt = np.array([0.0, 0.0])
        self.mymodel.NCPU = self.NCPU
        return True

    ############################################
    #
    #  READ DATA (IMAGES AND PSFs)
    #
    def readData(self, del_data=True):
        """ Function redefinition for the immultifit class."""
        if self.reinvert:

            os.system('rm -rf %s.fft.*' % self.residual)
            os.system('rm -rf %s.fft.*' % self.psf)

            self._printInfo("Inverting images into Fourier space")
            ia.open(self.residual)
            ia.fft(real=self.residual + '.fft.real', imag=self.residual + '.fft.imag')
            ia.close()
            ia.open(self.psf)
            ia.fft(amp=self.psf + '.fft.amp')
            ia.close()

        ia.open(self.psf + '.fft.amp')
        wgt = np.abs(ia.getchunk()[:, :, :, self.start:self.nchan])
        imdims = np.shape(wgt)
        npix2 = imdims[0] * imdims[1]
        nnu = imdims[3]

        ui = np.ones(npix2, dtype=np.float64)
        vi = np.ones(npix2, dtype=np.float64)
        wgti = np.zeros((npix2, nnu), dtype=np.float64)
        datare = np.zeros((npix2, nnu), dtype=np.float64)
        dataim = np.zeros((npix2, nnu), dtype=np.float64)
        # outpre = np.zeros((npix2, nnu), dtype=np.float64)
        # outpim = np.zeros((npix2, nnu), dtype=np.float64)
        # fixedre = np.zeros((npix2, nnu), dtype=np.float64)
        # fixedim = np.zeros((npix2, nnu), dtype=np.float64)
        freqs = np.zeros(nnu, dtype=np.float64)

        # zeros = []
        for i in range(imdims[3]):
            freqs[i] = ia.toworld([0, 0, self.stokes, i + self.start])['numeric'][3]

        self._printInfo("Reading gridded UV coordinates and weights")
        u0, v0, _, _ = ia.toworld([0, 0, self.stokes, self.start])['numeric']
        u1, v1, _, _ = ia.toworld([1, 1, self.stokes, self.start])['numeric']
        du = u1 - u0
        dv = v1 - v0
        for j in range(imdims[0]):
            self._printInfo("Reading row %i of %i" % (j, imdims[0]))
            for k in range(imdims[1]):
                cpix = imdims[1] * j + k
                ui[cpix] = (u0 + du * j)
                vi[cpix] = (v0 + dv * k)
                wgti[cpix, :] = np.copy(wgt[j, k, self.stokes, :])

        del wgt
        ia.close()
        # Maximum dynamic range set by "Window effect"
        maskwgt = wgti < np.max(wgti) * (10.**(self.dBcut / 10.))
        self._printInfo("Reading gridded visibilities")
        ia.open(self.residual + '.fft.real')
        datas = ia.getchunk()
        for j in range(imdims[0]):
            self._printInfo("Reading row %i of %i" % (j, imdims[0]))
            for k in range(imdims[1]):
                cpix = imdims[1] * j + k
                datare[cpix, :] = datas[j, k, self.stokes, self.start:self.nchan]

        del datas
        ia.close()

        ia.open(self.residual + '.fft.imag')
        datas = ia.getchunk()
        for j in range(imdims[0]):
            self._printInfo("Reading row %i of %i" % (j, imdims[0]))
            for k in range(imdims[1]):
                cpix = imdims[1] * j + k
                dataim[cpix, :] = -datas[j, k, self.stokes, self.start:self.nchan]

        del datas
        ia.close()

        toLamb = self.LtSpeed / freqs[0]
        self.u = [ui * toLamb]
        self.v = [vi * toLamb]

        self.averdata = [(datare + 1.j * dataim) / wgti]
        bads = np.logical_or(np.isnan(self.averdata[0].real), np.isnan(self.averdata[0].imag))
        badmsk = np.logical_or(bads, maskwgt)
        self.averdata[0][badmsk] = 0.0
        # self.averdata[0][1][badmsk] = 0.0
        wgti[badmsk] = 0.0
        self.averweights = [wgti]
        self.iscancoords = [[[0, -1, 0, 0]]]
        self.averfreqs = [freqs]
        self.t = [np.linspace(0., 1., npix2)]

        # DUMMY VALUES:
        self.RAshift = [np.zeros(npix2)]
        self.Decshift = [np.zeros(npix2)]
        self.ant1 = [np.zeros(npix2, dtype=np.int32)]
        self.ant2 = [np.ones(npix2, dtype=np.int32)]
        self.w = [np.zeros(npix2)]
        self.Stretch = [np.ones(npix2)]
        self.t0 = 0.0
        self.t1 = 1.e12
        self.tArr = [np.linspace(0., 1., npix2)]
        self.tIdx = [np.arange(npix2, dtype=np.int32)]

        self.Nants = 2

        self._printInfo("Done reading.")

        self.initData()

        return True
