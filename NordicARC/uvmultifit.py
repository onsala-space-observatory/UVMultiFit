# Copyright (c) Ivan Marti-Vidal 2012.
#               EU ALMA Regional Center. Nordic node.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# a. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# b. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
# c. Neither the name of the author nor the names of contributors may
#    be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""A CASA-based flexible visibility-fitting engine developed at the Nordic node
of the ALMA Regional Center.

Current development is partially founded by the *European Union's
Horizon 2020 research and innovation programme* under grant agreement
No **730562 [RadioNet]**, in the frame of the **RINGS Working
Package** (*Radio Interferometry New Generation Software*).

This program can be used to fit multi-component models to the
visibilities in one (or several) measurement set(s). The fits can be
performed to the continuum or in spectral-line mode. Advanced model
structures can be implemented easily and several different fits can be
performed with no need of re-loading the data.

===========
Basic Usage
===========

THESE ARE ALL THE KEYWORDS (AND THEIR DEFAULTS) OF THE CURRENT
UVMULTIFIT VERSION:

>>> uvmultifit(vis='', spw='0', column = 'data', field = 0, scans = [],
>>>            uniform=False, chanwidth = 1, timewidth = 1, stokes = 'I',
>>>            write='', MJDrange=[-1.0,-1.0], model=['delta'],
>>>            var=['p[0],p[1],p[2]'], p_ini=[0.0,0.0,1.0], phase_center = '',
>>>            fixed=[], fixedvar=[], scalefix='1.0', outfile = 'modelfit.dat',
>>>            NCPU=4, pbeam=False, ldfac = 1.22, dish_diameter=0.0, gridpix=0,
>>>            OneFitPerChannel=True, bounds=None, cov_return=False,
>>>            finetune=False, uvtaper=0.0, method='levenberg', wgt_power=1.0,
>>>            LMtune=[1.e-3,10.,1.e-5,200,1.e-3], SMPtune=[1.e-4,1.e-1,200],
>>>            only_flux=False, proper_motion = 0.0, HankelOrder = 80,
>>>            amp_gains = {}, phase_gains = {})

Let us assume that you have imported this module with the alias
``uvm`` (which is the case if you followed the installation
instructions). To fit the flux density of a point source (centered at
the coordinates of the field called ``M81``), using the continuum
emission in all the spws of a measurement set called ``myms.ms``, just
do:

>>> myfit = uvm.uvmultifit(vis="myms.ms", spw="", model="delta",
>>>                        var="0, 0, p[0]", field="M81",
>>>                        OneFitPerChannel=False,
>>>                        outfile="results.dat")

The program will write the fitting result (and other useful metadata)
into the file called ``results.dat``. But it will also return a
so-called *uvmultifit instance* (which, in this example, has been
called ``myfit``), that can be used to do a lot of things with your
modelling. For instance, the fitting result (together with
uncertainties, reduced Chi Square and covariance matrix, if the fit
would have more than one fitting parameter) is stored in a dictionary
that can be accessed as:

>>> myfit.result

This dictionary has several keys worth noticing:

* ``myfit.result['Parameters']`` The fitting parameters. If the fit is done in
  spectral-line mode, these are organized per spw and channel.

* ``myfit.result['Frequency']`` The frequency corresponding to each
  set of fitting parameters (useful for cases of fits in spectral-line
  mode).

* ``myfit.result['Uncertainties']`` The uncertainties of the fitted parameters.

  .. note:: The uncertainties are estimated from the Jacobian matrix,
          and scaled so that the reduced Chi squared equals
          unity. Null (or ridicously small) uncertainties may indicate
          an unsuccessful fit. It is always a good idea to take a look
          at the post-fit residuals to assess the quality of your fit.
          .. note:: The Jacobian is computed using numerical
          approximations of the derivatives of the fitting functions
          w.r.t. the parameters.

* ``myfit.result['Reduced Chi Square']`` The reduced Chi Squared,
  before rescaling the uncertainties.

  .. note:: Notice that this is the reduced post-fit Chi squared,
            computed using the *original* visibility uncertainties
            (i.e., those read from the data). In other words, this is
            the reduced Chi Square computed *before* UVMultiFit scales
            the uncertainties to make the reduced Chi squared equal to
            unity.  The user can check with this value whether the
            visibility uncertainties are well scaled to the natural
            scatter of the data.

  .. note:: Quite large values of the reduced Chi Squared are likely
            indicative of too high data weights, which should be
            harmless to the fit, as long as the *relative* weights
            among visibilities are OK (in other words, a fit with
            *UVMultiFit* is insensitive to any global scaling factor
            of the weights).

  .. note:: Too high values of the reduced Chi Squared could also be
            indicative of a bad fit. In such a case, though, the
            a-posteriori parameter uncertainties would also be
            high. The user should check the parameter uncertainties as
            an extra assessment of the quality of the fit.

* ``myfit.result['covariance']`` The full covariance matrix of the
  parameters (in case the user asked for it).

.. note:: There are other elements in your *uvmultifit instance* that
   may be interesting to you, if you want to do more advanced
   stuff. Look at the section **Some Useful Methods** for details.

If you made a fit in *spectral-line mode* and want to plot the first
fitted parameter (i.e., p[0], see below for the syntax details)
against frequency, for the third spectral window in the measurement
set, the command could be (se assume that ``pylab`` has been
imported):

>>> pylab.plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3][:, 0])

To plot the second fitted parameter (i.e., p[1]), just execute:

>>> pylab.plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3][:, 1])

and so on. If there is only 1 fitting parameter:

>>> pylab.plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3])

.. note:: The fitted parameters are ordered by spectral window *in the
     order given by the 'spw' parameter*. Hence, if spw='2, 3', then
     the first element of the list of parameters (i.e.,
     ``myfit.result['Parameters'][0]``) will be the parameters fitted
     for the spw number 2.

.. note:: For fits to the continuum (i.e., all channels together, see
     below), ``myfit.result['Parameters']`` will only have *one entry*
     (i.e., the parameter values fitted to all channels at once)

==============
Model Examples
==============

In the following examples, we only give the most important keywords to
setup the fit.  The rest of keywords (e.g., the name of the
measurement set, the fields, the range of spectral windows, etc.) are
not important to understand the model setup and are thus avoided for
clarity.

Details about the variables of the models and how to set fitting
parameters are given in the sections below.

- **EXAMPLE 0**:  A point source with free position and free flux density:

  >>> model = ['delta']
  >>> var = ['p[0], p[1], p[2]']

  In this case, the RA shift w.r.t. the image center is ``p[0]`` (in
  arcsec), the Dec shift is ``p[1]`` (also in arcsec), and the flux
  density is ``p[2]`` (in Jy).  If we know that the source has to be
  located within 1 arcsec of the phase center, we can add this
  information as a set of bounds, i.e.:

  >>> bounds = [[-1, 1], [-1, 1], None]

  If we also want to force the flux density to be positive:

  >>> bounds = [[-1, 1], [-1, 1], [0, None]]

- **EXAMPLE 1**: Two deltas, being the position of the second one
  fixed w.r.t. the position of the first one. Let's say that the
  second delta is shifted at 0.5 and 0.6 arcsec, in RA and Dec
  (respectively), from the first delta, and we want the position of
  the first delta to be free in our fit. Then:

  >>> model = ['delta', 'delta']
  >>> var = ['p[0], p[1], p[2]', 'p[0]+0.5, p[1]+0.6, p[3]']

  In this case, ``p[0]`` and ``p[1]`` are the RA and Dec of the first
  delta; p[2] is the flux density of the first delta; and ``p[3]`` is
  the flux density of the second delta. Notice that the RA and Dec
  position of the second delta is equal to that of the first delta
  plus a fixed shift.

- **EXAMPLE 2**: A ring plus a delta at its center. The absolute
  position of the compound source is free, and the ring is circular.

  >>> model = ['ring', 'delta']
  >>> var = ['p[0], p[1], p[2], p[3], 1.0, 0.0', 'p[0], p[1], p[4]']

  In this case, ``p[0]`` and ``p[1]`` are the RA and Dec of both
  components (i.e., the ring and the delta); ``p[2]`` is the total
  flux density of the ring and ``p[3]`` is its diameter; ``p[4]`` is
  the flux density of the delta. Notice that the axes Ratio of the
  ring is set constant (and unity) and the position angle is also set
  constant (although it's meaningless for ``Ratio=1``).  For extended
  models, it is a good idea to bound the size to have positive values.
  Hence, in this case:

  >>> bounds = [None, None, None, [0, None]]

  In case we also want to bound the fitted fluxes to be positive, we
  would have:

  >>> bounds = [None, None, [0, None], [0, None], [0, None]]

- **EXAMPLE 3**: Like Example 1, but fixing also the flux-density of
  the second delta to be 2.5 times that of the first delta:

  >>> model = ['delta', 'delta']
  >>> var = ['p[0], p[1], p[2]', 'p[0]+0.5, p[3]+0.6, p[2]*2.5']

- **EXAMPLE 4**: A circularly-symmetric disc with a hole (i.e., with
  its inner half subtracted):

  >>> model = ['disc', 'disc']
  >>> var = ['p[0], p[1], 4/3*p[2], p[3], 1, 0',
  >>>        'p[0], p[1], -p[2]/3, p[3]/2, 1, 0']

  In this case, ``p[0]`` and ``p[1]`` are the RA and Dec shifts,
  respectively; ``p[2]`` is the flux density of the disc with the hole
  (i.e., with its inner half subtracted), and ``p[3]`` is the disc
  diameter. Notice that the hole in the disc has been created by
  adding a *negative* disc of size equals to 1/2 of the size of the
  larger disc, and flux density equals to -1/4 of that of the larger
  disc. The overall effect of both discs is that of one single disc
  with a hole (i.e., with no emission at radii < 1/2 of the outer
  radius).

- **EXAMPLE 5**: A delta component with a spectral index, fitted to
  the whole dataset (i.e., setting ``OneFitPerChannel=False``):

  >>> model = ['delta']
  >>> var = ['p[0], p[1], p[2]*(nu/1.e9)**p[3]']

  In this case, the fitted spectral index will be ``p[3]``, and
  ``p[2]`` will be the flux density at 1 GHz. Notice that p_ini (the
  list of initial values for the parameters) must also have the a
  priori value of the spectral index. For this example, p_ini could be
  (for a source close to the center, with an approximate flux of 2.3
  Jy at 1 GHz and an a priori spectral index of -0.7):

  >>> p_ini = [0.0, 0.0, 2.3, -0.7]

  If the spectral index is well known, it can of course be fixed in
  the fit. In such a case:

  >>> model = ['delta']
  >>> var = ['p[0], p[1], p[2]*(nu/1.e9)**(-0.7)']
  >>> p_ini = [0.0, 0.0, 2.3]

  .. note:: Fitting sizes *and* spectral indices at the same time may
      produce crazy results, since both quantities are quite coupled
      in Fourier space. Hence, some care should be taken when fitting
      to all frequencies together. Notice also that, unless your
      bandwidth is quite wide and/or your SNR is quite high, any fit
      to the spectral index may not be reliable.

- **EXAMPLE 6**: A filled sphere with its inner half (i.e., the core)
  removed. The sphere is fixed at the image center:

  >>> model = ['sphere', 'sphere']
  >>> var = ['0, 0, 9/8*p[0], p[1], 1, 0', '0, 0, -p[0]/8, p[1]/2, 1, 0']

  In this case, ``p[0]`` is the total flux density and ``p[1]`` is the
  outer diameter. This example is similar to that of the disc with a
  hole, but in this case the model is a sphere, where we have removed
  all the emission at radii < 0.5 times the outer radius.

- **EXAMPLE 7**: A disc with a hole of variable size. In this case,
  the smaller disc with negative flux density (which is used to
  generate the hole) must *always* have the same surface brightness
  (in absolute value) that the larger positive disc. Hence, and if we
  fix the disc position at the image center (for simplicity), we have:

  >>> model = ['disc', 'disc']
  >>> var = ['0, 0, p[0]*p[1]**2./(p[1]**2.-1), p[2], 1, 0',
  >>>        '0, 0, -p[0]/(p[1]**2.-1), p[2]/p[1], 1, 0']

  In this case, the flux density of the disc with the hole is ``p[0]``,
  the outer dimaeter is ``p[2]``, and ``p[1]`` is the ratio of the outer
  size to the inner size.

- **EXAMPLE 8**: Simultaneous self-calibration and source fitting. We
  fit a source flux density and perform Global Fringe Fitting (i.e.,
  we fit for phases, delays, and delay rates) in one shot. Let's fit
  stokes RR and assume that there are 3 antennas (for simplicity):

  >>> stokes = 'RR'
  >>> model = ['delta']
  >>> phase_gains = {1:'p[0] + 6.2832*(p[1]*(nu-nu0) + p[2]*t)',
  >>>                2:'p[3] + 6.2832*(p[4]*(nu-nu0) + p[5]*t)'}
  >>> var = '0, 0, p[6]'

  Good initial estimates for the phase gains can be obtained using the
  "QuinnFF" method.


=============
Model Details
=============

Each model depends on a list of variables that can be written as *any*
algebraic combination of fitting parameters. Some explanatory examples
are given in the section **Model Examples** (see also ``help(uvm)``,
if you imported the module with the ``uvm`` alias).

Here is the list of currently-supported models, and their variables:

- ``delta``    -> Variables: RA, Dec, Flux

- ``disc``     -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle

- ``ring``     -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle

- ``Gaussian`` -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle

- ``sphere``   -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle

- ``bubble``   -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle

- ``expo``     -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle

- ``power-2``  -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle

- ``power-3``  -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle

- ``GaussianRing``  -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle, Sigma

These models have the following meaning:

* ``sphere`` stands for an optically-thin uniform filled sphere.

* ``bubble`` stands for a uniform spherical surface.

* ``expo`` stands for an exponential radial flux decay.

* ``power-2`` stands for a decay as :math:`1/(r^2 + r_0^2)` (notice
  that in this case, the flux is the integral from ``r=0`` to ``r=r0``)

* ``power-3`` stands for a decay as :math:`1/(1 + (2^{2/3} - 1)(r/r_0)^2)^{3/2}`.

* ``GaussianRing`` stands for a radial profile following the expression:
  :math:`F*exp(-(Q - R)^2/(2*Sigma^2))`

  .. note:: The GaussianRing is built by a series expansion of the
    Bessel function in the computation of the Hankel transform. The
    order of the expansion (default is 80) can be set in the new
    HankelOrder keyword. The larger the source (in beam units) the
    higher HankelOrder should be, to keep the numerical precision.

The meaning of the model variables is:

- *RA* and *Dec* are the shifts w.r.t. the phase center (in arcsec)

- *Flux* is the total flux density of the component (in Jy)

- *Major* is the diameter along the major axis

- *Ratio* is the size ratio between the reference axis and the other
  axes (i.e., it is set to 1.0 for circularly-symmetric sources).

  .. warning:: If *Ratio* is forced to be higher than 1.0, then
     *Major* (see above) will indeed become the **minor** axis! If you
     want to avoid this issue, you should bound *Ratio* to be lower
     than (or equal to) 1.0. Otherwise, you will have to remember this
     issue when you interprete your fitted parameters.

- *PositionAngle* is the angle of the reference axis, from North to East (in deg.)

- *Sigma* (whenever it is used) is an auxiliary variable for models
  that need more than one 'size-like' parameter to be defined (e.g.,
  the 'GaussianRing', where we have the size of the ring and its
  width).

==========================
UVMultiFit Data Properties
==========================

*UVMultiFit* will return an instance with many properties and
methods. Some examples are listed below.


- ``averdata``: A set of arrays (one per spectral window) with the
  visibilities, averaged in time and/or frequency (according to the
  ``chanwidth`` and ``timewidth`` parameters).

- ``averweights``: The (square root of the) weights used in the fit.

- ``u``, ``v``, and ``w``: the coordinates in Fourier space of the
  visibilities (one array per spectral window).

- ``t``: the observing times (Modified Julian Dates). One array per spw.

.. note:: A UVMultiFit instance can occupy quite a bit of physical
     memory. It may be a good idea to delete the instance from memory
     (i.e., using ``del myfit`` and ``gc.collect()``) or restart CASA
     once the user has got the desired results from the fit (this is
     the recommended approach).

The most important methods of *UVMultiFit* are described below.

"""

__license__ = 'GPL-v3'
__revision__ = " $Id: 3.0.0-p2 2019-03-28 15:00:00 marti-vidal $ "
__docformat__ = 'reStructuredText'

import shutil
import re
import gc
import os
import time
import sys

import numpy as np
import scipy.special as spec

from casatools import ms
from casatools import table
from casatools import coordsys as cs
from casatools import image as ia

import uvmultimodel as uvmod
# global ms, tb

__version__ = "3.0-p2"
date = 'JAN 2021'

################
# Import all necessary modules.

print("\nC++ shared library loaded successfully\n")


tb = table()
ms = ms()

# if True:
#     from taskinit import gentools
#     from clearcal_cli import clearcal_cli as clearcal
#     ms = gentools(['ms'])[0]
#     tb = gentools(['tb'])[0]
#     ia = gentools(['ia'])[0]
#     cs = gentools(['cs'])[0]

greetings = '\n ##########################################################################\n'
greetings += ' # UVMULTIFIT --  ' + date + '. EUROPEAN ALMA REGIONAL CENTER (NORDIC NODE).  #\n'
greetings += ' #       Please, add the UVMULTIFIT reference to your publications:       #\n'
greetings += ' #                                                                        #\n'
greetings += ' #      Marti-Vidal, Vlemmings, Muller, & Casey 2014, A&A, 563, A136      #\n'
greetings += ' #                                                                        #\n'
greetings += ' ##########################################################################\n\n'


class uvmultifit():
    """This is the main UVMultiFit class.

   It includes pointers to the data and many different methods
   to perform fits, retrieve Chi square values, covariance matrices, calibrate data, write either model
   or residual visibilities to your measurement sets, etc.

   Fits with different models can be performed on the same data, (or on different time ranges within the
   same data to, e.g., perform variability studies) with no need to re-read the same data every time.

   The multi-threaded C++ engine, coupled to the flexibility of Python and CASA, makes UVMultiFit a fast,
   flexible and powerful tool for your advanced visibility-based analysis.

   :Parameters:
   ------------
   **vis** : `str`
     Name of the measurement set. It can also be a list of MS names.
   **spw** : `str`
     String (or list of strings). These are the spectral window(s) and channel selection(s)
     to be fitted. For each spw, the user can select many channel ranges in the usual CASA
     way. If one string is given, all ms listed in 'vis' are supposed to have the same
     spectral configuration. If a list of strings is given, one per ms, this restriction
     doesn't apply.
   **column** : `str`
     The data column. It can be either 'data' or 'corrected'
   **field** : `str`
     The id number (or name) of the target source in the measurement set(s).
   **pbeam** : `bool`
     If false, the primary-beam correction is not applied. This is *very important* for
     fitting mosaic data.
   **dish_diameter** : `double` or `dict`
     In case that the antenna diameters cannot be read from the datasets, the user must
     provide the antenna diameters (in meters). This can be given as a single float (so
     the array is assumed to be homogeneous) or as a dictionary, whose keys are
     antenna names (or *regular expressions*, that match antenna-name patterns) and whose
     elements are the antenna diameters (in meters).

     .. note:: For PB correction of heterogeneous, please use either one concatenated
               Measurement Set (MS) for all your data or several MSs with *the same*
               antenna table.

   **ldfac** : `double`
     Proportionality constant between the ratio 'lambda/Diameter' and the FWHM of the
     primary beam (assumed to be a Gaussian!). I.e.: FWHM = ldfac*lambda/Diameter.
     Normally, ``ldfac = 1.22`` should be fine, although 1.00 works better with data coming
     from simobserve.
   **scans** : `list`
     List of integers; default []. The id numbers of the scans to load. Default means to
     load all the scans of the selected source. If multiple measurement sets are selected,
     this should be a list of lists (i.e., one list of scans per measurement set).
     For instance, if ``vis = ["ms1.ms","ms2.ms"]``, then ``scans = [[1, 2, 3],[]]`` would select the
     scans 1, 2, and 3 from the ``ms1.ms`` dataset and all the scans from the ``ms2.ms``
     dataset.
   **chanwidth** : `int`
     Number of spectral channels to average in each chunk of data **before the fit**. BEWARE of
     the bandwidth-smearing effects if your baselines are long (and/or your source is large).
   **timewidth** : `int`
     Number of time channels (i.e., integration times) to average in each chunk of data
     **before the fit**. The averaging is always cut at the scan boundaries (so a very large
     timewidth means to average down to one visibility per scan). BEWARE of the
     time-smearing effects!
   **MJDrange** : `list`
     List of two floats. These are the initial and final Modified Julian Dates of the data
     to be used in the fitting. Notice that all the scans asked by the user are loaded
     a-priori, and the MJDrange condition is applied *afterwards*. This way, variability
     studies can be performed efficiently, by loading all data at once and then setting
     different MJDranges iteratively. Default (i.e., <= 0.0) means not to select data based
     on JD time range.
   **stokes** : `str`
     Polarization product. Can be any of ``PI, I, Q, U, V``. It also accepts individual
     correlation products: ``XX, YY, XY, YX, RR, LL, LR``, or ``LR``. Default is ``I``.
     If ``PI`` (which stands for *Polarization Independent*) is given, the program will
     compute ``I`` whenever possible and use either ``XX, YY, RR`` or ``LL`` otherwise.
     This way, the amount of data used in polarization-independent fits is maximized.
   **model** : `list`
     List of strings (i.e., model components to fit). Each component is given as a string.
     Possible models are: ``delta, disc, Gaussian, ring, sphere, bubble, expo, power-2,
     power-3``, and ``GaussianRing``.
     If only one model component is being fitted, the **model** keyword can also be a string.
   **var** : `list`
     List of strings (or just one string, if only one model component is being fitted).
     These are the variables for each model. The variables can be set to *any* algebraic
     expression involving the fitting parameters (being the ith parameter represented by
     ``p[i]``) and the observing frequency in Hz (represented by ``nu``). Any numpy function
     can also be called, using the prefix ``np`` (e.g., ``p[0]*np.power(nu, p[1])``).
     See some examples below.
   **p_ini** : `list`
     List of the initial values of the fitting parameters. This is expected to be a list
     of floats.
   **phase_center** : `str`
     The sky position where all components are referenced to. If an empty string is given,
     the phase center will be that of the first field id that is being read (i.e., if a
     mosaic is being read, the first pointing will be set as the phase center).
     If the string is not empty, the program expects to find a coordinate in CASA format
     (i.e., ``J2000 RA Dec``, where **RA** is in format ``00h00m00.0s`` and **Dec** is in
     format ``00d00m00.0s``).
   **fixed** : `list`
     Like **model**, but defines model components with completely fixed variables (i.e.,
     whose variables are defined only by numbers; not fitting parameters). This model will
     be computed only once (i.e., just before the fit), hence making the code execution
     faster. The user can load the model column of the measurement set(s) as a fixed model,
     by setting ``fixed='model_column'``.
   **fixedvar** : `list`
     Like **var**, but refers to the **fixed** model. Hence, it is expected to be either a
     list of numbers or a list of strings representing numbers. This is not needed if
     ``fixed = 'model_column'`` (since, in that case, the model column is read *as is* from
     the measurement set).
   **scalefix** : `str`
     String representing a function that sets a scaling factor for the fixed-model's total
     flux density. It *can be* a function of the fitting parameters (e.g., ``scalefix='p[0]'``
     will multiply the overall flux density of the fixed model by ``p[0]``) and can also be
     a function of the observing frequency ``nu``.
   **OneFitPerChannel** : `bool`
     If True, independent fits are performed to the different frequency channels, one by
     one. If False, one common fit is performed to all data. In this case, the user may
     want to fit for the spectral indices of the components, if the fractional bandwidth
     is wide.
   **outfile** : `str`
     Name of the output file to store results (i.e., fitting parameters, uncertainties,
     and metadata, in ascii format).
   **bounds** : `list`
     List of boundaries (i.e., minimum and maximum allowed values) for the fitting
     parameters. 'None' means that no bound is defined. If the list is empty, no bounds
     are assumed for any parameter (see examples below).
   **cov_return** : `bool`
     If True, the covariance matrix for each fit is added to the dictionary returning from
     the ``fit()`` method (the dictionary key will have the name ``'covariance'``, see
     the **Returns** section below).
   **uvtaper** : `double`
     Default is 0.0. If not 0.0, the weights of the visibilities are multiplied by a
     Gaussian in Fourier space, whose HWHM is the value of **uvtaper**, *in meters*.
   **uniform** : `bool`
     Default is False. If True, the weights of all data are made equal. Notice that a
     uvtaper can still be applied (i.e., by setting the **uvtaper** parameter).
   **finetune** : `bool`
     Default is False. If set to True, the fit is not performed, but only a ``uvmultifit``
     instance is created with the data properly read and the models compiled and ready.
     The user can then run different methods of the UVMultiFit class by him/herself (see
     the help text for each method) before actually fitting the data. This can be useful,
     for instanse, if the user wants to try many different models (and/or subtract many
     different fixed models), apply ad-hoc calibrations, perform an *MCMC* exploration of
     the parameter space, etc., without having to reload the data every time after each
     step.
   **wgt_power** : `double`
     Default is 1. Power index of the visibility weights in the computation of the Chi
     square. ``wgt_power = 1`` would be the *statistically justified* value, although other
     values may be tested if the user suspects that the visibility uncertainties are not
     well estimated in his/her dataset.
   **method** : `str`
     Method to use in the chi-square minimization. Default is ``simplex``. Possible values
     are ``simplex`` and ``levenberg``. Sometimes, the least-squares minimization may not
     converge well with the *Levenberg-Marquardt* method (if the model is complicated and/or
     the uv coverage is quite sparse). *Levenberg-Marquardt* also requires a lot of memory,
     so may not be very convenient if the datasets are very large. In these cases, a Chi
     square minimization using the *simplex* algorithm may work beter. However, *simplex*
     may require more function evaluations to find the minimum.

     .. warning:: The SIMPLEX method does not currently provide parameter uncertainties.
                  It just returns the reduced Chi squared!

   **write** : `str`
     The kind of information to store in the measurement set(s) after the fit.
     Default is '' (i.e., does not change anything in the datasets).

   * If it is set to ``model``, the best-fit model is saved in the *model column* of the
     measurement set.

   * If it is set to ``residuals``, the fit residuals are saved in the *corrected' column*
     of the measurement set.

   * If it is set to ``calibrated`` (this option will be available in the next release),
     the gains defined in the **amp_gains** and **phase_gains** dictionaries (see below)
     will be applied, and the calibrated data will be saved in the *corrected column*
     of the ms.

     Currently, this keyword is only activated if **stokes** is set to either ``PI, I`` or
     an individual correlation product (like ``XX`` or ``XY``) *and* if both **timewidth**
     and **chanwidth** are set to 1.
   **NCPU** : `int`
     Default is 4. Number of threads allowed to run in parallel.
   **SMPtune** : `list`
     Used to fine-tune the Simplex algorithm. Change only if you really know what you
     are doing. The meaning of the list elements is:

   * ``SMPtune[0]`` -> Maximum allowed error in the parameters, from the search of the Chi
     Square minimum. Default is 1.e-4.

   * ``SMPtune[1]`` -> Relative size of the first step in the parameter search.
     Default is 1.e-1.

   * ``SMPtune[2]`` -> Maximum number of iterations allowed per fitting parameter.
     Default is 200
   **LMtune** : `list`
     Used to fine-tune the Levenberg-Marquardt algorithm. Change only if you really know
     what you are doing. The meaning of the list elements is:

   * ``LMtune[0]`` -> The Lambda factor to weight up the Newton component
     (i.e., *H + Lambda*H_diag = Res*), where *H* is the Hessian and *Res* the vector of
     residuals. Default is 1.e-3

   * ``LMtune[1]`` -> The factor to multiply (divide) Lambda if the iterator worsened
     (improved) the Chi Squared. Default: 10.

   * ``LMtune[2]`` -> The maximum relative error allowed for the ChiSq. Default: 1.e-5

   * ``LMtune[3]`` -> Maximum number of iterations allowed per fitting parameter.
     Default: 200

   * ``LMtune[4]`` -> (optional) Maximum relative error allowed for the parameters.
     Default: 1.e-3. If it is not provided, then ``LMtune[4] = LMtune[2]``
   **only_flux** : `bool`
     Default is False. If True, the program assumes that only the flux densities of all
     components are going to be fitted. Furthermore, the fitting parameters shall be just
     equal to the flux densities being fitted, and should be given in the same order as
     the list of model components. In these cases, using ``only_flux = True`` can speed
     up the fit quite a bit, especially if there are many components to be fitted. This
     option may be useful to implement simple *compressed sensing*-like algorithms.
   **proper_motion** : `list`
     List of 2-element lists of numbers: Each element (i.e., each list of two numbers)
     is the proper motion, in RA and Dec, of each model component. The units are
     arc-seconds per year. Proper motions cannot be fitted yet, but may be fittable in
     future versions of UVMultiFit. Default is a float = 0.0, meaning that all proper
     motions are null. If the proper motions are not null, the position used as reference
     for the fit will correspond to that of the first integration time of the first scan
     of the observed field.
   **HankelOrder** : `int`
     Only used for models without an analytic expression (i.e., at the moment, only the
     ``GaussianRing`` model). In these cases, UVMultiFit performs the Hankel transform by
     using the series expansion of the Bessel function J0. The order of this expansion is
     set by this keyword. The larger the distance in Fourier space (i.e., the larger the
     source, in beam units), the higher HankelOrder should be, to keep the numerical
     precision. Default is 80. For cases of sources with sizes similar to the synthesized
     beam, HankelOrder=80 should suffice. The use of too high values may cause Overflow
     errors!
   **amp_gains**: `dict`
     Dictionary (default empty). Controls whether to solve for antenna amplitude gains
     and source-model parameters *simultaneously*. The keys of the dictionary are the
     indices of the antennas to self-calibrate. The values of the dictionary are strings,
     whose elements represent functions of the fitting parameters. See examples below.
   **phase_gains** : `dict`
     Same as for amp_gains, but dealing with phase gains.

     .. note:: For the gain solver to work, please use either one concatenated Measurement
               Set (MS) for all your data or ensure that all your MSs have the same
               antenna table.

    """

    # global sys, ms, uvmod, goodclib, time, np, sp, spopt, spec, os, gc, gentools, re, clearcal, ms, tb, cs

    ############################################
    #
    #  FREE MEMORY
    #
    def _deleteData(self, delmodel=True):
        """ Delete pointers to the data.
        Hopefully, this will release memory when gc.collect() is run."""

        for i in range(len(self.averdata) - 1, -1, -1):
            del self.averdata[i]
            if self.takeModel:
                try:
                    del self.avermod[i]
                except:
                    pass
            del self.averweights[i]
            del self.averfreqs[i]
            del self.u[i]
            del self.v[i]
            del self.w[i]
            del self.t[i]
            del self.tArr[i]
            del self.tIdx[i]
            del self.ant1[i]
            del self.ant2[i]
            del self.RAshift[i]
            del self.Decshift[i]
            del self.Stretch[i]

        try:
            for mdi in self.iscancoords[::-1]:
                Npar = len(mdi)
                for mdp in range(Npar - 1, -1, -1):
                    del mdi[mdp]
                del mdi
            del self.iscancoords
        except:
            pass

        del self.averdata, self.avermod, self.averweights, self.averfreqs, self.v, self.u, self.w, self.t
        del self.ant1, self.ant2, self.RAshift, self.Decshift, self.Stretch

        if delmodel:
            self.mymodel._deleteData()

    def __del__(self):
        # Delete the model first, so that C++ has green light to free pointers:
        del self.mymodel

        # Clear all data and C++ pointers:
        self._deleteData(delmodel=False)
        uvmod.clearPointers(0)
        uvmod.clearPointers(1)

    ############################################
    #
    #  CREATE INSTANCE
    #
    def __init__(self, vis='', spw='0', column='data', field=0, scans=[], uniform=False,
                 chanwidth=1, timewidth=1, stokes='I', write='', MJDrange=[-1.0, -1.0], ldfac=1.22,
                 model=['delta'], var=['p[0], p[1], p[2]'], p_ini=[0.0, 0.0, 1.0], phase_center='',
                 fixed=[], fixedvar=[], scalefix='1.0', outfile='modelfit.dat', NCPU=4, pbeam=False,
                 dish_diameter=0.0, gridpix=0,
                 OneFitPerChannel=True, bounds=None, cov_return=False, finetune=False, uvtaper=0.0,
                 method='levenberg', wgt_power=1.0,
                 LMtune=[1.e-3, 10., 1.e-5, 200, 1.e-3], SMPtune=[1.e-4, 1.e-1, 200], only_flux=False,
                 proper_motion=0.0, HankelOrder=80, phase_gains={}, amp_gains={}):
        """ Just the constructor method, for class instantiation."""
        self.first_time = True
        uvmod.clearPointers(2)

        self.implemented = ['delta', 'disc', 'ring', 'Gaussian', 'sphere',
                            'bubble', 'expo', 'power-2', 'power-3', 'GaussianRing']
        # Number of variables for the models
        numvars = [3, 6, 6, 6, 6, 6, 6, 6, 6, 7]
        # Models that need a series expansion for J0.
        self.isNumerical = ['GaussianRing']
        self.uniform = uniform
        self.vis = vis
        self.column = column
        self.averdata = []
        self.avermod = []
        self.averweights = []
        self.averfreqs = []
        self.u = []
        self.v = []
        self.model = model
        self.var = var
        self.fixed = fixed
        self.wgt_power = wgt_power
        self.fixedvar = fixedvar
        self.scalefix = scalefix
        self.bounds = bounds
        self.outfile = outfile
        self.ranges = []
        self.OneFitPerChannel = OneFitPerChannel
        self.spw = spw
        self.field = field
        self.chanwidth = chanwidth
        #  self.timewidth = timewidth
        self.stokes = stokes
        self.p_ini = p_ini
        self.nspws = []
        self.cov_return = cov_return
        self.uvtaper = uvtaper
        self.times = []
        self.method = method
        self.field_id = []
        self.pointing = []
        self.sourscans = []
        self.LtSpeed = 2.99792458e+8
        self.FouFac = (2. * np.pi) * np.pi / 180. / 3600.
        self.variables = []
        self.NCPU = NCPU
        self.MJDrange = MJDrange
        self.scans = scans
        self.finetune = finetune
        self.minDp = np.sqrt(np.finfo(np.float64).eps)
        self.LMtune = LMtune
        self.SMPtune = SMPtune
        self.only_flux = only_flux
        self.proper_motion = proper_motion
        self.phase_center = phase_center
        self.HankelOrder = HankelOrder
        self.amp_gains = amp_gains
        self.phase_gains = phase_gains
        self.ldfac = ldfac

        # Some variables to tune the mosaic fitting:
        self.pbeam = pbeam
        self.dish_diameter = dish_diameter

        # Number of variables for each model
        self.numvar = {}
        for i, component in enumerate(self.implemented):
            self.numvar[component] = numvars[i]
        self.maxNvar = np.max(numvars)

        try:
            self.write_model = {'': 0, 'model': 1, 'residuals': 2, 'calibrated': 3}[write]
        except IndexError:
            self._printError("ERROR: keyword 'write' should be set to either '', 'model', 'residuals' or 'calibrated'")

        #######################
        # TIMEWIDTH IS NOT CURRENTLY USED (ISSUES WITH NEW MS TOOL):
        self.timewidth = 1
        if timewidth != 1:
            self._printInfo("WARNING! timewdith>1 cannot be currently set, due to issues with new MS tool")

        # Start instance:
        self._startUp()

    ############################################
    #
    #  STARTERS
    #
    # This method will be overriden in the GUI, to avoid execution of CheckInputs() and the fit:
    def _startUp(self):
        """ This is run each time the class is instantiated."""
#         if not goodclib:
#             self._printError("ERROR: C++ library cannot be loaded! Please, contact the Nordic ARC node.")
#             return False

        self._printInfo(greetings)

        self.mymodel = modeler()
        self.mymodel.Ccompmodel = uvmod.modelcomp

        # Check parameters and read the data in:
        goodread = self.checkInputs()
        if goodread:
            goodread = self.readData(del_data=False)
        else:
            self._printError("\n\n ABORTING UVMULTIFIT! BAD INPUTS! CHECK INPUTS!\n")

            # Compile Model:
        if goodread:
            goodread = self.initModel()
        else:
            self._printError("\n\n ABORTING UVMULTIFIT! BAD DATA! CHECK INPUTS!\n")

        # Fit if finetune == False:
        if goodread:
            if not self.finetune:
                goodfit = self.fit()
                if not goodfit:
                    self._printError("\n FAILED FIT!\n")
                elif self.write_model in [1, 2, 3]:
                    if self.timewidth == 1 and self.stokes not in ['Q', 'U', 'V']:
                        self._printInfo("\nWriting into measurement set(s)\n")
                        self.writeModel()
                    else:
                        msg = "\nCANNOT write into measurement set!\n"
                        msg += "\nIf you want to fill-in the model (or corrected) column:\n"
                        msg += "    1.- 'timewidth' and 'chanwidth' should be set to 1\n"
                        msg += "    2.- 'stokes' should be set to either 'I', 'PI', or a corr. product."
                        self._printError(msg)

                if goodfit:
                    self._printInfo("\n\n\nFit done!! And UVMULTIFIT class instantiated successfully!\n")
            else:
                self._printInfo("\n\n\nUVMULTIFIT class instantiated successfully!\n")

        else:
            self._printError("\n\n ABORTING UVMULTIFIT! BAD MODEL! CHECK INPUTS!\n")

    ############################################
    #
    #  PRINT MESSAGES AND ERRORS
    #
    # Functions overriden in GUI mode:
    def _printError(self, message):
        """ Prints a message and raises an exception."""
        sys.stdout.write(message)
        sys.stdout.flush()
        raise Exception(message)

    def _printInfo(self, message):
        """ Prints a message."""
        sys.stdout.write(message)
        sys.stdout.flush()

    ##################################
    # FRINGE FITTER:
    #
    def QuinnFF(self, IF, refant, doModel, doGlobal):
        """Perform a *Global Fringe Fitting* search in delay-rate space.

      It uses the Quinn estimator for the fringe peak.

      :Parameters:
      ------------
      **IF** : `int`
          The spectral window to fringe-fit (index referred to the data already read).
      **refant** : `int`
          The index of the antenna to use as reference. If this antenna is missing or has
          bad data, another refant will be picked automatically.
      **doModel** : `bool`
          If True, deconvolve the fixed model before the GFF.
      **doGlobal** : `bool`
          If True, globalize the solutions (i.e., use all available baselines to estimate
          the antenna gains).

      :Returns:
      ---------
      Returns a `list` of 5 elements:

      * An array with the delays (in seconds). One value per antenna (the reference
        antenna, as well as any other antenna with failed solutions, will have null
        values).

      * An array with the rates (in 1/s), one value per antenna.

      * An array with the phases (in radians), one value per antenna.

      * A double (the delay equivalent to one pixel in the delay/rate matrices).

      * A double (the rate equivalent to one pixel in the delay/rate matrices).

      """
        return uvmod.QuinnFF(IF, refant, doModel, doGlobal)

    ############################################
    #
    #  WRITE MODEL (OR RESIDUALS)
    #
    def writeModel(self):
        """ Writes the requested information into the measurement sets.

        The information can be either the predictions of the compiled model(s) (i.e., they are
        written into the *model column* of the measurement set(s)), or the post-fit residuals
        (they are written into the *corrected column*) or the calibrated data (they are written
        into the *corrected column* as well). The actual information to write is set by the value
        of the ``write`` keyword of the *UVMultiFit* instance when the ``fit()`` method was called.

        This function is executed only if the ``stokes`` keyword is set to either ``PI``, ``I``
        or an individual correlation product (like, e.g., ``XX`` or ``XY``) *and* if no averaging
        has been performed neither in time nor frequency (i.e., if both ``timewidth`` and
        ``chanwidth`` are set to 1). This function should be called AFTER having run ``fit()``.
        """

        self._printInfo("\n\n WARNING: WRITING TO MOSAICS IS EXPERIMENTAL AND MAY NOT WORK!\n")

        for vi, v in enumerate(self.vis):
            # Get the columns of parallel-hand correlations:
            success = ms.open(v)
            if not success:
                self._printError("ERROR! %s cannot be openned in write mode" % (v))
                return False

            spws = list(map(int, self.iscan[v].keys()))

            ms.selectinit(datadescid=spws[0])
            polprods = [x[0] for x in list(ms.range(['corr_names'])['corr_names'])]
            ms.close()

            if self.stokes in polprods:
                polii = [polprods.index(self.stokes)]
            elif self.stokes in ['PI', 'I']:
                if 'XX' in polprods:
                    polii = [polprods.index('XX'), polprods.index('YY')]
                elif 'RR' in polprods:
                    polii = [polprods.index('RR'), polprods.index('LL')]
                else:
                    self._printError("Stokes not understood for %s! Will not update the model column!" % (v))
                    # ms.close()
                    return False
            else:
                self._printError("Stokes not understood for %s! Will not update the model column!" % (v))
                # ms.close()
                return False

            if self.write_model == 1:
                column = 'MODEL_DATA'
            elif self.write_model in [2, 3]:
                column = 'CORRECTED_DATA'

            # NEW CODE TO WRITE MODEL, BASED ON TB TOOL:
            for sp in spws:
                for scan in self.iscan[v][sp].keys():
                    for field in range(len(self.iscan[v][sp][scan])):
                        sys.stdout.write("\r Doing %s: spw %i, scan_id %i" % (v, sp, scan))
                        sys.stdout.flush()
                        tb.open(v, nomodify=False)
                        select = self.iscan[v][sp][scan][field]
                        tb2 = tb.selectrows(select[-1])
                        moddata = tb2.getcol(column)
                        re = np.transpose(self.mymodel.output[select[0]][select[1]:select[1] + select[2], :])
                        for nui, r in enumerate(select[3]):
                            for poli in polii:
                                moddata[poli, r, :] = re[nui, :]
                        tb2.putcol(column, moddata)
                        tb.close()

            self._printInfo("\n %s written successfully!\n" % column)

    ############################################
    #
    #  SANITY CHECKS AND DEFINING BASIC PARAMETERS
    #
    def _checkOrdinaryInputs(self):
        """ Performs some sanity checks on the input parameters.
      This function should not be called directly by the user."""

        if isinstance(self.model, str):
            self.model = list([self.model])
        else:
            self.model = list(self.model)

        if isinstance(self.var, str):
            self.var = list([self.var])
        else:
            self.var = list(self.var)

        if type(self.fixedvar) in [str, float]:
            self.fixedvar = list([self.fixedvar])
        else:
            self.fixedvar = list(self.fixedvar)

        if isinstance(self.fixed, str):
            self.fixed = list([self.fixed])
        else:
            self.fixed = list(self.fixed)

        if type(self.p_ini) in [float, str]:
            self.p_ini = list([self.p_ini])
        else:
            self.p_ini = list(self.p_ini)

        # Check if numerical approximations are needed:
        isNum = False
        for mods in self.isNumerical:
            if mods in self.model:
                isNum = True
        if isNum:
            self.applyHankel = self.HankelOrder
        else:
            self.applyHankel = 0

        self.isMixed = False
        self.phase_gainsNu = {}
        self.phase_gainsT = {}
        self.amp_gainsNu = {}
        self.amp_gainsT = {}
        self.phase_gainsNuT = {}
        self.amp_gainsNuT = {}

        if not isinstance(self.phase_gains, dict):
            self._printError("\n 'phase_gains' must be a dictionary!")

        if not isinstance(self.amp_gains, dict):
            self._printError("\n 'amp_gains' must be a dictionary!")

        for key in self.phase_gains.keys():
            if isinstance(key, int):
                self.isMixed = True
                break

        for key in self.amp_gains.keys():
            if isinstance(key, int):
                self.isMixed = True
                # self._printError("\n Inconsistent 'amp_gains' and 'phase_gains'!")
                break

        for key in self.phase_gains.keys():
            if key == 'nuG':
                if self.isMixed:
                    self._printError(
                        "You cannot define split gains and mixed gains in the same fit.\n It makes no sense!\n")
                self.phase_gainsNu = self.phase_gains['nuG']
            elif key == 'tG':
                if self.isMixed:
                    self._printError(
                        "You cannot define split gains and mixed gains in the same fit.\n It makes no sense!\n")
                self.phase_gainsT = self.phase_gains['tG']
            elif not isinstance(key, int):
                self._printError("\n The keys of 'phase_gains' must be integers or 'nuG/tG'!")
            else:
                self.phase_gainsNuT[key] = self.phase_gains[key]

        for key in self.phase_gainsNu.keys():
            if not isinstance(key, int):
                self._printError("\n The keys of 'phase_gains[nuG]' must be integers!")
        for key in self.phase_gainsT.keys():
            if not isinstance(key, int):
                self._printError("\n The keys of 'phase_gains[tG]' must be integers!")

        for key in self.amp_gains.keys():
            if key == 'nuG':
                if self.isMixed:
                    self._printError(
                        "You cannot define split gains and mixed gains in the same fit.\n It makes no sense!\n")
                self.amp_gainsNu = self.amp_gains['nuG']
            elif key == 'tG':
                if self.isMixed:
                    self._printError(
                        "You cannot define split gains and mixed gains in the same fit.\n It makes no sense!\n")
                self.amp_gainsT = self.amp_gains['tG']
            elif not isinstance(key, int):
                self._printError("\n The keys of 'amp_gains' must be integers or 'nuG/tG'!")
            else:
                self.amp_gainsNuT[key] = self.amp_gains[key]

        for key in self.amp_gainsNu.keys():
            if not isinstance(key, int):
                self._printError("\n The keys of 'amp_gains[nuG]' must be integers!")
        for key in self.amp_gainsT.keys():
            if not isinstance(key, int):
                self._printError("\n The keys of 'amp_gains[tG]' must be integers!")

        self.useGains = set(list(self.phase_gainsNuT) + list(self.amp_gainsNuT) +
                            list(self.phase_gainsNu) + list(self.amp_gainsNu) +
                            list(self.phase_gainsT) + list(self.amp_gainsT))

        # Check phase center:

        if not isinstance(self.phase_center, str):
            self._printError("\n 'phase_center' must be a string!")
        else:
            if len(self.phase_center) == 0:
                self.phrefset = False
            else:
                if True:
                    #  try:
                    self.phrefset = True
                    csys = cs.newcoordsys(direction=True)
                    dirstr = self.phase_center.split()
                    if len(dirstr) == 2:
                        csys.setdirection(refcode="J2000", refval=self.phase_center)
                    else:
                        csys.setdirection(refcode=dirstr[0], refval=" ".join(dirstr[1:]))
                    csys.convertdirection("J2000")
                    self.refpos = np.copy(csys.torecord()['direction0']['crval'])
                else:
                    #  except:
                    self._printError("\n 'phase_center' is not a CASA-formatted sky coordinate!")

        # Did the user forget how does this task work? :D
        for param in self.var:
            if not isinstance(param, str):
                self._printError("\n 'var' must be a list of strings!")
                return False

        # Get the number of parameters and check the model consistency:
        maxpar = 0
        # Typical function that can be used.
        lf = ['GaussLine', 'GaussLine', 'LorentzLine',
              'LorentzLine', 'power', 'maximum', 'minimum']
        for i, component in enumerate(self.model):
            vic = self.var[i].count
            checkpars = self.var[i].split(',')
            nfuncs = sum(list(map(vic, lf)))  # 2*(self.var[i].count('GaussLine')+self.var[i].count('LorentzLine'))
            paridx = zip([m.start() + 1 for m in re.finditer('\[', self.var[i])],
                         [m.start() for m in re.finditer('\]', self.var[i])])
            maxpar = max([maxpar] + list(map(float, [self.var[i][ss[0]:ss[1]] for ss in paridx])))
            if component not in self.implemented:
                msg = "\nModel component '" + str(component) + "' is not known!\n"
                fmt = "Supported models are:" + " '%s' " * len(self.implemented)
                msg += fmt % tuple(self.implemented)
                self._printError(msg)
                return False
            if (component in self.implemented and self.numvar[component] != (len(checkpars) - nfuncs) != 6):
                self._printError("Wrong number of variables (%i) in '%s' model." % (len(checkpars) - nfuncs, component))
                return False

        # Scalefix must be a string representing a function:
        if not isinstance(self.scalefix, str):
            self._printError("'scalefix' should be a string!")
            return False

        # Get the overal number of parameters (i.e., from the variables AND the scalefix string):
        paridx = zip([m.start() + 1 for m in re.finditer('\[', self.scalefix)],
                     [m.start() for m in re.finditer('\]', self.scalefix)])
        maxpar = max([maxpar] + list(map(float, [self.scalefix[ss[0]:ss[1]] for ss in paridx])))

        # Get the model parameters for the gains:
        #  self.maxGainTerm = 0
        #  self.NgainAnts = {}
        for key in self.phase_gainsNuT:
            term = self.phase_gainsNuT[key]
            paridx = zip([m.start() + 1 for m in re.finditer('\[', term)],
                         [m.start() for m in re.finditer('\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.amp_gainsNuT:
            term = self.amp_gainsNuT[key]
            paridx = zip([m.start() + 1 for m in re.finditer('\[', term)],
                         [m.start() for m in re.finditer('\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.phase_gainsNu:
            term = self.phase_gainsNu[key]
            paridx = zip([m.start() + 1 for m in re.finditer('\[', term)],
                         [m.start() for m in re.finditer('\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.amp_gainsNu:
            term = self.amp_gainsNu[key]
            paridx = zip([m.start() + 1 for m in re.finditer('\[', term)],
                         [m.start() for m in re.finditer('\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.phase_gainsT:
            term = self.phase_gainsT[key]
            paridx = zip([m.start() + 1 for m in re.finditer('\[', term)],
                         [m.start() for m in re.finditer('\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.amp_gainsT:
            term = self.amp_gainsT[key]
            paridx = zip([m.start() + 1 for m in re.finditer('\[', term)],
                         [m.start() for m in re.finditer('\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        # The fixed model should contain CONSTANT variables:
        for fcomp in self.fixedvar:
            for pi, param in enumerate(fcomp.split(',')):
                try:
                    temp = float(param)
                except ValueError:
                    self._printError("\n Fixed variables must be a list of floats (or strings\n representing floats)!")
                    return False

        # There should be as many p_inis as parameters!
        if len(self.p_ini) != maxpar + 1:
            self._printError("\n 'p_ini' is of length %i, but there are %i parameters used!\n ABORT!" %
                             (len(self.p_ini), maxpar + 1))
            return False

        # Check for consistency in the fixed model:
        if len(self.fixed) > 0:
            self.takeModel = 'model_column' in self.fixed
        else:
            self.takeModel = False

        if self.takeModel:
            self._printInfo("MODEL COLUMN will be taken as fixed model.\n" +
                            "Skipping the 'fixedvar' column and all other fixed components")
            self.fixed = ['model_column']
            self.fixedvar = []
        else:
            for i, component in enumerate(self.fixed):
                checkpars = self.fixedvar[i].split(',')
                if component not in self.implemented:
                    msg = "\nModel component '" + str(component) + "' is not known!\n"
                    fmt = "Supported components are " + " '%s' " * len(self.implemented)
                    msg += fmt % tuple(self.implemented)
                    self._printError(msg)
                    return False
                if (component in self.implemented and self.numvar[component] != len(checkpars) != 6):
                    self._printError("Wrong number of variables (%i) in '%s' fixed model." %
                                     (len(checkpars), component))
                    return False

        # Set outfile name:
        if self.outfile == '':
            self._printInfo("\nSetting 'outfile' to its default (i.e., 'modelfit.dat').\n")
            self.outfile = 'modelfit.dat'

        # Set (and check) the bounds in the fit:
        if (self.bounds is None) or (len(self.bounds) == 0):
            self.bounds = None
        else:
            for b, bound in enumerate(self.bounds):
                if bound is None:
                    self.bounds[b] = [None, None]
                if self.bounds[b][0] is not None and self.p_ini[b] <= self.bounds[b][0]:
                    self._printError(
                        "Initial value (%.2e) of parameter %i is lower (or equal) than its lower boundary (%.2e)!" %
                        (self.p_ini[b],
                         b, self.bounds[b][0]))
                    return False
                if self.bounds[b][1] is not None and self.p_ini[b] >= self.bounds[b][1]:
                    self._printError(
                        "Initial value (%.2e) of parameter %i is larger (or equal) than its upper boundary (%.2e)!" %
                        (self.p_ini[b],
                         b, self.bounds[b][1]))
                    return False

        if (self.bounds is not None) and (len(self.bounds) != len(self.p_ini)):
            self._printError("Length of 'bounds' list (%i) is not equal to number of parameters (%i)!" %
                             (len(self.bounds), len(self.p_ini)))
            return False
        return True

    ############################################
    #
    #  MORE SANITY CHECKS
    #
    def checkInputs(self):
        """ Reads all the inputs parameters and performs sanity checks.

        You have to run this *everytime* after changing any parameter of the UVMUltiFit
        instance.

        Many self-consistency checks are performed, related to data selection, number of
        internal variables reading, averaging, and formatting, etc. It is always a good
        idea to run this method before fitting, if the user has changed some parameters
        (to redo a fit) or is using ``finetune=True``.

        .. warning:: If you change any keyword related to the data selection (with the
        exception of ``MJDRange``) you **must** run ``readData()`` after ``checkInputs()``,
        to ensure that the new data are properly read and set before the fit.

        The keywords that need a rerun of ``readData()`` are:
        ``vis, spw, stokes, column, scan, field, uniform, timewidth, chanwidth, phase_center,
        wgt_power, uvtaper``.

        .. warning:: If you change the equations of your fitting model (ond/or those of
        the antenna gains) you **must** run ``initModel()`` before the fit.

        The keywords that need a rerun of ``initModel()`` are:
        ``model, var, fixed, fixedvar, scalefix``

        Notice that just changing ``p_ini`` or ``bounds`` is OK. You do not need to reinit
        the model to repeat the fit.
      """

        # Some preliminary (self-consistency) checks of the parameters:
        self.savemodel = True

        # As the function says, check ordinary inputs :)
        success = self._checkOrdinaryInputs()
        if not success:
            return False

        # We always work with lists here:
        if isinstance(self.vis, str):
            self.vis = list([self.vis])
        else:
            self.vis = list([str(ss) for ss in self.vis])

        # Set the right format for the spectral window(s):
        if isinstance(self.spw, str):
            self.spw = list([self.spw])
        elif len(self.spw) > 1:
            self.spw = list([str(ss) for ss in self.spw])
            if self.OneFitPerChannel:
                self._printInfo("\n\n SPW is a LIST! User BEWARE! Any fit in *spectral mode*\n" +
                                "WILL fit the model to each MS separately!\n\n")
        elif len(self.spw) > 0:
            self.spw = list([self.spw[0]])
        else:
            self._printError("\nBad formatted spw!\n")
            return False

        # Set list of scans:
        try:
            if isinstance(self.scans, list) and len(self.scans) == 0:
                self.scans = [[] for v in self.vis]
            # We work with lists of nvis elements:
            self.scans = list(self.scans)
            if isinstance(self.scans[0], int):
                if len(self.vis) == 1:
                    self.scans = [self.scans]
                else:
                    self._printError("\n 'scans' should be a list of integers (or a list of lists of integers,\n" +
                                     "if there are several measurement sets).\n ABORTING!")
                    return False
            if len(self.scans) != len(self.vis):
                self._printError("\n List of (lists of) scans does not have the same length " +
                                 "as the list of measurement sets!\n ABORTING!")
                return False

            for si, sc in enumerate(self.scans):
                if isinstance(sc, str):
                    self.scans[si] = list(map(int, sc.split(',')))
        except:
            self._printError("\n 'scans' should be a list of integers (or a list of lists of integers,\n" +
                             "if there are several measurement sets).\n ABORTING!")
            return False

        # Check dimensions of vis, spw, model, etc.:
        if len(self.vis) != len(self.spw) and len(self.spw) > 1:
            self._printError("\n\nThe length of 'spw' is not equal to the length of 'vis'!\n Aborting!\n")
            return False

        if not isinstance(self.stokes, str):
            self._printError("'stokes' must be a string!")
            return False

        if not isinstance(self.only_flux, bool):
            self._printError("'only_flux' must be a boolean!")

        if self.only_flux:
            if len(self.p_ini) != len(self.model):
                self._printError("If only_flux=True, number of parameters must be equal to " +
                                 "number of model components!\n Aborting!\n")

        if self.proper_motion == 0.0:
            self.proper_motion = [[0., 0.] for i in self.model]
        elif not isinstance(self.proper_motion, list):
            self._printError("'proper_motion' must be a list!")
        else:
            if len(self.proper_motion) != len(self.model):
                self._printError("The length of 'proper_motion' must be equal to the number of model components!")
            for pi in self.proper_motion:
                if not isinstance(pi, list):
                    self._printError("The elements of 'proper_motion' must be lists of two floats!")
                elif len(pi) != 2:
                    self._printError("The elements of 'proper_motion' must be lists of two floats!")
                elif not isinstance(pi[0], float) or not isinstance(pi[1], float):
                    self._printError("The elements of 'proper_motion' must be lists of two floats!")

        # Do the mss exist?
        for visi in self.vis:
            if not os.path.exists(visi):
                self._printError("\nMeasurement set %s does not exist!" % (visi))
                return False

        # Can the required column exist?
        if self.column not in ['data', 'corrected_data', 'corrected']:
            self._printError("\n'column' can only take values 'data' or 'corrected'!")
            return False
        if self.column == 'corrected':
            self.column = 'corrected_data'

        self.pointing = []
        self.sourscans = []
        phasedirs = [{} for v in self.vis]
        self.field_id = []

        # Get number of antennas:
        print(os.path.join(self.vis[0], 'ANTENNA'))
        print(tb)
        tb.open(os.path.join(self.vis[0], 'ANTENNA'))
        self.Nants = len(tb.getcol('NAME'))
        tb.close()

        # Open MS and look for the selected data:
        for vi, v in enumerate(self.vis):
            self.field_id.append([])
            success = ms.open(v)
            if not success:
                self._printError("Failed to open measurement set '+v+'!")
                return False

            allfields = list(ms.range('fields')['fields'])

            try:
                if type(self.field) in [str, int]:
                    fitest = int(self.field)
                else:
                    fitest = int(self.field[vi])
                phasedirs[vi][fitest] = ms.range('phase_dir')['phase_dir']['direction'][:, fitest]
                self.field_id[-1].append(fitest)
            except:
                if isinstance(self.field, str):
                    aux = str(self.field)
                    self.field = [aux for v in self.vis]

                for f, field in enumerate(allfields):
                    if self.field[vi] in field:
                        self.field_id[-1].append(f)
                        phasedirs[vi][f] = ms.range('phase_dir')['phase_dir']['direction'][:, f]
                if len(self.field_id[-1]) == 0:
                    self._printError("\nERROR! Field %s is not in %s" % (self.field[vi], v))
                    ms.close()
                    return False

            ms.close()
        self.phasedirs = phasedirs

        # Ref. position:
        if not self.phrefset:
            self._printInfo("\nSetting phase center on first scan\n")
            self.refpos = self.phasedirs[0][min(self.phasedirs[0].keys())]
        else:
            self._printInfo("\nSetting phase center on %s\n" % self.phase_center)

        # Find out all the scans where this field id is observed:
        for vi, v in enumerate(self.vis):
            success = ms.open(v)
            if not success:
                self._printError("Failed to open measurement set '+v+'!")
                return False

            info = ms.getscansummary()
            ms.close()

            self.sourscans.append([])
            self.pointing.append([])
            for key in info.keys():
                if info[key]['0']['FieldId'] in self.field_id[vi]:
                    self.sourscans[-1].append(int(key))
                for fieldid in info[key].keys():
                    myfield = info[key][fieldid]['FieldId']
                    if myfield in self.field_id[vi]:
                        fi = self.field_id[vi].index(myfield)
                        self.pointing[-1].append(phasedirs[vi][myfield])

            self.pointing[-1] = np.array(self.pointing[-1])

        for vi, v in enumerate(self.vis):
            if len(self.scans[vi]) > 0:
                goodscid = [x for x in self.scans[vi] if x in self.sourscans[vi]]
                if len(goodscid) != len(self.scans[vi]):
                    badscid = [x for x in self.scans[vi] if x not in self.sourscans[vi]]
                    msg = '\n ERROR! The following scans do NOT correspond to source %s:\n' % (str(self.field))
                    msg += str(badscid)
                    self._printError(msg)
                    return False
                self.sourscans[vi] = list(goodscid)

        # Get info on spectral configuration and polarization:
        # This code is more complicated than it should, since we
        # want to leave it ready to implement further features:
        self.spwlist = []
        self.pol2aver = []
        self.polmod = []
        self.polii = []

        spwi = 0
        for vi, v in enumerate(self.vis):
            j = {True: 0, False: vi}[len(self.spw) == 1]

            success = ms.open(v)
            if not success:
                self._printError("Failed to open measurement set '+v+'!")
                return False

            freqdic = ms.getspectralwindowinfo()
            spwchans = ms.range(['num_chan'])['num_chan']

            aux = channeler(self.spw[j], width=self.chanwidth, maxchans=spwchans)
            if aux[0]:
                ranges = list(aux[1])
            else:
                self._printError(aux[1] + "\nSomething seems to be wrong with the 'spw' number %i. \n ABORTING!" % (vi))
                return False

            nspws = range(len(spwchans))

            # These are the spws with selected channels for the fit:
            selspws = list([list([i, ranges[i]]) for i in nspws if len(ranges[i]) > 0])
            self.spwlist.append(list([j, vi, selspws, spwi]))
            # spwlist[vi][2] es una lista con [spwid, chanranges]
            # spwlist[vi][3]+index(spwlist[vi][2]) es la spw "local"
            spwi += {True: 0, False: len(selspws)}[len(self.spw) == 1]

            #################
            # Deal with polarization:
            ms.selectinit(datadescid=selspws[0][0])
            polprods = [x[0] for x in list(ms.range(['corr_names'])['corr_names'])]
            self.pol2aver.append(np.zeros(len(polprods)))
            # 0: normal,   1: multiply by i,   2: pol. independent, 3: just one product
            self.polmod.append(0)

            if self.stokes not in ['I', 'Q', 'U', 'V'] + polprods:
                self._printError("\n Bad Stokes parameter %s" % self.stokes)
                return False

            # User asks for one correlation product:
            if self.stokes in polprods:
                self.pol2aver[-1][polprods.index(self.stokes)] = 1.0
                self.polii.append([polprods.index(self.stokes)])
                self.polmod[-1] = 3

            # User asks for a Stokes parameter:
            # CASE 1: Circular feeds.
            elif 'RR' in polprods:
                try:
                    if self.stokes == 'I':
                        self.polii.append([polprods.index('RR'), polprods.index('LL')])
                        self.pol2aver[-1][polprods.index('RR')] = 0.5
                        self.pol2aver[-1][polprods.index('LL')] = 0.5
                    if self.stokes == 'PI':
                        self.polii.append([polprods.index('RR'), polprods.index('LL')])
                        self.pol2aver[-1][polprods.index('RR')] = 0.5
                        self.pol2aver[-1][polprods.index('LL')] = 0.5
                        self.polmod[-1] = 2
                    if self.stokes == 'Q':
                        self.polii.append([polprods.index('RL'), polprods.index('LR')])
                        self.pol2aver[-1][polprods.index('RL')] = 0.5
                        self.pol2aver[-1][polprods.index('LR')] = 0.5
                    if self.stokes == 'U':
                        self.polii.append([polprods.index('RL'), polprods.index('LR')])
                        self.pol2aver[-1][polprods.index('RL')] = 0.5
                        self.pol2aver[-1][polprods.index('LR')] = 0.5
                        self.polmod[-1] = 1
                    if self.stokes == 'V':
                        self.polii.append([polprods.index('RR'), polprods.index('LL')])
                        self.pol2aver[-1][polprods.index('RR')] = 0.5
                        self.pol2aver[-1][polprods.index('LL')] = -0.5
                except:
                    self._printError("Cannot convert to '+self.stokes+'!!")
                    return False
            #  CASE 2: Linear feeds.
            elif 'XX' in polprods:
                try:
                    if self.stokes == 'I':
                        self.polii.append([polprods.index('XX'), polprods.index('YY')])
                        self.pol2aver[-1][polprods.index('XX')] = 0.5
                        self.pol2aver[-1][polprods.index('YY')] = 0.5
                    if self.stokes == 'PI':
                        self.polii.append([polprods.index('XX'), polprods.index('YY')])
                        self.pol2aver[-1][polprods.index('XX')] = 0.5
                        self.pol2aver[-1][polprods.index('YY')] = 0.5
                        self.polmod[-1] = 2
                    if self.stokes == 'Q':
                        self.polii.append([polprods.index('XX'), polprods.index('YY')])
                        self.pol2aver[-1][polprods.index('XX')] = 0.5
                        self.pol2aver[-1][polprods.index('YY')] = -0.5
                    if self.stokes == 'U':
                        self.polii.append([polprods.index('XY'), polprods.index('YX')])
                        self.pol2aver[-1][polprods.index('XY')] = 0.5
                        self.pol2aver[-1][polprods.index('YX')] = 0.5
                    if self.stokes == 'V':
                        self.polii.append([polprods.index('YX'), polprods.index('XY')])
                        self.pol2aver[-1][polprods.index('YX')] = 0.5
                        self.pol2aver[-1][polprods.index('XY')] = -0.5
                        self.polmod[-1] = 1
                except:
                    self._printError("Cannot convert to '+self.stokes+'!!")
                    return False
            else:
                self._printError("Polarization " + self.stokes + " not understood.\n ABORTING!")
                return False

            ms.close()

        # Try to compile the equations for the variables:
        #  if data_changed:
        #  try:
        #    del self.mymodel
        #  except:
        #    self.printInfo("UVMULTIFIT model does not seem to exist yet.")
        #    pass

        #  self.mymodel = modeler(self.model, self.var, self.fixed, self.fixedvar, self.scalefix,
        #                         self.NCPU, self.only_flux, self.applyHankel, self.isNumerical,
        #                         self.useGains, [self.phase_gains, self.amp_gains])

        #  if self.mymodel.failed:
        #    self.printError(self.mymodel.resultstring)
        #    return False
        #
        #  if not data_changed:
        #    self.initData()
        self._setEngineWgt()
        return self.success

    ############################################
    #
    #  PREPARE THE UVMOD LIBRARY AND SET THE WEIGHTING
    #
    def _setEngineWgt(self):
        """Prepare the PB correction, the model and the number of cores for the modeler.
        Not to be called directly by the user."""

        # Load the C++ library:
        # if self.mymodel.Ccompmodel is None:
        #     self.mymodel.Ccompmodel = uvmod.modelcomp
        self.mymodel.NCPU = self.NCPU
        self.success = self._setWgtEq()

        return self.success

    ############################################
    #
    #  SET THE WEIGHTING AND PB-CORRECTION
    #
    def _setWgtEq(self):
        """Depends on setEngineWgt. Not to be called by the user."""

        tempfloat = 0.0

        if not isinstance(self.dish_diameter, dict):
            try:
                self.dish_diameter = float(self.dish_diameter)
            except ValueError:
                self._printError("\n The dish diameter must be a number! (in meters)\n")
                return False

        tb.open(os.path.join(self.vis[0], 'ANTENNA'))
        self.antnames = tb.getcol('NAME')

        if self.pbeam:
            self._printInfo("""\n
You selected to apply primary-beam correction.
PLEASE, remember that the beam is being approximated
with a Gaussian, so it may not be very accuracte far
from the pointing direction.\n\n""")
            if isinstance(self.dish_diameter, float):
                if self.dish_diameter == 0.0:
                    try:
                        tempfloat = np.copy(tb.getcol('DISH_DIAMETER'))
                    except:
                        self._printInfo("\n\n Dish Diameter column NOT found in antenna tables!\n")
                    tempfloat = np.zeros(len(self.antnames))
                else:
                    self._printInfo("\n\n An antenna diameter of %.3f m will be applied\n" % self.dish_diameter)
                    tempfloat = np.array([self.dish_diameter for a in self.antnames])

            elif isinstance(self.dish_diameter, dict):
                tempfloat = np.array([0.0 for a in self.antnames])
                for anam in self.dish_diameter.keys():
                    for anid in range(len(self.antnames)):
                        antids = re.search(anam, self.antnames[anid])
                        if 'start' in dir(antids):
                            tempfloat[anid] = self.dish_diameter[anam]
                self._printInfo("\n\nManual antenna-size setting:\n")
                for anid in range(len(self.antnames)):
                    self._printInfo("  Antenna %s has a diameter of %.2fm.\n" % (self.antnames[anid], tempfloat[anid]))

            else:
                self._printError("\n\n BAD dish_diameter! Should be a float or a dict!\n")

            if np.max(tempfloat) == 0.0:
                self._printError("\nThe antenna diameters are not set in the ms.\n" +
                                 "Please, set it manually or turn off primary-beam correction.\n")
                return False
            else:
                # Negative means not to apply PB corr for that antenna
                tempfloat[tempfloat == 0.0] = -1.0
                FWHM = self.ldfac / tempfloat * (2.99e8)
                sigma = FWHM / 2.35482 * (180. / np.pi) * 3600.
                self.mymodel.KfacWgt = 1. / (2. * sigma**2.) * (tempfloat > 0.0)  # (0.5*(tempfloat/1.17741)**2.)
                self.userDiameters = tempfloat
        else:
            self.mymodel.KfacWgt = np.zeros(len(self.antnames))

        tb.close()
        # May refine this function in future releases:
        #  self.mymodel.wgtEquation = lambda D, Kf: -D*Kf
        return True

    ############################################
    #
    #  READ THE DATA. ARRANGE ALL ARRAYS
    #
    def readData(self, del_data=True):
        """ Reads the data, according to the properties ``vis, column, chanwidth``, etc.

        It then fills in the properties ``averdata, averfreqs, averweights, u, v, w``, etc.

        Each one of these properties is a list with the data (one list item per spectral window/scan).

        A previous successful run of function ``checkInputs()`` is assumed.

        .. note:: Instead of re-reading data from scratch, using the same uvmultifit instance, a
        better approach may be to restart CASA and create a fresh uvmultifit instance with the
        new data, avoiding some memory leakage related to potential hidden references to the
        data in the IPython's *recall* prompt."""

        tic = time.time()

        if del_data:  # data_changed:
            self._deleteData()
            #    self.clearPointers(0)
            #    self.clearPointers(1)

        self.success = False

        # Initiate the lists and arrays where the data will be read-in:
        ntotspw = self.spwlist[-1][3] + len(self.spwlist[-1][2])
        nsprang = range(ntotspw)
        self.u = [[] for sp in nsprang]
        self.v = [[] for sp in nsprang]
        self.w = [[] for sp in nsprang]
        self.t = [[] for sp in nsprang]
        self.tArr = [[] for sp in nsprang]
        self.tIdx = [[] for sp in nsprang]
        self.ant1 = [[] for sp in nsprang]
        self.ant2 = [[] for sp in nsprang]
        self.RAshift = [[] for sp in nsprang]
        self.Decshift = [[] for sp in nsprang]
        self.Stretch = [[] for sp in nsprang]
        self.averdata = [[] for sp in nsprang]
        self.avermod = [[] for sp in nsprang]
        self.averweights = [[] for sp in nsprang]
        self.averfreqs = [[] for sp in nsprang]
        self.iscancoords = [[] for sp in nsprang]

        maxDist = 0.0

        # Read data for each spectral window:
        self.iscan = {}
        for vi in self.vis:
            self.iscan[vi] = {}

        for si in nsprang:
            # These are temporary lists of arrays that will be later concatenated:
            datascanAv = []
            modelscanAv = []
            #   datascanim = []
            weightscan = []
            uscan = []
            vscan = []
            wscan = []
            tscan = []
            tArray = []
            tIndex = []
            ant1scan = []
            ant2scan = []
            RAscan = []
            Decscan = []
            Stretchscan = []

            i0scan = 0

            # BEWARE! was sp
            for vidx, vis in enumerate(filter(lambda x: x[3] <= si, self.spwlist)):

                msname = self.vis[vis[1]]

                self._printInfo("\n\n Opening measurement set " + msname + ".\n")
                tb.open(msname)
                SPW = tb.getcol('DATA_DESC_ID')
                crosscorr = tb.getcol('ANTENNA1') != tb.getcol('ANTENNA2')
                tb.close()

                tb.open(os.path.join(msname, 'DATA_DESCRIPTION'))
                DDSC = tb.getcol('SPECTRAL_WINDOW_ID')
                tb.close()

                for spidx, spi in enumerate(vis[2]):
                    if vis[3] + spidx == si:

                        sp = spi[0]
                        rang = spi[1]
                        self.iscan[msname][sp] = {}

                        DDs = np.where(DDSC == sp)[0]

                        if len(DDs) > 1:
                            self._printInfo("WARNING! spw %i has more than one Data Description ID!" % sp)

                        maskspw = np.zeros(np.shape(crosscorr), dtype=np.bool)

                        for ddi in DDs:
                            maskspw = np.logical_or(maskspw, SPW == ddi)

                        maskspw *= crosscorr

                        # For the first ms in the list, read the frequencies of the spw.
                        # All the other mss will be assumed to have the same frequencies:
                        if True:
                            tb.open(os.path.join(msname, 'SPECTRAL_WINDOW'))
                            origfreqs = tb.getcol('CHAN_FREQ')
                            self.averfreqs[si] = np.array([np.average(origfreqs[r]) for r in rang])
                            nfreq = len(rang)
                            tb.close()
                        self._printInfo("\n\n Reading scans for spw %i \n" % sp)

                        # Read all scans for this field id:
                        for sc, scan in enumerate(self.sourscans[vis[1]]):
                            tb.open(msname)

                            masksc = maskspw * (tb.getcol('SCAN_NUMBER') == int(scan))
                            fieldids = list(np.sort(np.unique(tb.getcol('FIELD_ID')[masksc])))

                            tb.close()

                            for fieldid in fieldids:
                                self._printInfo("\r Reading scan #%i (%i of %i). Field: %i" %
                                                (scan, sc + 1, len(self.sourscans[vis[1]]), fieldid))
                                tb.open(msname)
                                maskfld = np.where(masksc * (tb.getcol('FIELD_ID') == int(fieldid)))[0]

                                if len(maskfld) == 0:
                                    tb.close()
                                else:
                                    tb2 = tb.selectrows(maskfld)
                                    uvscan = {'uvw': tb2.getcol('UVW'), 'antenna1': tb2.getcol('ANTENNA1'),
                                              'antenna2': tb2.getcol('ANTENNA2'), 'time': tb2.getcol('TIME')}
                                    times = np.unique(uvscan['time'])
                                    if self.takeModel:
                                        # datascan = ms.getdata([self.column, 'model_data', 'weight', 'flag'],
                                        #                       ifraxis=True)
                                        datascan = {self.column: tb2.getcol((self.column).upper()),
                                                    'model_data': tb2.getcol('MODEL_DATA'),
                                                    'weight': tb2.getcol('WEIGHT'),
                                                    'flag': tb2.getcol('FLAG')}
                                    else:
                                        # datascan = ms.getdata([self.column, 'weight', 'flag'], ifraxis=True)
                                        datascan = {self.column: tb2.getcol((self.column).upper()),
                                                    'weight': tb2.getcol('WEIGHT'), 'flag': tb2.getcol('FLAG')}

                                    tb.close()

                                    # NOTE: There is a bug in np.ma.array that casts complex
                                    # to float under certain operations (e.g., np.ma.average).
                                    # That's why we average real and imag separately.

                                    # Compute the polarization product:

                                    # Bad data has zero weight:
                                    datascan['weight'][np.logical_not(np.isfinite(datascan['weight']))] = 0.0

                                    # All unflagged weights set to equal (if uniform):
                                    if self.uniform:
                                        datascan['weight'][datascan['weight'] > 0.0] = 1.0

                                    datascan['weight'][datascan['weight'] < 0.0] = 0.0
                                    copyweight = np.copy(datascan['weight'])

                                    totalmask = datascan['flag']
                                    origmasked = np.ma.array(datascan[self.column], mask=totalmask, dtype=np.complex128)

                                    if self.takeModel:
                                        origmodmasked = np.ma.array(
                                            datascan['model_data'],
                                            mask=totalmask, dtype=np.complex128)

                                    # The weights are weighting the RESIDUALS, and not the ChiSq terms.
                                    # Hence, we divide wgt_power by 2.:
                                    origweight = np.power(copyweight, self.wgt_power / 2.)

                                    # Completely flagged times/baselines:
                                    origweight[np.sum(np.logical_not(totalmask), axis=1) == 0] = 0.0

                                    datamask = 0.0
                                    weightmask = 0.0
                                    flagmask = 0.0
                                    modelmask = 0.0

                                    # Construct the required polarization from the correlation products:
                                    polavg = [pol != 0.0 for pol in self.pol2aver[vis[1]]]

                                    if self.polmod[vis[1]] == 2:
                                        flagmask = np.ma.logical_and(
                                            totalmask[self.polii[vis[1]][0], :, :],
                                            totalmask[self.polii[vis[1]][1], :, :])
                                        datamask = np.ma.average(origmasked[polavg, :].real, axis=0) + \
                                            1.j * np.ma.average(origmasked[polavg, :].imag, axis=0)

                                        if self.takeModel:
                                            modelmask = np.ma.average(origmodmasked[polavg, :].real, axis=0) + \
                                                1.j * np.ma.average(origmasked[polavg, :].imag, axis=0)

                                        weightmask = np.sum(origweight[polavg, :], axis=0)
                                    else:
                                        if self.polmod[vis[1]] == 3:
                                            flagmask = totalmask[self.polii[vis[1]][0], :, :]
                                        else:
                                            flagmask = np.ma.logical_or(totalmask[self.polii[vis[1]][0], :, :],
                                                                        totalmask[self.polii[vis[1]][1], :, :])

                                        for pol in self.polii[vis[1]]:
                                            datamask += origmasked[pol, :] * self.pol2aver[vis[1]][pol]
                                            if self.takeModel:
                                                modelmask += origmodmasked[pol, :] * self.pol2aver[vis[1]][pol]
                                            weightmask += origweight[pol, :]

                                        if self.polmod[vis[1]] == 1:
                                            datamask *= 1.j
                                            if self.takeModel:
                                                modelmask *= 1.j
                                                #   weightmask[flagmask] = 0.0

                                    # Free some memory:
                                    for key in list(datascan):
                                        del datascan[key]
                                    del datascan, origmasked, origweight
                                    del copyweight, totalmask, flagmask
                                    if self.takeModel:
                                        del origmodmasked

                                    # Mosaic-related corrections:
                                    phshift = 3600. * 180. / np.pi * (self.phasedirs[vis[1]][fieldid] - self.refpos)
                                    strcos = np.cos(self.phasedirs[vis[1]][fieldid][1])

                                    if phshift[0] != 0.0 or phshift[1] != 0.0:
                                        self._printInfo(
                                            "\r\t\t\t\t Offset: %.2e RA (tsec) %.2e Dec (asec)." %
                                            (phshift[0] / 15., phshift[1]))

                                    # Average spectral channels:
                                    nchannel, ndata = np.shape(datamask)
                                    ntimes = len(times)
                                    ntav = int(max([1, round(float(ntimes) / self.timewidth)]))
                                    datatemp = np.ma.zeros((nfreq, ndata), dtype=np.complex128)
                                    if self.takeModel:
                                        modeltemp = np.ma.zeros((nfreq, ndata), dtype=np.complex128)
                                    weighttemp = np.ma.zeros((nfreq, ndata))

                                    if self.chanwidth == 1:
                                        concRan = [c[0] for c in rang]
                                        datatemp[:, :] = datamask[concRan, :]
                                        if self.takeModel:
                                            modeltemp[:, :] = modelmask[concRan, :]
                                        weighttemp[:, :] = weightmask[np.newaxis, :]
                                    else:
                                        for nu in range(nfreq):
                                            datatemp[nu, :] = np.ma.average(
                                                datamask[rang[nu], :].real, axis=0) + \
                                                1.j * np.ma.average(datamask[rang[nu], :].imag, axis=0)
                                            if self.takeModel:
                                                modeltemp[nu, :] = np.ma.average(
                                                    modelmask[rang[nu], :].real, axis=0) + \
                                                    1.j * np.ma.average(modelmask[rang[nu], :].imag, axis=0)
                                            weighttemp[nu, :] = weightmask

                                    # Average in time and apply uvtaper:
                                    GaussWidth = 2. * (self.uvtaper / 1.17741)**2.
                                    RAoffi = np.zeros(np.shape(uvscan['time']))
                                    Decoffi = np.copy(RAoffi)
                                    Stretchi = np.copy(RAoffi)

                                    #  if self.timewidth ==1:

                                    RAoffi[:] = float(phshift[0])
                                    Decoffi[:] = float(phshift[1])
                                    Stretchi[:] = float(strcos)

                #########################
                # CODE FOR TIMEWIDTH>1 HAS TO BE BASED ON TB TOOL. WORK IN PROGRESS
                #     else:

                #       ant1s = uvscan['antenna1'][crosscorr]
                #       ant2s = uvscan['antenna2'][crosscorr]
                #       maskan1 = np.zeros(np.shape(ant1s), dtype=np.bool)
                #       mask = np.copy(maskan1)
                #       mask2 = np.copy(mask)
                #       # Antennas participating in this scan:
                #       allants1 = np.unique(ant1s)
                #       allants2 = np.unique(ant2s)
                #       # Fill in time-averaged visibs:
                #       for nt in range(ntav):
                #         t0 = nt*self.timewidth ; t1 = min([ntimes,(nt+1)*self.timewidth])
                #         mask[:] = (uvscan['time']>=times[t0])*(uvscan['time']<times[t1])*crosscorr
                #         for an1 in allants1:
                #           maskan1[:] = (ant1s==an1)*mask
                #           for an2 in allants2:
                #             if an2>an1:
                #                 mask2[:] = maskan1*(ant2s==an2)
                #                 uu = uvscan['uvw'][0, mask2]
                #                 vv = uvscan['uvw'][1, mask2]
                #                 ww = uvscan['uvw'][2, mask2]
                #                 ui[:, nt] = np.average(uu, axis=1)  # Baseline dimension
                #                 vi[:, nt] = np.average(vv, axis=1)  # Baseline dimension
                #                 wi[:, nt] = np.average(ww, axis=1)  # Baseline dimension
                #                 ant1i[:, nt] = ant1s  # Baseline dimension
                #                 ant2i[:, nt] = ant2s  # Baseline dimension
                #                 timei[:, nt] = np.average(uvscan['time'][t0:t1])/86400.
                #                 tArrayi[nt] = time[0, nt]  # np.average(uvscan['time'][t0:t1])/86400.
                #                 tIndexi[:, nt] = nt
                #                 RAoffi[:, nt] = float(phshift[0])
                #                 Decoffi[:, nt] = float(phshift[1])
                #                 Stretchi[:, nt] = float(strcos)
                #
                #         if self.uvtaper > 0.0:
                #           GaussFact = np.exp(-(uu*uu + vv*vv)/GaussWidth)
                #         else:
                #           GaussFact = np.ones(np.shape(uu))
                #
                #       broadwgt = weighttemp[:, :, t0:t1]
                #       avercompl[:, :, nt] = np.ma.average(datatemp[:, :, t0:t1].real, axis=2, weights=broadwgt)+
                #                             1.j*np.ma.average(datatemp[:, :, t0:t1].imag, axis=2, weights=broadwgt)
                #       if self.takeModel:
                #         avermodl[:, :, nt] = np.ma.average(modeltemp[:, :, t0:t1].real, axis=2, weights=broadwgt)+
                #                              1.j*np.ma.average(modeltemp[:, :, t0:t1].imag, axis=2, weights=broadwgt)
                #
                #       if self.uniform:
                #         averwgt[:, :, nt] = np.ma.sum(np.ones(np.shape(broadwgt))*GaussFact[np.newaxis, :, :], axis=2)
                #       else:
                #         averwgt[:, :, nt] = np.ma.sum(broadwgt*GaussFact[np.newaxis, :, :], axis=2)
                #########################

                                    ant1scan.append(np.copy(uvscan['antenna1'][:]))
                                    ant2scan.append(np.copy(uvscan['antenna2'][:]))
                                    uscan.append(np.copy(uvscan['uvw'][0, :]))
                                    vscan.append(np.copy(uvscan['uvw'][1, :]))
                                    wscan.append(np.copy(uvscan['uvw'][2, :]))
                                    tscan.append(np.copy(uvscan['time'][:]))
                                    tArray.append(np.copy(times))

                                    tIndexi = np.zeros(np.shape(uvscan['time']), dtype=np.int32)
                                    for tid, tiii in enumerate(times):
                                        tIndexi[uvscan['time'] == tiii] = tid

                                    if len(tIndex) > 1:
                                        tIndexi += np.max(tIndex[-1]) + 1
                                    tIndex.append(tIndexi)

                                    RAscan.append(RAoffi)
                                    Decscan.append(Decoffi)
                                    Stretchscan.append(Stretchi)

                                    datascanAv.append(np.transpose(datatemp))
                                    if self.takeModel:
                                        modelscanAv.append(np.transpose(modeltemp))
                                    if self.uniform:
                                        weightscan.append(np.transpose(np.ones(np.shape(weighttemp))))
                                    else:
                                        weightscan.append(np.transpose(weighttemp))

                                    # Useful info for function writeModel() and for pointing correction:
                                    if scan not in self.iscan[msname][sp].keys():
                                        self.iscan[msname][sp][scan] = []

                                    self.iscan[msname][sp][scan].append(
                                        [int(si), int(i0scan), int(ndata), list(rang), np.copy(maskfld)])
                                    self.iscancoords[si].append([i0scan, i0scan + ndata, phshift[0], phshift[1]])

                                    i0scan += ndata

            # Concatenate all the scans in one single array. Notice that we separate real and imag and save them
            # as floats. This is because ctypes doesn' t handle complex128.
            self.averdata[si] = np.require(np.concatenate(datascanAv, axis=0),
                                           requirements=['C', 'A'])  # , np.concatenate(datascanim, axis=0)]

            if self.takeModel:
                self.avermod[si] = np.require(np.concatenate(modelscanAv, axis=0), requirements=['C', 'A'])

            self.averweights[si] = np.require(np.concatenate(weightscan, axis=0), requirements=['C', 'A'])
            self.u[si] = np.require(np.concatenate(uscan, axis=0), requirements=['C', 'A'])
            self.v[si] = np.require(np.concatenate(vscan, axis=0), requirements=['C', 'A'])
            self.w[si] = np.require(np.concatenate(wscan, axis=0), requirements=['C', 'A'])
            self.t[si] = np.require(np.concatenate(tscan, axis=0), requirements=['C', 'A'])
            self.tArr[si] = np.require(np.concatenate(tArray, axis=0), requirements=['C', 'A'])
            self.tIdx[si] = np.require(np.concatenate(tIndex, axis=0), requirements=['C', 'A'])

            self.RAshift[si] = np.require(np.concatenate(RAscan, axis=0), requirements=['C', 'A'])
            self.Decshift[si] = np.require(np.concatenate(Decscan, axis=0), requirements=['C', 'A'])
            self.Stretch[si] = np.require(np.concatenate(Stretchscan, axis=0), requirements=['C', 'A'])
            self.ant1[si] = np.require(np.concatenate(ant1scan, axis=0), requirements=['C', 'A'])
            self.ant2[si] = np.require(np.concatenate(ant2scan, axis=0), requirements=['C', 'A'])

            # Free some memory:
            #   del uu, vv, ww, avercomplflat, weightflat
            del datatemp, weighttemp, uvscan  # , avercompl, averwgt
            for dda in datascanAv:
                del dda
            #   for dda in datascanim:
            #     del dda
            if self.takeModel:
                for dda in modelscanAv:
                    del dda
            for dda in weightscan:
                del dda
            for dda in tscan:
                del dda
            for dda in uscan:
                del dda
            for dda in vscan:
                del dda
            for dda in wscan:
                del dda
            for dda in RAscan:
                del dda
            for dda in Decscan:
                del dda
            for dda in Stretchscan:
                del dda
            for dda in ant1scan:
                del dda
            for dda in ant2scan:
                del dda
            for dda in tArray:
                del dda
            for dda in tIndex:
                del dda

            del datascanAv  # , datascanim
            del weightscan, tscan, uscan, vscan, wscan, tArray, tIndex
            del RAscan, Decscan, Stretchscan, ant1scan, ant2scan
            if self.takeModel:
                del modelscanAv

            try:
                del GaussFact
            except:
                pass
            gc.collect()

        # Initial and final times of observations (time reference for proper motions):
        self.t0 = np.min([np.min(ti) for ti in self.t])
        self.t1 = np.max([np.max(ti) for ti in self.t])

        self._printInfo("\n\n Done reading\n")
        tac = time.time()
        self._printInfo("\nReading took %.2f seconds.\n" % (tac - tic))

        self.success = True

        NIFs = len(self.averdata)
        self.Nspw = NIFs

        for i in range(self.Nspw):
            self._printInfo("There are %i integrations (%i visibs.) in spw %i\n" %
                            (len(self.tArr[i]), len(self.t[i]), i))

        #  import pickle as pk
        #  kk = open('DEBUG.dat', 'w')
        #  pk.dump([self.tIdx, self.tArr, self.t, self.ant1, self.ant2, self.u, self.v, self.w, self.averdata], kk)
        #  kk.close()
        #  raw_input('INPUT')

        #  self.clearPointers(0)

        # Set pointers to data, model, etc.:
        self.initData(del_data=del_data)

        return True

    ############################################
    #
    #  COMPUTE MODELS TO BE DIRECTLY SUBTRACTED FROM THE DATA
    #
    def computeFixedModel(self):
        """ Computes the value of the fixed model on all the u-v data points.

        It saves the results in the 'mymodel.output' property, which should have been zeroed previously.
        Notice that the fixed model is recomputed for each spectral channel (i.e., if
        OneFitPerChannel=True), and is computed only once if OneFitPerChannel==False."""

        if self.takeModel:
            self._printInfo("\nFixed model was taken from model column.\n" +
                            "Notice that if you are RE-FITTING, you'll need to RELOAD the model column!\n\n")
        else:
            self._printInfo("\nGoing to compute fixed model (may need quite a bit of time)\n")
            self.mymodel.residuals([0], mode=0)

    ############################################
    #
    #  COMPUTE CHI SQUARE FOR A PARTICULAR REALIZATION OF THE FIT
    #
    def chiSquare(self, p):
        """ Returns a list with 2 items: The Chi Square value, computed for the parameter values
        given in the ``p`` list) and the number of degrees of freedom.

        :Parameters:
        ------------
        **p** : `list`
        List of parameter values where the Chi square will be computed."""

        return self.mymodel.residuals(p, mode=-2)

    ############################################
    #
    #  SET SOME DATA (AND MODEL) ARRAYS
    #
    def initData(self, isUVfit=True, del_data=True):
        """ Initiates the data pointers of the 'modeler' instance.

        The 'modeler' instance stores pointers to the data and metadata, the compiled model (and
        fixedmodel), the parameter values, and all the methods to compute residuals, Chi Squared, etc."""

        # Reset pointers for the modeler:
        self.mymodel._deleteData()

        #  if del_data:  # data_changed:
        #    self.clearPointers(0)

        # Set number of spectral windows:
        gooduvm = uvmod.setNspw(int(self.Nspw))

        if gooduvm != 0:
            self._printError("\nError in the C++ extension!\n")
            return False

        # Maximum number of frequency channels (i.e., the maximum from all the selected spws):
        self.maxnfreq = 0

        # Fill in data pointers for each spw:
        for spidx in range(self.Nspw):
            self.mymodel.data.append([])
            self.mymodel.dt.append([])
            self.mymodel.dtArr.append([])
            self.mymodel.dtIdx.append([])
            self.mymodel.wgt.append([])
            self.mymodel.wgtcorr.append([])
            self.mymodel.uv.append([])
            self.mymodel.offset.append([])
            self.mymodel.output.append([])
            self.mymodel.fixedmodel.append([])
            self.mymodel.ants.append([])

            self.mymodel.fittable.append([])
            self.mymodel.fittablebool.append([])
            self.mymodel.isGain.append([])

            self.mymodel.iscancoords = self.iscancoords

            self.mymodel.freqs.append(np.require(self.averfreqs[spidx], requirements=['C', 'A']))
            self.maxnfreq = max(self.maxnfreq, len(self.averfreqs[spidx]))

            # Only have to multiply by freq, to convert these into lambda units:
            ulambda = np.require(self.FouFac * self.u[spidx] / self.LtSpeed, requirements=['C', 'A'])
            vlambda = np.require(self.FouFac * self.v[spidx] / self.LtSpeed, requirements=['C', 'A'])
            wlambda = np.require(self.FouFac * self.w[spidx] / self.LtSpeed, requirements=['C', 'A'])

            # Data, uv coordinates, and weights:

            # Data are saved in two float arrays (i.e., for real and imag):
            self.mymodel.data[-1] = self.averdata[spidx]  # [0], self.averdata[spidx][1]]
            self.mymodel.wgt[-1] = self.averweights[spidx]
            self.mymodel.uv[-1] = list([ulambda, vlambda, wlambda])
            self.mymodel.offset[-1] = list([self.RAshift[spidx], self.Decshift[spidx], self.Stretch[spidx]])
            self.mymodel.ants[-1] = list([self.ant1[spidx], self.ant2[spidx]])

            PBFactor = -np.sqrt(self.mymodel.KfacWgt[self.ant1[spidx]] * self.mymodel.KfacWgt[self.ant2[spidx]])

            # Set number of antennas and whether each one has a fittable gain:
            self.mymodel.Nants = self.Nants
            self.mymodel.isGain[-1] = np.require(np.zeros(len(self.t[spidx]), dtype=np.int8), requirements=['C', 'A'])
            for i in self.useGains:
                mask0 = self.mymodel.ants[-1][0] == i
                mask1 = self.mymodel.ants[-1][1] == i
                self.mymodel.isGain[-1][np.logical_or(mask0, mask1)] = True

            # Release memory:
            del ulambda, vlambda, wlambda

            # Time spent on observation:
            self.mymodel.dt[-1] = np.require(self.t[spidx] - self.t0, requirements=['C', 'A'])
            self.mymodel.t0 = self.t0
            self.mymodel.t1 = self.t1
            self.mymodel.dtArr[-1] = np.require(self.tArr[spidx] - self.t0, requirements=['C', 'A'])
            self.mymodel.dtIdx[-1] = np.require(self.tIdx[spidx], requirements=['C', 'A'])

            # Array to save the residuals (or model, or any output from the C++ library):
            self.mymodel.output[-1] = np.require(np.zeros(np.shape(self.averdata[spidx]),
                                                          dtype=np.complex128), requirements=['C', 'A'])
            if self.takeModel:
                try:
                    self.mymodel.output[-1][:] = self.avermod[spidx]
                except:
                    self._printError("You already used the model column! Should read it again!\n")

            ########
            # Array of booleans, to determine if a datum enters the fit:
            self.mymodel.fittable[-1] = self.mymodel.wgt[-1][:, 0] > 0.0
            # self.mymodel.wgtcorr[-1] = np.require(np.copy(np.tile(PBFactor,(len(self.model), 1))),
            #                                       requirements=['C', 'A'])
            self.mymodel.wgtcorr[-1] = np.require(np.copy(PBFactor), requirements=['C', 'A'])

            del PBFactor

            self.mymodel.fittablebool[-1] = np.require(np.copy(self.mymodel.fittable[-1]
                                                               ).astype(np.bool), requirements=['C', 'A'])

            gooduvm = uvmod.setData(spidx, self.mymodel.uv[-1][0], self.mymodel.uv[-1][1], self.mymodel.uv[-1][2],
                                    self.mymodel.wgt[-1], self.mymodel.data[-1],
                                    self.mymodel.output[-1], self.mymodel.freqs[-1],
                                    self.mymodel.fittable[-1], self.mymodel.wgtcorr[-1], self.mymodel.dt[-1],
                                    self.mymodel.dtArr[-1], self.mymodel.dtIdx[-1],
                                    self.mymodel.offset[-1][0], self.mymodel.offset[-1][1], self.mymodel.offset[-1][2],
                                    self.mymodel.ants[-1][0], self.mymodel.ants[-1][1],
                                    self.mymodel.isGain[-1], self.Nants)

            if gooduvm != 10:
                self._printError("\nError in the C++ extension!\n")
                return False

        try:
            for spidx in range(self.Nspw - 1, -1, -1):
                del self.avermod[spidx]
        except:
            pass
        gc.collect()

    ############################################
    #
    #  SET MODEL ARRAYS AND FILL MODEL DATA TO BE SENT TO MODELER
    #
    def initModel(self):
        """ Allocates memory for the modeler data, which will be used by the C++ extension.

        Also compiles the model variables. It is a good idea to run this method everytime that
        the user reads new data to be fitted (i.e., after a call to readData())

        It is indeed MANDATORY to run this function if the model to be fitted has changed."""

        if self.mymodel.initiated:
            #     self.clearPointers(1)
            self.mymodel._deleteModel()
        else:
            self._printInfo("UVMULTIFIT compiled model(s) do not seem to exist yet.\n")

        # Set number of threads:
        gooduvm = uvmod.setNCPU(int(self.NCPU))
        if gooduvm != 0:
            self._printError("\nError in the C++ extension!\n")
            return False

        # Initial setup modeler:
        self.mymodel._setup(self.model, self.var, self.fixed, self.fixedvar, self.scalefix,
                            self.NCPU, self.only_flux, self.applyHankel, self.isNumerical, self.useGains,
                            [self.phase_gainsNuT, self.amp_gainsNuT, self.phase_gainsNu,
                             self.amp_gainsNu, self.phase_gainsT, self.amp_gainsT],
                            self.isMixed)

        if self.mymodel.failed:
            self._printError(self.mymodel.resultstring)
            return False

        # Fill-in proper motions (in as/day)
        for i in range(len(self.model)):
            self.mymodel.propRA[i] = self.proper_motion[i][0] / 365.
            self.mymodel.propDec[i] = self.proper_motion[i][1] / 365.

        #####
        # Array to save the variables of the model and the 'scalefix' value,
        # all of them as a function of frequency. The C++ library will read the variables from here
        # at each iteration and for each model component:
        self.mymodel.initiated = True
        self.mymodel.varbuffer = [np.zeros((len(self.model),
                                            self.maxNvar + self.applyHankel, self.maxnfreq))
                                  for i in range(len(self.p_ini) + 1)]
        self.mymodel.varfixed = [np.zeros(self.maxnfreq) for i in range(len(self.p_ini) + 1)]
        self.mymodel.dpar = np.zeros(len(self.p_ini), dtype=np.float64)
        self.mymodel.par2 = np.zeros((3, len(self.p_ini)), dtype=np.float64)

        self.mymodel.Hessian = np.zeros((len(self.p_ini), len(self.p_ini)), dtype=np.float64)
        self.mymodel.Gradient = np.zeros(len(self.p_ini), dtype=np.float64)
        self.mymodel.imod = np.zeros(len(self.model), dtype=np.int32)

        # Compile equations that set the model variables and gains:
        self._printInfo("Going to compile models\n")
        self.mymodel._compileAllModels()

        self.mymodel.GainBuffer = [[[] for AI in range(self.Nants)] for spidx in range(self.Nspw)]

        for spidx in range(self.Nspw):
            for AI in range(self.Nants):
                if self.isMixed:
                    self.mymodel.GainBuffer[spidx][AI] = [np.ones(
                        (len(self.tArr[spidx]),
                         len(self.averfreqs[spidx])),
                        dtype=np.complex128) for i in range(len(self.mymodel.parDependence[AI]))]
                else:
                    self.mymodel.GainBuffer[spidx][AI] = [np.ones(
                        len(self.tArr[spidx]) + len(self.averfreqs[spidx]),
                        dtype=np.complex128) for i in range(len(self.mymodel.parDependence[AI]))]

        compFixed = len(self.fixed) > 0

        # Fill in all the information and data pointers for the modeler:
        self._printInfo("Going to run setModel\n")
        gooduvm = uvmod.setModel(self.mymodel.imod, self.mymodel.Hessian, self.mymodel.Gradient,
                                 self.mymodel.varbuffer, self.mymodel.varfixed, self.mymodel.dpar,
                                 self.mymodel.propRA, self.mymodel.propDec, self.refpos,
                                 self.mymodel.parDependence, self.mymodel.GainBuffer, compFixed, self.isMixed)

        # Delta to estimate numerical derivatives:
        self.mymodel.minnum = self.minDp

        if gooduvm != 10:
            self._printError("\nError in the C++ extension!\n")
            return False
        else:
            self._printInfo("Success!\n")

        goodinit = uvmod.setWork()
        if goodinit != 10:
            self._printError("Memory allocation error!")

        self.mymodel.bounds = self.bounds
        self.mymodel.LMtune = [float(lmi) for lmi in self.LMtune]

        return True

    ############################################
    #
    #  PERFORM (AND SAVE) THE FIT
    #
    def fit(self, redo_fixed=True, reinit_model=False, save_file=True, nuidx=-1, reset_flags=False):
        """ Fits the data, using the models previously compiled with ``initModel()``.

        Parameters
        ----------

        **reinit_model** : `bool`
          It is False by default. If False, the models used are those already compiled. If True,
          the models are *recompiled*, according to the contents of the ``model, var, fixed,
          fixedvar`` properties, and all the references to the data arrays are refreshed.

        **redo_fixed** : `bool`
          It is True by default. If False, the fixed model will **not** be recomputed throughout
          the fit.

          .. warning:: Setting ``redo_fixed=False`` can be dangerous if you are fitting in
             spectral-line mode and have channels with very different frequencies (i.e., a wide
             fractional bandwidth), since the UV coordinates will **not** be re-projected in that
             case.

        **save_file** : `bool`
          It is True by default. If False, the external ascii file with the fitting results will
          not be created.

        **nuidx** : `int`
          A helper parameter for internal use (if -1, fit to the continuum. Otherwise, fit only
          the ``nui`` channel). The user should NOT need to redefine (or force) this parameter when
          calling this function.

        **reset_flags** : `bool`
          Default is False. This is used to clear out the status of *bad data* that may have been
          set by spetial routines of **UVMultiFit** (e.g., the Quinn Fringe Fitter). Default is
          to NOT reset flags (this option should be OK most of the time)."""

        tic = time.time()

        self.mymodel.bounds = self.bounds
        self.mymodel.LMtune = [float(lmi) for lmi in self.LMtune]

        npars = len(self.p_ini)
        nspwtot = self.spwlist[-1][3] + len(self.spwlist[-1][2])

        notfit = [[] for si in range(nspwtot)]

        # Select data according to time range:
        datatot = 0
        self.allflagged = False

        for si in range(nspwtot):
            if reset_flags:
                self.mymodel.fittablebool[si][:] = True
            if self.MJDrange[0] > 0.0 and self.MJDrange[1] > 0.0:
                self._printInfo("\nSelecting data by Modified Julian Date\n")
                timeWindow = np.logical_and(np.logical_and(
                    self.t[si] >= self.MJDrange[0],
                    self.t[si] <= self.MJDrange[1]))
                self.mymodel.fittablebool[si][:] = np.logical_and(timeWindow, self.mymodel.fittablebool[si])
                del timeWindow

            self.mymodel.fittable[si][:] = self.mymodel.fittablebool[si][:].astype(np.int8)

        #    else:
        #
        #      self.mymodel.fittable[si][:] = 1
        #      self.mymodel.fittablebool[si][:] = 1

        # Check if there is data available:
            unflagged = np.sum(self.mymodel.wgt[si][self.mymodel.fittablebool[si], :] != 0.0, axis=0)
            ntot = np.sum(self.mymodel.fittablebool[si])
            if self.OneFitPerChannel:
                if np.sum(unflagged == 0.0) > 0:
                    self._printInfo(
                        "ERROR: NOT ENOUGH DATA FOR THIS TIME RANGE! \n CHANNELS: " +
                        str(list(np.where(unflagged == 0.0))))
                    self.allflagged = True
                    notfit[si] = list(np.where(unflagged == 0.0)[0])
                    #  return False
            else:
                datatot += np.sum(unflagged)

        #   self.mymodel.output[si][:] = 0.0
        if datatot == 0 and not self.OneFitPerChannel:
            self._printInfo("ERROR: NOT ENOUGH DATA FOR THIS TIME RANGE!\n")
            self.allflagged = True
            return False

        #  for si in range(nspwtot):
        #    self.mymodel.wgtcorr[si][:] = -self.mymodel.KfacWgt
        self._printInfo("\nNow, fitting model \n")

        # Initialize model:
        if reinit_model:
            if self.mymodel.initiated:
                goodinit = self.initModel()
            if self.mymodel.failed or goodinit is False:
                self._printError("ERROR: Bad model (re)initialization! \n%s" % self.mymodel.resultstring)
                return False

        for i in range(len(self.p_ini) + 1):
            self.mymodel.varbuffer[i][:] = 0.0
            self.mymodel.varfixed[i][:] = 0.0

        # FIT!!!
        ndata = 0.0

        ##################
        # CASE OF SPECTRAL-MODE FIT:
        if self.OneFitPerChannel:
            fitparams = [[] for j in range(nspwtot)]
            fiterrors = [[] for j in range(nspwtot)]
            covariance = [[] for j in range(nspwtot)]
            ChiSq = [[] for j in range(nspwtot)]
            Nvis = [[] for j in range(nspwtot)]

            # Fit channel-wise for each spw:
            for si in range(nspwtot):
                rang = np.shape(self.mymodel.wgt[si])[1]
                print("")

                for nuidx in range(rang):
                    self._printInfo("\r Fitting channel " + str(nuidx+1) + " of " + str(rang) + " in spw " + str(si))

                    self.mymodel.currspw = si
                    self.mymodel.currchan = nuidx

                    # Compute fixed model:
                    if redo_fixed and len(self.mymodel.fixed) > 0 and not self.takeModel:
                        self.computeFixedModel()

                    # Fit with simplex (if asked for):

                    if nuidx not in notfit[si]:

                        if self.method == 'simplex':
                            fitsimp = _mod_simplex(self.mymodel.ChiSquare, self.p_ini,
                                                   args=(self.bounds, self.p_ini), disp=False,
                                                   relxtol=self.SMPtune[0], relstep=self.SMPtune[1],
                                                   maxiter=self.SMPtune[2] * len(self.p_ini))
                            fit = [fitsimp[0], np.zeros((len(self.p_ini), len(self.p_ini))), fitsimp[1]]
                        else:
                            fit = self.mymodel.LMMin(self.p_ini)
                            if not fit:
                                return False
                        # Estimate the parameter uncertainties and save the model in the output array:
                        if self.savemodel:
                            if self.write_model == 1:
                                Chi2t = self.mymodel.residuals(fit[0], mode=-3)
                            elif self.write_model == 2:
                                Chi2t = self.mymodel.residuals(fit[0], mode=-4)
                            elif self.write_model == 3:
                                Chi2t = self.mymodel.residuals(fit[0], mode=-5)

                    # DON'T FIT THIS CHANNEL (ALL DATA FLAGGED)
                    else:
                        fit = [[0.0 for pi in self.p_ini], np.zeros((len(self.p_ini), len(self.p_ini))), 0]

                    fitparams[si].append([float(f) for f in fit[0]])

                    # Only add the unflagged data to compute the DoF
                    ndata = float(np.sum(self.mymodel.wgt[si][:, nuidx] > 0.0))

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
            self._printInfo("\nFitting to all frequencies at once.\n")

            # This will tell the modeller to solve in continuum mode:
            self.mymodel.currspw = -1
            self.mymodel.currchan = -1

            # Compute fixed model:
            if redo_fixed and len(self.mymodel.fixed) > 0 and not self.takeModel:
                print("Generating fixed model. May take some time")
                self.computeFixedModel()
                print("Done!")

            # Pre-fit with simplex:
            if self.method == 'simplex':
                fitsimp = _mod_simplex(self.mymodel.ChiSquare, self.p_ini,
                                       args=(self.bounds, self.p_ini), disp=False, relxtol=self.SMPtune[0],
                                       relstep=self.SMPtune[1], maxiter=self.SMPtune[2] * len(self.p_ini))
                fit = [fitsimp[0], np.zeros((len(self.p_ini), len(self.p_ini))), fitsimp[1]]
            else:

                # Bound least-squares fitting:
                fit = self.mymodel.LMMin(self.p_ini)
                if not fit:
                    return False

            # Estimate the parameter uncertainties and save the model in the output array:
            if self.savemodel:
                if self.write_model == 1:
                    Chi2t = self.mymodel.residuals(fit[0], mode=-3)
                elif self.write_model == 2:
                    Chi2t = self.mymodel.residuals(fit[0], mode=-4)
                elif self.write_model == 3:
                    Chi2t = self.mymodel.residuals(fit[0], mode=-5)

            fitparams = fit[0]

            for si in range(nspwtot):
                # Only add the unflagged data to compute the DoF
                ndata += float(np.sum(self.mymodel.wgt[si] > 0.0))

            if fit[2] > 0.0:
                ChiSq = fit[2] / ndata  # Reduced chi squared.
            else:
                # There are 0 'really-free' parameters?! Watch out!!
                ChiSq = fit[2]
            Nvis = ndata

            fiterrors = [np.sqrt(fit[1][i, i] * ChiSq) for i in range(npars)]
            covariance = fit[1] * ChiSq

        self._printInfo("\n The reduced Chi Squared will be set to 1 by re-scaling the visibility weights.\n")
        # Free some memory:
        gc.collect()

        #####
        # Set the 'result' property:
        if not self.OneFitPerChannel:
            to_return = {'Frequency': self.averfreqs[0][0], 'Parameters': np.array(fitparams),
                         'Uncertainties': np.array(fiterrors), 'Reduced Chi squared': ChiSq,
                         'Fit': fit, 'Degrees of Freedom': Nvis}
            prtpars = []
            for pp, ppar in enumerate(fitparams):
                prtpars.append(ppar)
                prtpars.append(fiterrors[pp])

        else:
            Freq = []  # {}
            Par = []  # {}
            Err = []  # {}
            Chi2 = []  # {}
            prtpars = [[[] for spi in range(len(fitparams[sp]))] for sp in range(nspwtot)]

            for sp in range(nspwtot):
                Freq.append(self.averfreqs[sp])
                Par.append(np.array(fitparams[sp]))
                Err.append(np.array(fiterrors[sp]))
                Chi2.append(np.array(ChiSq[sp]))
                for pp, ppar in enumerate(fitparams[sp]):
                    for ppi, ppari in enumerate(ppar):
                        prtpars[sp][pp].append(ppari)
                        prtpars[sp][pp].append(fiterrors[sp][pp][ppi])

            to_return = {'Frequency': Freq, 'Parameters': Par, 'Uncertainties': Err,
                         'Reduced Chi squared': Chi2, 'Fit': fit, 'Degrees of Freedom': Nvis}

        if self.cov_return:
            to_return['covariance'] = covariance

        self.result = to_return

        # Save results in external file:
        if not save_file:
            return True

        self._printInfo("\n DONE FITTING. Now, saving to file.\n")

        outf = open(self.outfile, 'w')
        outf.write("# MODELFITTING RESULTS FOR MS: " + ','.join(self.vis))
        outf.write("\n\n# NOTICE THAT THE REDUCED CHI SQUARE SHOWN HERE IS THE VALUE\n")
        outf.write("# *BEFORE* THE RE-SCALING OF THE VISIBILITY WEIGHTS.\n")
        outf.write("# HOWEVER, THE PARAMETER UNCERTAINTIES ARE *ALREADY* GIVEN\n")
        outf.write("# FOR A REDUCED CHI SQUARE OF 1.\n")
        if isinstance(Nvis, list):
            DOG = int(np.average(Nvis))
        else:
            DOG = Nvis
        outf.write("# AVG. NUMBER OF DEGREES OF FREEDOM: %i" % DOG)
        if self.pbeam:
            if isinstance(self.dish_diameter, float):
                outf.write("# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. USING A DISH DIAMETER OF: %.3f METERS" %
                           self.dish_diameter)
            else:
                outf.write("# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. USING: \n")
                for iant in range(len(self.antnames)):
                    outf.write("#    FOR ANTENNA %s A DIAMETER OF %.2fm \n" %
                               (self.antnames[iant], self.userDiameters[iant]))

        else:
            outf.write("# PRIMARY-BEAM CORRECTION HAS NOT BEEN APPLIED.")

        outf.write("\n\n###########################################\n")
        outf.write("#\n# MODEL CONSISTS OF:\n#")
        for m, mod in enumerate(self.model):
            outf.write("\n# '" + mod + "' with variables: " + self.var[m])
        if len(self.fixed) > 0:
            outf.write("\n#\n#\n# FIXED MODEL CONSISTS OF:\n#")
            if self.takeModel:
                outf.write("\n# THE MODEL COLUMN FOUND IN THE DATA")
            else:
                for m, mod in enumerate(self.fixed):
                    var2print = list(map(float, self.fixedvar[m].split(',')))
                    outf.write(
                        "\n# '" + mod + "' with variables: " + ' '.join(['%.3e'] * len(var2print)) % tuple(var2print))
                outf.write("\n#\n#  - AND SCALING FACTOR: %s" % self.scalefix)

        outf.write("\n#\n#\n# INITIAL PARAMETER VALUES:\n#")
        for p0i, p0 in enumerate(self.p_ini):
            if self.bounds is not None:
                outf.write("\n#  p[%i] = %.5e with bounds: %s " % (p0i, p0, str(self.bounds[p0i])))
            else:
                outf.write("\n#  p[%i] = %.5e with no bounds" % (p0i, p0))

        outf.write("\n#\n##########################################\n\n")
        parshead = []

        for pp in range(len(self.p_ini)):
            parshead.append(pp)
            parshead.append(pp)

        headstr = (
            "# Frequency (Hz)   " + "p[%i]  error(p[%i])   " * len(self.p_ini) + "Red. Chi Sq.\n") % tuple(parshead)
        outf.write(headstr)

        if not self.OneFitPerChannel:
            formatting = "%.12e   " + "%.4e " * (2 * npars) + "   %.4e \n"
            toprint = tuple([np.average(self.averfreqs)] + prtpars + [ChiSq])
            outf.write(formatting % toprint)

        else:
            formatting = "%.12e    " + "%.4e " * (2 * npars) + "  %.4e \n"

            for spwvi in self.spwlist:
                for r, rr in enumerate(spwvi[2]):
                    k = spwvi[3] + r
                    rang = rr[1]
                    for nu, freq in enumerate(self.averfreqs[k]):
                        toprint = tuple([freq] + prtpars[k][nu] + [ChiSq[k][nu]])
                        outf.write(formatting % toprint)
        outf.close()

        tac = time.time()

        self._printInfo("\n Fit took %.2f seconds.\n\n END!\n" % (tac - tic))
        #  uvmod.unsetWork()

        return True


############################################
#
# MODELER CLASS (CALLED BY UVMULTIFIT)
#
class modeler():
    """ Class to deal with model equations and fitting procedures.

    If convert strings representing models and parameters into compiled equations, to be used in a ChiSq
    visibility fitting. It also interacts with the C++ extension of UVMultiFit.

    This class should NOT be instantiated by the user (it is called from the UVMultiFit class."""

    ############################################
    #
    #  FREE MEMORY
    #
    def __del__(self):
        self._deleteData()
        self._deleteModel()

    ############################################
    #
    #  FREE MEMORY JUST FOR THE MODEL-RELATED DATA:
    #
    def _deleteModel(self):
        """ Free pointers to the model-related arrays and parameters."""

        for mdi in range(len(self.varbuffer) - 1, -1, -1):  # [::-1]:
            del self.varbuffer[mdi]
        #    del self.varbuffer

        for mdi in range(len(self.varfixed) - 1, -1, -1):  # [::-1]:
            del self.varfixed[mdi]
        #    del self.varfixed

        #    for mdi in range(len(self.dpar)-1, -1, -1):  # -[::-1]:
        #      del self.dpar[mdi]
        del self.dpar

        #    for mdi in range(len(self.par2)-1, -1, -1):  # [::-1]:
        #      del self.par2[mdi]
        del self.par2

        del self.Hessian, self.Gradient, self.imod

    def _deleteData(self):
        """ Free pointers to the data-related arrays and gain buffers."""

        for mdi in range(len(self.data) - 1, -1, -1):  # [::-1]:
            del self.data[mdi]
        #    del self.data

        for mdi in range(len(self.wgt) - 1, -1, -1):  # [::-1]:
            del self.wgt[mdi]
        #    del self.wgt

        for mdi in range(len(self.uv) - 1, -1, -1):  # [::-1]:
            try:
                del self.uv[mdi][2], self.uv[mdi][1], self.uv[mdi][0]
                del self.uv[mdi]
            except:
                pass
        #    del self.uv

        for mdi in range(len(self.offset) - 1, -1, -1):  # [::-1]:
            try:
                del self.offset[mdi][2], self.offset[mdi][1], self.offset[mdi][0]
                del self.offset[mdi]
            except:
                pass

        #    del self.offset
        for mdi in range(len(self.ants) - 1, -1, -1):  # [::-1]:
            try:
                del self.ants[mdi][1], self.ants[mdi][0]
                del self.ants[mdi]
            except:
                pass
        #    del self.ants
        for mdspw in range(len(self.GainBuffer) - 1, -1, -1):  # [::-1]:
            NA = len(self.GainBuffer[mdspw])
            for a in range(NA - 1, -1, -1):
                NP = len(self.GainBuffer[mdspw][a])
                for mdp in range(NP - 1, -1, -1):
                    del self.GainBuffer[mdspw][a][mdp]
                del self.GainBuffer[mdspw][a]
            del self.GainBuffer[mdspw]
        #    del self.GainBuffer

        for mdi in range(len(self.dt) - 1, -1, -1):  # [::-1]:
            del self.dt[mdi]
        #    del self.dt

        for mdi in range(len(self.dtArr) - 1, -1, -1):  # [::-1]:
            del self.dtArr[mdi]
        #    del self.dtArr

        for mdi in range(len(self.dtIdx) - 1, -1, -1):  # [::-1]:
            del self.dtIdx[mdi]
        #    del self.dtIdx

        for mdi in range(len(self.output) - 1, -1, -1):  # [::-1]:
            del self.output[mdi]
        #    del self.output

        for mdi in range(len(self.freqs) - 1, -1, -1):  # [::-1]:
            del self.freqs[mdi]

        for mdi in range(len(self.fittable) - 1, -1, -1):  # self.fittable[::-1]:
            del self.fittable[mdi]
        #    del self.fittable

        for mdi in range(len(self.wgtcorr) - 1, -1, -1):  # [::-1]:
            del self.wgtcorr[mdi]
        #    del self.wgtcorr

        for mdi in range(len(self.fittablebool) - 1, -1, -1):  # [::-1]:
            del self.fittablebool[mdi]
        #    del self.fittablebool

        for mdi in range(len(self.isGain) - 1, -1, -1):  # [::-1]:
            del self.isGain[mdi]
        #    del self.isGain

        for mdi in range(len(self.iscancoords) - 1, -1, -1):  # [::-1]:
            Npar = len(self.iscancoords[mdi])
            for mdp in range(Npar - 1, -1, -1):
                del self.iscancoords[mdi][mdp]
            del self.iscancoords[mdi]
        #    del self.iscancoords

    ############################################
    #
    #  CREATE INSTANCE
    #
    def __init__(self):
        """ Just the constructor of the 'modeler' class."""
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
        # Tells us if the fixedmodel array exists and should be used.
        self.removeFixed = False
        self.varfixed = []
        self.dpar = []
        self.Hessian = []
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
        self.Ccompmodel = None
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
        self.LorentzLine = lambda nu, nu0, P, G: P * 0.25 * G * G / (np.power(nu - nu0, 2.) + (0.5 * G)**2.)
        self.GaussLine = lambda nu, nu0, P, G: P * np.exp(-np.power((nu - nu0) / (G * KGaus), 2.))

        self.pieceWise = lambda t, p0, p1, t0, t1: np.clip(p0 + (p1 - p0) * (t - t0) / (t1 - t0), p0, p1)

        self.wgtEquation = lambda D, Kf: -D * Kf
        self.KfacWgt = 1.0

        # List of currently-supported model components:
        self.allowedmod = ['delta', 'Gaussian', 'disc', 'ring', 'sphere',
                           'bubble', 'expo', 'power-2', 'power-3', 'GaussianRing']
        self.resultstring = ''

    def _setup(self, model, parameters, fixedmodel, fixedparameters, scalefix, NCPU, only_flux, HankelOrder,
               isNumerical, useGains, gainFunction, isMixed):
        """ Setup the model components, fitting parameters, and gain functions. Compile the equations.

        Not to be called by the user."""

        import numpy as np

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
    def _compileAllModels(self):
        """ Compile all models (fixed, variable, and fixed-scale. Not to be called directly by the user. """
        self.resultstring = ''
        self.failed = False
        self._compileModel()
        self._compileFixedModel()
        self._compileScaleFixed()
        self._compileGains()

    #  print self.gainFunction
    def _compileGains(self):
        """ Compile the functions related to the antenna gains."""
        if self.isMixed:
            self.phaseAntsFunc = [lambda t, nu, p: 0. for i in range(self.Nants)]
            self.ampAntsFunc = [lambda t, nu, p: 1. for i in range(self.Nants)]
        else:
            self.phaseAntsFuncNu = [lambda nu, p: 0. for i in range(self.Nants)]
            self.ampAntsFuncNu = [lambda nu, p: 1. for i in range(self.Nants)]
            self.phaseAntsFuncT = [lambda t, p: 0. for i in range(self.Nants)]
            self.ampAntsFuncT = [lambda t, p: 1. for i in range(self.Nants)]

        self.parDependence = [[0] for i in range(self.Nants)]

        if self.isMixed:
            for ni in self.gainFunction[0].keys():
                tempstr = self.gainFunction[0][ni].replace(
                    'pieceWise(', 'self.pieceWise(t, ').replace(
                        't', 't[:,__import__(\'numpy\').newaxis]').replace(
                            'nu0', '%.12f' % self.freqs[0][0])
                modstr = 'self.phaseAntsFunc[' + str(ni) + '] = lambda t, nu, p: ' + tempstr

                parpos = [x.start() for x in re.finditer('p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                try:
                    exec(modstr, locals())
                except:
                    self.failed = True
                    self.resultstring = 'Syntax error in phase gain of antenna %i' % (ni)
                    return

            for ni in self.gainFunction[1].keys():
                tempstr = self.gainFunction[1][ni].replace(
                    'pieceWise(', 'self.pieceWise(t,').replace(
                        't', 't[:,__import__(\'numpy\').newaxis]').replace(
                            'nu0', '%.12f' % self.freqs[0][0])

                modstr = 'self.ampAntsFunc[' + str(ni) + '] = lambda t, nu, p: ' + tempstr

                parpos = [x.start() for x in re.finditer('p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                try:
                    exec(modstr, locals())
                except:
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

                parpos = [x.start() for x in re.finditer('p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                try:
                    exec(modstr, locals())
                except:
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

                parpos = [x.start() for x in re.finditer('p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                try:
                    exec(modstr, locals())
                except:
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

                parpos = [x.start() for x in re.finditer('p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                try:
                    exec(modstr, locals())
                except:
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

                parpos = [x.start() for x in re.finditer('p\[', tempstr)]
                for p0 in parpos:
                    p1 = tempstr[p0:].find(']') + p0
                    self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

                try:
                    exec(modstr, locals())
                except:
                    self.failed = True
                    self.resultstring = 'Syntax error in amp gain of antenna %i' % (ni)
                    return

        for ni in range(self.Nants):
            self.parDependence[ni] = np.unique(self.parDependence[ni]).astype(np.int32)

    def _compileModel(self):
        """ Compiles the variable model, according to the contents of the 'model' and 'parameters' lists."""

        # Define variable model:
        self.varfunc = [0.0 for component in self.model]

        for ii, component in enumerate(self.model):
            tempstr = self.var[ii].replace(
                'LorentzLine(', 'self.LorentLine(nu,').replace(
                    'GaussLine(', 'self.GaussLine(nu, ').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            if self.only_flux:
                try:
                    tempstr2 = tempstr.split(',')
                    if len(tempstr2) > 3:
                        self.strucvar.append(list(map(float, tempstr2[:2] + tempstr2[3:])))
                    else:
                        self.strucvar.append(list(map(float, tempstr2[:2])))
                except:
                    print(tempstr.split(','))
                    self.resultstring = '\n If only_flux=True, all variables but the flux must be constants! Aborting!'
                    self.failed = True

            modstr = 'self.varfunc[' + str(ii) + '] = lambda p, nu: [' + tempstr + ']'
            try:
                exec(modstr, locals())
            except:
                self.failed = True
                self.resultstring = 'Syntax error in component number %i of the variable model' % (ii)
                return

            if component not in self.allowedmod:
                self.resultstring = "\n Component '" + component + "' is unknown. Aborting!"
                self.failed = True
                return
            else:
                self.imod[ii] = self.allowedmod.index(component)

    def _compileFixedModel(self):
        """ Compiles the fixed model, according to the contents of the 'fixed' and 'fixedpars' lists."""

        if len(self.fixed) > 0 and 'model_column' in self.fixed:
            return

        # Define fixed model:
        import numpy as np
        self.fixedvarfunc = [0.0 for component in self.fixed]

        self.ifixmod = np.zeros(len(self.fixed), dtype=np.int32)  # []

        for ii, component in enumerate(self.fixed):
            tempstr = self.fixedvar[ii].replace(
                'LorentzLine(', 'self.LorentzLine(nu,').replace(
                    'GaussLine(', 'self.GaussLine(nu, ').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            modstr = 'self.fixedvarfunc[' + str(ii) + '] = lambda p, nu: [' + tempstr + ']'
            try:
                # if True:
                exec(modstr, locals())
            except:
                self.resultstring = 'Syntax error in component number %i of the fixed model' % (ii)
                self.failed = True
                return

            if component not in self.allowedmod:
                self.resultstring = "\n Component '" + component + "' is unknown. Aborting!"
                self.failed = True
                return
            else:
                self.ifixmod[ii] = self.allowedmod.index(component)

    def _compileScaleFixed(self):
        """ Compiles the scaling factor for the fixed model """

        tempstr = self.scalefix.replace(
            'LorentzLine(', 'self.LorentzLine(nu,').replace(
                'GaussLine(', 'self.GaussLine(nu, ').replace(
                    'nu0', '%.12f' % self.freqs[0][0])

        scalefixedstr = 'self._compiledScaleFixed = lambda p, nu: ' + tempstr + ' + 0.0'
        try:
            exec(scalefixedstr, locals())
        except:
            self.resultstring = 'Syntax error in the flux-scale equation'
            self.failed = True
            return

    ############################################
    #
    #  GET UNBOUND PARAMETER SPACE FROM BOUND PARAMETERS
    #
    def getPar2(self, mode=0):
        """ Function to change fitting parameters to/from the unbound space from/to the bound space.

        It also comptues the gradient for the Jacobian matrix. Unbound values are in index 0, gradient is
        in index 1, and bound values are in index 2. The equations for the changes of variables are taken
        from the MINUIT package."""

        if self.bounds is None:
            self.par2[1, :] = 1.0
            if mode == 0:
                self.par2[0, :] = self.par2[2, :]
            else:
                self.par2[2, :] = self.par2[0, :]
            return

        if mode == 0:
            p = self.par2[2, :]
            for p2 in range(len(p)):
                if self.bounds[p2] is None or (self.bounds[p2][0] is None and self.bounds[p2][1] is None):
                    self.par2[:, p2] = [p[p2], 1.0, p[p2]]
                elif self.bounds[p2][0] is None and self.bounds[p2][1] is not None:
                    self.par2[0, p2] = np.sqrt(np.power(self.bounds[p2][1] - p[p2] + 1, 2.) - 1.)
                elif self.bounds[p2][0] is not None and self.bounds[p2][1] is None:
                    self.par2[0, p2] = np.sqrt(np.power(p[p2] - self.bounds[p2][0] + 1, 2.) - 1.)
                else:
                    Kbound = (self.bounds[p2][1] - self.bounds[p2][0]) / 2.
                    self.par2[0, p2] = np.arcsin((p[p2] - self.bounds[p2][0]) / Kbound - 1.0)

        else:
            p = self.par2[0, :]
            for p2 in range(len(p)):
                if self.bounds[p2] is None or (self.bounds[p2][0] is None and self.bounds[p2][1] is None):
                    self.par2[2, p2] = p[p2]
                elif self.bounds[p2][0] is None and self.bounds[p2][1] is not None:
                    self.par2[2, p2] = self.bounds[p2][1] + 1. - np.sqrt(p[p2]**2. + 1.)
                elif self.bounds[p2][0] is not None and self.bounds[p2][1] is None:
                    self.par2[2, p2] = self.bounds[p2][0] - 1. + np.sqrt(p[p2]**2. + 1.)
                else:
                    Kbound = (self.bounds[p2][1] - self.bounds[p2][0]) / 2.
                    self.par2[2, p2] = self.bounds[p2][0] + Kbound * (np.sin(p[p2]) + 1.)

        for p2 in range(len(p)):
            if self.bounds[p2] is None or (self.bounds[p2][0] is None and self.bounds[p2][1] is None):
                self.par2[1, p2] = 1.0
            elif self.bounds[p2][0] is None and self.bounds[p2][1] is not None:
                self.par2[1, p2] = -self.par2[0, p2] / np.sqrt(self.par2[0, p2]**2. + 1.)
            elif self.bounds[p2][0] is not None and self.bounds[p2][1] is None:
                self.par2[1, p2] = self.par2[0, p2] / np.sqrt(self.par2[0, p2]**2. + 1.)
            else:
                Kbound = (self.bounds[p2][1] - self.bounds[p2][0]) / 2.
                self.par2[1, p2] = Kbound * np.cos(self.par2[0, p2])

    ############################################
    #
    #  LEVENBERG-MARQUARDT LEAST-SQUARE MINIMIZATION
    #
    def LMMin(self, p):
        """ Implementation of the Levenberg-Marquardt algorithm. Not to be called directly by the user. """

        NITER = int(self.LMtune[3] * len(p))
        self.calls = 0

        Lambda = float(self.LMtune[0])
        kfac = float(self.LMtune[1])
        functol = float(self.LMtune[2])
        if len(self.LMtune) > 4:
            partol = float(self.LMtune[4])
        else:
            partol = functol
        Chi2 = 0
        self.par2[2, :] = p
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
            print("\nComputing structure-only parameters")

            for midx in range(len(p)):
                if len(self.strucvar[midx]) > 3:
                    tempvar = self.strucvar[midx][:2] + [p[midx]] + self.strucvar[midx][2:]
                else:
                    tempvar = self.strucvar[midx] + [p[midx]]
                for i in range(len(tempvar)):
                    for j in range(len(p) + 1):
                        self.varbuffer[j][midx, i, :nnu] = tempvar[i]

        CurrChi = self.residuals(self.par2[2, :], -1, dof=False)
        Hessian2[:, :] = self.Hessian * (self.par2[1, :])[np.newaxis, :] * (self.par2[1, :])[:, np.newaxis]
        Gradient2[:] = self.Gradient * self.par2[1, :]
        backupHess[:, :] = Hessian2
        backupGrad[:] = Gradient2

        controlIter = 0
        for i in range(NITER):
            controlIter += 1
            for n in range(len(p)):
                HessianDiag[n, n] = Hessian2[n, n]
            try:
                goodsol = True
                Inverse[:] = np.linalg.pinv(Hessian2 + Lambda * HessianDiag)
                Dpar = np.dot(Inverse, Gradient2)
                DirDer = sum([Hessian2[n, n] * Dpar[n] * Dpar[n] for n in range(len(p))])
                DirDer2 = np.sum(Gradient2 * Dpar)
                TheorImpr = DirDer - 2. * DirDer2
            except:
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
            sys.stdout.write("\n REACHED MAXIMUM NUMBER OF ITERATIONS!\n" +
                             "The algorithm may not have converged!\n" +
                             "Please, check if the parameter values are meaningful.\n")
            sys.stdout.flush()

        self.getPar2(mode=1)

        try:
            return [self.par2[2, :], np.linalg.pinv(self.Hessian), Chi2]
        except:
            return False

    ############################################
    #
    #  NON-ALGEBRAIC MODELS
    #
    # Compute elements of Taylor expansion of the source's Hankel transform:
    def gridModel(self, imod, tempvar):
        """ Compute elements of Taylor expansion of the source's Hankel transform."""

        n = self.HankelOrder - 1

        if imod == 'GaussianRing':   # Gaussian Ring

            a = 1.0
            k = 1. / (2. * np.power(tempvar[6], 2.)) * tempvar[3] * tempvar[3] / 4.
            m = 2 * n + 1

            # Recurrence relation:
            merf = (1. - spec.erf(-np.sqrt(k) * a))

            m0 = np.sqrt(np.pi / k) / 2. * merf
            m1 = np.exp(-k * a * a) / 2. / k + np.sqrt(np.pi / k) * a / 2. * merf

            tempvar.append(np.ones(np.shape(m1)))  # term for n=0, normalized

            if m in [0, 1]:
                return tempvar  # Only order n=0.

            else:  # Higher orders.
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
                        tempvar.append(
                            res / m1 * np.power(-1., (mi - 1) / 2) / np.power(np.math.factorial((mi - 1) / 2),
                                                                              2.))

            return tempvar

        else:

            raise ValueError("\n\nModel %i was not correctly interpreted!\n\n" % imod)

    ############################################
    #
    #  COMPUTE RESIDUALS FOR A MODEL REALIZATION
    #
    def residuals(self, p, mode=-1, dof=True):
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
        for j in range(len(p)):
            #  self.dpar[j] = self.minnum
            if np.abs(p[j]) < 1.0:
                self.dpar[j] = self.minnum
            else:
                self.dpar[j] = np.abs(p[j]) * self.minnum

        # Compute antenna-gain corrections:
        if len(self.useGains) > 0 and mode != 0:
            for spw in spwrange:
                Nchan = len(self.freqs[spw])
                Nint = len(self.dtArr[spw])
                for i in self.useGains:
                    for j in range(len(self.parDependence[i])):
                        ptemp = [pi for pi in p]
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
            isfixed = True
            modbackup = self.imod[0]
            for spw in spwrange:
                self.output[spw][:] = 0.0
                for midx, mi in enumerate(self.ifixmod):
                    tempvar = self.fixedvarfunc[midx](p, self.freqs[spw])
                    self.imod[0] = self.ifixmod[midx]
                    if self.imod[0] in self.isNumerical:
                        tempvar = self.gridModel(self.imod[0], tempvar)
                    for i in range(len(tempvar)):
                        nnu = len(self.freqs[spw])
                        self.varbuffer[0][0, i, :nnu] = tempvar[i]
                    self.Ccompmodel(spw, nui, 0)

            self.imod[0] = modbackup
            return 0

        else:  # Compute the variable model:
            currmod = self.model
            currvar = self.varfunc

            currimod = self.imod
            isfixed = False

        ChiSq = 0.0
        ndata = 0
        for spw in spwrange:
            nt, nnu = np.shape(self.output[spw])
            if nui == -1:
                scalefx = self._compiledScaleFixed(p, self.freqs[spw])
                self.varfixed[0][:nnu] = scalefx
                for j in range(len(p)):
                    ptemp = [pi for pi in p]
                    ptemp[j] += self.dpar[j]
                    # Variables of current component
                    scalefx = self._compiledScaleFixed(ptemp, self.freqs[spw])
                    self.varfixed[j + 1][:nnu] = scalefx
            else:
                scalefx = self._compiledScaleFixed(p, self.freqs[spw][nui])
                self.varfixed[0][0] = scalefx
                for j in range(len(p)):
                    ptemp = [pi for pi in p]
                    ptemp[j] += self.dpar[j]
                    # Variables of current component
                    scalefx = self._compiledScaleFixed(ptemp, self.freqs[spw][nui])
                    self.varfixed[j + 1][0] = scalefx

            for midx, modi in enumerate(currmod):
                ptemp = [float(pi) for pi in p]

                if nui == -1:  # Continuum

                    if self.only_flux:
                        currflux = p[midx]
                        nstrucpars = len(self.strucvar[midx])
                        self.varbuffer[0][midx, 2, :nnu] = currflux

                    else:
                        # Variables of current component
                        tempvar = currvar[midx](p, self.freqs[spw])
                        if modi in self.isNumerical:
                            tempvar = self.gridModel(modi, tempvar)

                        for i in range(len(tempvar)):
                            self.varbuffer[0][midx, i, :nnu] = tempvar[i]

                    if self.only_flux:
                        ptemp[midx] += self.dpar[midx]
                        self.varbuffer[midx + 1][midx, 2, :nnu] = ptemp[midx]

                    else:

                        for j in range(len(p)):
                            ptemp = [pi for pi in p]
                            ptemp[j] += self.dpar[j]
                        # Variables of current component
                            tempvar = currvar[midx](ptemp, self.freqs[spw])
                            if modi in self.isNumerical:
                                tempvar = self.gridModel(modi, tempvar)
                            for i in range(len(tempvar)):
                                self.varbuffer[j + 1][midx, i, :nnu] = tempvar[i]

                else:  # Spectral mode

                    # Variables of current component
                    tempvar = currvar[midx](p, self.freqs[spw][nui])
                    if modi in self.isNumerical:
                        tempvar = self.gridModel(modi, tempvar)
                    for i in range(len(tempvar)):
                        self.varbuffer[0][midx, i, 0] = tempvar[i]

                    for j in range(len(p)):
                        ptemp = [pi for pi in p]
                        ptemp[j] += self.dpar[j]  # self.minnum
                        # Variables of current component
                        tempvar = currvar[midx](ptemp, self.freqs[spw][nui])
                        if modi in self.isNumerical:
                            tempvar = self.gridModel(modi, tempvar)
                        for i in range(len(tempvar)):
                            self.varbuffer[j + 1][midx, i, 0] = tempvar[i]

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
                sys.stdout.write("\r Iteration # %i. " % (self.calls))
                sys.stdout.write("\r \t\t\t\t Achieved ChiSq:  %.8e" % (ChiSq))
                sys.stdout.flush()

        if ChiSq <= 0.0:
            raise ValueError("Invalid Chi Square!" +
                             "Maybe the current fitted value of flux (and/or size) is negative!" +
                             "Please, set BOUNDS to the fit!")

        if mode in [-1, -2, -3]:
            if dof:
                self.calls = 0
                return [ChiSq, ndata]
            else:
                return ChiSq

    ############################################
    #
    #  COMPUTE QUI SQUARE FOR A MODEL REALIZATION (NO DERIVATIVES)
    #
    def ChiSquare(self, p, bounds=None, p_ini=[]):
        """ Just a wrapper of the 'residuals' function, usually called by simplex."""

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


############################################
#
#  HELPER FUNCTIONS
#
def channeler(spw, width=1, maxchans=[3840, 3840, 3840, 3840]):
    """ Function to convert a string with spw selection into lists of channels to select/average.

    It follows the CASA syntax."""

    if spw == '':
        spw = ','.join(list(map(str, range(len(maxchans)))))

    entries = spw.split(',')
    output = [[] for i in maxchans]

    for entry in entries:
        check = entry.split(':')

        spws = list(map(int, check[0].split('~')))
        if len(spws) == 1:
            selspw = [spws[0]]
        else:
            selspw = range(spws[0], spws[1] + 1)

        for sp in selspw:

            if sp + 1 > len(maxchans):
                errstr = "There are only %i spw in the data!\n Please, revise the 'spw' parameter" % (len(maxchans))
                return [False, errstr]

            if len(check) == 1:
                chranges = ['0~' + str(maxchans[sp] - 1)]
            else:
                chans = check[1]
                chranges = chans.split(';')
            ranges = []
            for chran in chranges:
                ch1, ch2 = list(map(int, chran.split('~')))
                if ch1 > ch2:
                    errstr = '%i is larger than %i. Revise channels for spw %i' % (ch1, ch2, sp)
                    return [False, errstr]
                ch2 = min([ch2, maxchans[sp] - 1])
                for i in range(ch1, ch2 + 1, width):
                    ranges.append(range(i, min([(i + width), ch2 + 1])))

            output[sp] = ranges

    return [True, output]


def modelFromClean(imname, ichan=0, echan=0):
    """Reads the '*.model' data created by task 'clean' and returns lists compatible with UVMultiFit.

    The user can use this function to define the 'model' and  'var' lists (or the 'fixed' and
    'fixedvar' lists) for UVMultiFit.

    Parameters
    ----------
    ichan : `int`
      First channel of the image (cube) to read.
    echan : `int`
      Last channel to read.

    Returns
    -------
    Model_components : `list`
      The flux densities that the task will return for each pixel will be the *frequency average*
      between ichan and echan."""

    # Read model:
    success = ia.open(str(imname))
    if not success:
        errstr = 'Could not open ', str(imname)
        return [False, errstr]

    echan += 1

    try:
        modarray = np.average(ia.getchunk()[:, :, ichan:echan, 0], axis=2)
        deltaxy = ia.summary()['incr']
        refpix = ia.summary()['refpix']
        cleanlocs = np.where(modarray != 0.0)
        cleans = np.transpose([cleanlocs[0], cleanlocs[1]] + [modarray[cleanlocs]])
        totflux = np.sum(cleans[:, 2])

        # Put shifts in arcseconds w.r.t. image center:
        cleans[:, 0] = (cleans[:, 0] - refpix[0]) * deltaxy[0] * 180. / np.pi * 3600.
        cleans[:, 1] = (cleans[:, 1] - refpix[1]) * deltaxy[1] * 180. / np.pi * 3600.

        # Set model as a set of deltas:
        model = ['delta' for cl in cleans]
        var = ['%.3e, %.3e, %.3e' % tuple(cl) for cl in cleans]

        # Send results:
        ia.close()
        return [True, model, var, [cleanlocs[0], cleanlocs[1]]]
    except:
        ia.close()
        errstr = 'Problem processing the CASA image. Are you sure this is an image (and not an image cube)?'
        return [False, errstr]

######################################
# THE FOLLOWING CODE IS A MODIFICATION OF THE SIMPLEX ALGORITHM IMPLEMENTED
# IN SCIPY. THIS NEW CODE ALLOWS US TO SET DIFFERENT ERRORS IN THE DIFFERENT
# PARAMETERS, INSTEAD OF ONE ERROR FOR ALL, WHICH MAY OF COURSE BE MUCH
# AFFECTED BY THE DIFFERENT POSSIBLE ORDERS OF MAGNITUDE IN THE PARAMETERS).
######################################

# Copyright (c) 2001, 2002 Enthought, Inc.
# All rights reserved.
#
# Copyright (c) 2003-2012 SciPy Developers.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# a. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# b. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
# c. Neither the name of the author nor the names of contributors may
#    be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


def wrap_function(function, args):
    ncalls = [0]

    def function_wrapper(x):
        ncalls[0] += 1
        return function(x, *args)
    return ncalls, function_wrapper


def _mod_simplex(func, x0, args=(), callback=None, relstep=1.e-1,
                 relxtol=1.e-4, relftol=1.e-3, maxiter=None, maxfev=None,
                 disp=False, return_all=False,
                 **unknown_options):
    """
    Minimization of scalar function of one or more variables using the
    Nelder-Mead algorithm.

    Options for the Nelder-Mead algorithm are:
        disp : bool
            Set to True to print("onvergence messages")
       relxtol : float
            Relative error in solution `xopt` acceptable for convergence.
       relftol : float
            Relative error in ``fun(xopt)`` acceptable for convergence.
       relstep : float
            Relative size of the first step in building the simplex.
        maxiter : int
            Maximum number of iterations to perform.
        maxfev : int
            Maximum number of function evaluations to make.

    """

    maxfun = maxfev
    retall = return_all

    fcalls, func = wrap_function(func, args)
    x0 = np.asfarray(x0).flatten()
    N = len(x0)
    rank = len(x0.shape)
    if not -1 < rank < 2:
        raise ValueError("Initial guess must be a scalar or rank-1 sequence.")
    if maxiter is None:
        maxiter = N * 200
    if maxfun is None:
        maxfun = N * 200

    rho = 1
    chi = 2
    psi = 0.5
    sigma = 0.5
    one2np1 = list(range(1, N + 1))

    if rank == 0:
        sim = np.zeros((N + 1,), dtype=x0.dtype)
    else:
        sim = np.zeros((N + 1, N), dtype=x0.dtype)
    fsim = np.zeros((N + 1,), float)
    sim[0] = x0
    if retall:
        allvecs = [sim[0]]
    fsim[0] = func(x0)
    nonzdelt = relstep
    zdelt = 0.00025
    for k in range(0, N):
        y = np.array(x0, copy=True)
        if y[k] != 0:
            y[k] = (1 + nonzdelt) * y[k]
        else:
            y[k] = zdelt

        sim[k + 1] = y
        f = func(y)
        fsim[k + 1] = f

    ind = np.argsort(fsim)
    fsim = np.take(fsim, ind, 0)
    # sort so sim[0, :] has the lowest function value
    sim = np.take(sim, ind, 0)

    iterations = 1

    minx = np.sqrt(np.finfo(np.float64).eps)
    if isinstance(relxtol, float):
        relxtol = relxtol * np.ones(len(x0))
    else:
        relxtol = np.array(relxtol)

    while (fcalls[0] < maxfun and iterations < maxiter):
        testsim = np.copy(sim[0])
        testsim[np.abs(testsim) < minx] = minx
        testfsim = max(fsim[0], minx)
        deltas = np.abs(sim[1:] - sim[0])
        smaller = len(deltas[deltas > relxtol])

        if (smaller == 0 and np.max(np.abs((fsim[0] - fsim[1:]) / testfsim)) <= relftol):
            break

        xbar = np.add.reduce(sim[:-1], 0) / N
        xr = (1 + rho) * xbar - rho * sim[-1]
        fxr = func(xr)
        doshrink = 0

        if fxr < fsim[0]:
            xe = (1 + rho * chi) * xbar - rho * chi * sim[-1]
            fxe = func(xe)

            if fxe < fxr:
                sim[-1] = xe
                fsim[-1] = fxe
            else:
                sim[-1] = xr
                fsim[-1] = fxr
        else:  # fsim[0] <= fxr
            if fxr < fsim[-2]:
                sim[-1] = xr
                fsim[-1] = fxr
            else:  # fxr >= fsim[-2]
                # Perform contraction
                if fxr < fsim[-1]:
                    xc = (1 + psi * rho) * xbar - psi * rho * sim[-1]
                    fxc = func(xc)

                    if fxc <= fxr:
                        sim[-1] = xc
                        fsim[-1] = fxc
                    else:
                        doshrink = 1
                else:
                    # Perform an inside contraction
                    xcc = (1 - psi) * xbar + psi * sim[-1]
                    fxcc = func(xcc)

                    if fxcc < fsim[-1]:
                        sim[-1] = xcc
                        fsim[-1] = fxcc
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        sim[j] = sim[0] + sigma * (sim[j] - sim[0])
                        print('simj', sim[j])
                        fsim[j] = func(sim[j])

        ind = np.argsort(fsim)
        sim = np.take(sim, ind, 0)
        fsim = np.take(fsim, ind, 0)
        if callback is not None:
            callback(sim[0])
        iterations += 1
        if retall:
            allvecs.append(sim[0])

    x = sim[0]
    fval = np.min(fsim)
    warnflag = 0

    return [x, fval]


############################################
#
#  UVMULTIFIT WORKING ON DIRTY IMAGES: IMMULTIFIT
#
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

    def __init__(self, parent=None, reinvert=True, start=0, nchan=-1, psf='', residual='', dBcut=-30., **kwargs):
        """ Constructor."""

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

        self._printInfo("\nThere is no measurement set (i.e., you are running 'immultifit').\n" +
                        "Will generate an image, instead.\n")
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
        npix = float(imdims[0] * imdims[1])
        temparr = np.zeros((imdims[0], imdims[1]), dtype=np.complex128)
        for i in range(0, self.start):
            resim[:, :, :, i] = 0.0
        for i in range(self.nchan, imdims[-1]):
            resim[:, :, :, i] = 0.0

        for i in range(self.start, self.nchan):
            self._printInfo("\r Doing frequency channel %i of %i" % (i + 1 - self.start, self.nchan - self.start))
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

        self._printInfo("\n\n Will now save the unconvolved model image.\n")
        modname = '.'.join(self.residual.split('.')[:-1]) + '.fitmodel.immultifit'
        if os.path.exists(modname):
            shutil.rmtree(modname)
        os.system('cp -r %s %s' % (self.residual, modname))

        ia.open(modname)
        resim = ia.getchunk()
        imdims = np.shape(resim)
        npix = float(imdims[0] * imdims[1])
        temparr = np.zeros((imdims[0], imdims[1]), dtype=np.complex128)
        for i in range(self.start, self.nchan):
            self._printInfo("\r Doing frequency channel %i of %i" % (i + 1 - self.start, self.nchan - self.start))
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
    def checkInputs(self, data_changed=False):
        """ Function re-definition for the immultifit class."""

        self._printInfo("\nIn image mode, the whole image is taken.\n" +
                        "Will override the 'spw', 'field', 'MJDrange', and 'scan' parameters.\n\n")
        self._printInfo("In image mode, the channelization used is that of the image.\n" +
                        "Will override the 'chanwidth' and 'timewidth' parameters.\n")

        self.savemodel = True
        success = self._checkOrdinaryInputs()

        if not success:
            return False

        try:
            self.stokes = abs(int(self.stokes))
        except:
            self._printInfo("\n\nIn image mode, 'stokes' must be an integer\n(i.e., the Stokes column of the image).\n")
            polimag = ['I', 'Q', 'U', 'V']
            if self.stokes in polimag:
                stchan = polimag.index(self.stokes)
                self._printInfo("Since stokes is '%s', will try with image column %i\n" % (self.stokes, stchan))
                self.stokes = stchan
            else:
                self._printError("\nPlease, set the 'stokes' keyword to an image pol. channel\n")
                return False

        try:
            ia.open(self.residual)
            imdims = np.array(ia.shape())
            imsum = ia.summary()
            ia.close()
        except:
            self._printError("The residual image is not an image (or does not exist)!\n\n")
            return False

        try:
            ia.open(self.psf)
            imdims2 = np.array(ia.shape())
            ia.close()
        except:
            self._printError("The PSF image is not an image (or does not exist)!\n\n")
            return False

        if max(np.abs(imdims - imdims2)) > 0:
            self._printError("The dimensions of the PSF image and the image of residuals do NOT coindice!\n\n")
            return False

        if imdims[2] > self.stokes:
            self._printInfo("Selecting stokes column %i" % self.stokes)
        else:
            self._printError("\n\nThe images only have %i stokes columns, but you selected stokes=%i\n" %
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
#       except:
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
    def readData(self, data_changed=True, del_data=True):
        """ Function redefinition for the immultifit class."""
        if self.reinvert:

            os.system('rm -rf %s.fft.*' % self.residual)
            os.system('rm -rf %s.fft.*' % self.psf)

            self._printInfo("\nInverting images into Fourier space\n")
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
        outpre = np.zeros((npix2, nnu), dtype=np.float64)
        outpim = np.zeros((npix2, nnu), dtype=np.float64)
        fixedre = np.zeros((npix2, nnu), dtype=np.float64)
        fixedim = np.zeros((npix2, nnu), dtype=np.float64)
        freqs = np.zeros(nnu, dtype=np.float64)

        zeros = []
        for i in range(imdims[3]):
            freqs[i] = ia.toworld([0, 0, self.stokes, i + self.start])['numeric'][3]

        self._printInfo("\nReading gridded UV coordinates and weights\n")
        u0, v0, s0, nu0 = ia.toworld([0, 0, self.stokes, self.start])['numeric']
        u1, v1, s0, nu1 = ia.toworld([1, 1, self.stokes, self.start])['numeric']
        du = u1 - u0
        dv = v1 - v0
        for j in range(imdims[0]):
            self._printInfo("\r Reading row %i of %i" % (j, imdims[0]))
            for k in range(imdims[1]):
                cpix = imdims[1] * j + k
                ui[cpix] = (u0 + du * j)
                vi[cpix] = (v0 + dv * k)
                wgti[cpix, :] = np.copy(wgt[j, k, self.stokes, :])

        del wgt
        ia.close()
        # Maximum dynamic range set by "Window effect"
        maskwgt = wgti < np.max(wgti) * (10.**(self.dBcut / 10.))
        self._printInfo("\n\nReading gridded visibilities\n")
        ia.open(self.residual + '.fft.real')
        datas = ia.getchunk()
        for j in range(imdims[0]):
            self._printInfo("\r Reading row %i of %i" % (j, imdims[0]))
            for k in range(imdims[1]):
                cpix = imdims[1] * j + k
                datare[cpix, :] = datas[j, k, self.stokes, self.start:self.nchan]

        del datas
        ia.close()

        self._printInfo(" ")
        ia.open(self.residual + '.fft.imag')
        datas = ia.getchunk()
        for j in range(imdims[0]):
            self._printInfo("\r Reading row %i of %i" % (j, imdims[0]))
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

        self._printInfo("\nDone reading.\n")

        self.initData()

        return True
