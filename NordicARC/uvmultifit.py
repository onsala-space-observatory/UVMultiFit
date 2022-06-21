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
>>>            NCPU=4, pbeam=False, ldfac = 1.22, dish_diameter=0.0,
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
__revision__ = "$Id: 3.0.0-p2 2019-03-28 15:00:00 marti-vidal $"
__docformat__ = 'reStructuredText'

import logging
import re
import gc
import os
import time

import numpy as np

from casatools import ms
from casatools import table
from casatools import coordsys
from casatools import image

import uvmultimodel as uvmod
from .modeler import modeler
# global ms, tb

__version__ = "3.0-p2"
date = 'JUN 2022'

print("C++ shared library loaded successfully\n")

tb = table()
ms = ms()
ia = image()
cs = coordsys()

# if True:
#     from taskinit import gentools
#     from clearcal_cli import clearcal_cli as clearcal
#     ms = gentools(['ms'])[0]
#     tb = gentools(['tb'])[0]
#     ia = gentools(['ia'])[0]
#     cs = gentools(['cs'])[0]

greetings = '#######################################################################\n'
greetings += '# UVMULTIFIT -- ' + date + '. EUROPEAN ALMA REGIONAL CENTER (NORDIC NODE) #\n'
greetings += '#      Please, add the UVMULTIFIT reference to your publications:     #\n'
greetings += '#     Marti-Vidal, Vlemmings, Muller & Casey 2014, A&A, 563, A136     #\n'
greetings += '#######################################################################\n'


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
    # def deleteData(self, delmodel=True):
    #     """Delete pointers to the data.
    #     Hopefully, this will release memory when gc.collect() is run."""
    #
    #     for i in range(len(self.averdata) - 1, -1, -1):
    #         del self.averdata[i]
    #         if self.takeModel:
    #             try:
    #                 del self.avermod[i]
    #             except Exception:
    #                 pass
    #         del self.averweights[i]
    #         del self.averfreqs[i]
    #         del self.u[i]
    #         del self.v[i]
    #         del self.w[i]
    #         del self.t[i]
    #         del self.tArr[i]
    #         del self.tIdx[i]
    #         del self.ant1[i]
    #         del self.ant2[i]
    #         del self.RAshift[i]
    #         del self.Decshift[i]
    #         del self.Stretch[i]
    #
    #     try:
    #         for mdi in self.iscancoords[::-1]:
    #             Npar = len(mdi)
    #             for mdp in range(Npar - 1, -1, -1):
    #                 del mdi[mdp]
    #             del mdi
    #         del self.iscancoords
    #     except Exception:
    #         pass
    #
    #     del self.averdata, self.avermod, self.averweights, self.averfreqs, self.v, self.u, self.w, self.t
    #     del self.ant1, self.ant2, self.RAshift, self.Decshift, self.Stretch
    #
    #     if delmodel:
    #         self.mymodel.deleteData()
    #
    # def __del__(self):
    #     # Delete the model first, so that C++ has green light to free pointers:
    #     # del self.mymodel
    #     #
    #     # # Clear all data and C++ pointers:
    #     # self.deleteData(delmodel=False)
    #     # uvmod.clearPointers(0)
    #     # uvmod.clearPointers(1)
    #     pass

    ############################################
    #
    #  CREATE INSTANCE
    #
    def __init__(self, uniform=False,
                 chanwidth=1, timewidth=1, stokes='I', write='', MJDrange=[-1.0, -1.0], ldfac=1.22,
                 phase_center='',
                 fixed=[], fixedvar=[], scalefix='1.0', outfile='modelfit.dat', NCPU=4, pbeam=False,
                 dish_diameter=0.0, cov_return=False, finetune=False, uvtaper=0.0,
                 method='levenberg', wgt_power=1.0,
                 LMtune=[1.e-3, 10., 1.e-5, 200, 1.e-3], SMPtune=[1.e-4, 1.e-1, 200], only_flux=False,
                 proper_motion=0.0, HankelOrder=80, phase_gains={}, amp_gains={}):
        """Just the constructor method, for class instantiation."""
        logging.basicConfig(level=logging.INFO,
                            format='%(name)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger("uvmultifit")
        self.logger.debug("uvmultifit::__init__")
        self._printInfo(greetings)
        # self._printWarning("beware")
        # self._printError("oops!")
        self.first_time = True
        uvmod.clearPointers(2)

        self.implemented = ['delta', 'disc', 'ring', 'Gaussian', 'sphere',
                            'bubble', 'expo', 'power-2', 'power-3', 'GaussianRing']
        # Number of variables for the models
        numvars = [3, 6, 6, 6, 6, 6, 6, 6, 6, 7]
        # Models that need a series expansion for J0.
        self.isNumerical = ['GaussianRing']
        self.uniform = uniform
        self.averdata = []
        self.avermod = []
        self.averweights = []
        self.averfreqs = []
        self.u = []
        self.v = []
        self.fixed = fixed
        self.wgt_power = wgt_power
        self.fixedvar = fixedvar
        self.scalefix = scalefix
        self.outfile = outfile
        self.ranges = []
        self.chanwidth = chanwidth
        #  self.timewidth = timewidth
        self.stokes = stokes
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
            self._printError("keyword 'write' should be set to either '', 'model', 'residuals' or 'calibrated'")

        #######################
        # TIMEWIDTH IS NOT CURRENTLY USED (ISSUES WITH NEW MS TOOL):
        self.timewidth = 1
        if timewidth != 1:
            self._printWarning("timewdith>1 cannot be currently set, due to issues with new MS tool")

        # Start instance:
        # self._startUp()

    def __repr__(self):
        txt = "uvmultifit("
        txt += f"vis = '{self.vis}', "
        txt += f"spw = {self.spw}, "
        txt += f"column = '{self.column}', "
        txt += f"field = {self.field}, "
        txt += f"scans = {self.scans}, "
        txt += f"uniform = {self.uniform}, "
        txt += f"chanwidth = {self.chanwidth}, "
        txt += f"timewidth = {self.timewidth}, "
        txt += f"stokes = '{self.stokes}', "
        txt += f"write_model = {self.write_model}, "
        txt += f"MJDrange = {self.MJDrange}, "
        txt += f"ldfac = {self.ldfac}, "
        txt += f"model = {self.model}, "
        txt += f"var = {self.var}, "
        txt += f"p_ini = {self.p_ini}, "
        txt += f"phase_center = '{self.phase_center}', "
        txt += f"fixed = {self.fixed}, "
        txt += f"fixedvar = {self.fixedvar}, "
        txt += f"scalefix = {self.scalefix}, "
        txt += f"outfile = {self.outfile}, "
        txt += f"NCPU = {self.NCPU}, "
        txt += f"pbeam = {self.pbeam}, "
        txt += f"dish_diameter = {self.dish_diameter}, "
        txt += f"OneFitPerChannel = {self.OneFitPerChannel}, "
        txt += f"bounds = {self.bounds}, "
        txt += f"cov_return = {self.cov_return}, "
        txt += f"finetune = {self.finetune}, "
        txt += f"uvtaper = {self.uvtaper}, "
        txt += f"method = '{self.method}', "
        txt += f"wgt_power = {self.wgt_power}, "
        txt += f"LMtune = {self.LMtune}, "
        txt += f"SMPtune = {self.SMPtune}, "
        txt += f"only_flux = {self.only_flux}, "
        txt += f"proper_motion = {self.proper_motion}, "
        txt += f"HankelOrder = {self.HankelOrder}, "
        txt += f"phase_gains = {self.phase_gains}, "
        txt += f"amp_gains = {self.amp_gains})"
        return txt

    def select_data(self, vis, spw='0', column='data', field=0, scans=[]):
        self.vis = vis
        self.spw = spw
        self.column = column
        self.field = field
        self.scans = scans

    def select_model(self, model=['delta'], var=['p[0], p[1], p[2]'], p_ini=[0.0, 0.0, 1.0],
                     bounds=None, OneFitPerChannel=False):
        for i, m in enumerate(model):
            self._printInfo("model[%d]: %s" % (i, m))

        if isinstance(var, list):
            for i, v in enumerate(var):
                self._printInfo("variables[%d]: %s" % (i, v))
        else:
            self._printInfo("variables[0]: %s" % (var))

        for i, p in enumerate(p_ini):
            self._printInfo("initial p[%d]: %s" % (i, p))

        self.model = model
        self.var = var
        self.p_ini = p_ini
        self.bounds = bounds
        self.OneFitPerChannel = OneFitPerChannel

    ############################################
    #
    #  STARTERS
    #
    # This method will be overriden in the GUI, to avoid execution of CheckInputs() and the fit:
    def start_fit(self):
        """This is run each time the class is instantiated."""
        # if not goodclib:
        #     self._printError("C++ library cannot be loaded! Please, contact the Nordic ARC node.")
        #     return False

        self.logger.debug("uvmultifit::_startup")
        # self._printInfo(greetings)

        self.mymodel = modeler()
        self.mymodel.Ccompmodel = uvmod.modelcomp

        # Check parameters and read the data in:
        if not self.checkInputs():
            self._printError("aborting UVMultiFit, checkInputs failed!")
            return False

        if not self.readData(del_data=False):
            self._printError("aborting UVMultiFit, readData failed!")
            return False

            # Compile Model:
        if not self.initData():
            self._printError("aborting UVMultiFit, initData failed!")
            return False

        if not self.initModel():
            self._printError("aborting UVMultiFit, initModel failed!")
            return False

        if not self.finetune:
            if not self.fit():
                self._printError("failed fit!\n")
                return False
            if self.write_model in [1, 2, 3]:
                if self.timewidth == 1 and self.stokes not in ['Q', 'U', 'V']:
                    self._printInfo("writing into measurement set(s)")
                    self.writeModel()
                else:
                    msg = "cannot write into measurement set!\n"
                    msg += "If you want to fill-in the model (or corrected) column:\n"
                    msg += "    1.- 'timewidth' and 'chanwidth' should be set to 1\n"
                    msg += "    2.- 'stokes' should be set to either 'I', 'PI', or a corr. product."
                    self._printError(msg)

            self._printInfo("fit done!! And UVMULTIFIT class instantiated successfully!")
        else:
            self._printInfo("UVMultiFit class successfully instantiated")
        return True

    ############################################
    #
    #  PRINT MESSAGES AND ERRORS
    #
    # Functions overriden in GUI mode:
    def _printError(self, message):
        """Prints a message and raises an exception."""
        for part in message.split('\n'):
            if len(part) > 0:
                self.logger.critical(part)
        raise Exception(message)

    def _printWarning(self, message):
        """Prints a warning message."""
        self.logger.warning(message)

    def _printInfo(self, message):
        """Prints a message."""
        for part in message.split('\n'):
            if len(part) > 0:
                self.logger.info(part)

    def _printDebug(self, message):
        """Prints a debug message."""
        self.logger.debug(message)

    ############################################
    #
    #  WRITE MODEL (OR RESIDUALS)
    #
    def writeModel(self):
        """Writes the requested information into the measurement sets.

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

        self._printDebug("uvmultifit::writeModel")
        self._printWarning("writing to mosaics is experimental and may not work!\n")

        for v in self.vis:
            # Get the columns of parallel-hand correlations:
            success = ms.open(v)
            if not success:
                self._printError("%s cannot be openned in write mode" % (v))
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
                        self._printInfo("Doing %s: spw %i, scan_id %i" % (v, sp, scan))
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

            self._printInfo("%s written successfully" % (column))

        return True

    ############################################
    #
    #  SANITY CHECKS AND DEFINING BASIC PARAMETERS
    #
    def _checkOrdinaryInputs(self):
        """Performs some sanity checks on the input parameters.
      This function should not be called directly by the user."""

        self._printDebug("uvmultifit::_checkOrdinaryInputs")
        if isinstance(self.model, str):
            self.model = list([self.model])
        else:
            self.model = list(self.model)

        # print(self.var)
        if isinstance(self.var, str):
            self.var = list([self.var])
        else:
            self.var = list(self.var)

        if isinstance(self.fixedvar, (float, str)):
            self.fixedvar = list([self.fixedvar])
        else:
            self.fixedvar = list(self.fixedvar)

        if isinstance(self.fixed, str):
            self.fixed = list([self.fixed])
        else:
            self.fixed = list(self.fixed)

        if isinstance(self.p_ini, (float, str)):
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
            self._printError("'phase_gains' must be a dictionary!")

        if not isinstance(self.amp_gains, dict):
            self._printError("'amp_gains' must be a dictionary!")

        for key in self.phase_gains.keys():
            if isinstance(key, int):
                self.isMixed = True
                break

        for key in self.amp_gains.keys():
            if isinstance(key, int):
                self.isMixed = True
                # self._printError("Inconsistent 'amp_gains' and 'phase_gains'!")
                break

        for key in self.phase_gains.keys():
            if key == 'nuG':
                if self.isMixed:
                    self._printError(
                        "You cannot define split gains and mixed gains in the same fit. It makes no sense!\n")
                self.phase_gainsNu = self.phase_gains['nuG']
            elif key == 'tG':
                if self.isMixed:
                    self._printError(
                        "You cannot define split gains and mixed gains in the same fit. It makes no sense!\n")
                self.phase_gainsT = self.phase_gains['tG']
            elif not isinstance(key, int):
                self._printError("The keys of 'phase_gains' must be integers or 'nuG/tG'!")
            else:
                self.phase_gainsNuT[key] = self.phase_gains[key]

        for key in self.phase_gainsNu.keys():
            if not isinstance(key, int):
                self._printError("The keys of 'phase_gains[nuG]' must be integers!")
        for key in self.phase_gainsT.keys():
            if not isinstance(key, int):
                self._printError("The keys of 'phase_gains[tG]' must be integers!")

        for key in self.amp_gains.keys():
            if key == 'nuG':
                if self.isMixed:
                    self._printError(
                        "You cannot define split gains and mixed gains in the same fit. It makes no sense!\n")
                self.amp_gainsNu = self.amp_gains['nuG']
            elif key == 'tG':
                if self.isMixed:
                    self._printError(
                        "You cannot define split gains and mixed gains in the same fit. It makes no sense!\n")
                self.amp_gainsT = self.amp_gains['tG']
            elif not isinstance(key, int):
                self._printError("The keys of 'amp_gains' must be integers or 'nuG/tG'!")
            else:
                self.amp_gainsNuT[key] = self.amp_gains[key]

        for key in self.amp_gainsNu.keys():
            if not isinstance(key, int):
                self._printError("The keys of 'amp_gains[nuG]' must be integers!")
        for key in self.amp_gainsT.keys():
            if not isinstance(key, int):
                self._printError("The keys of 'amp_gains[tG]' must be integers!")

        self.useGains = set(list(self.phase_gainsNuT) + list(self.amp_gainsNuT)
                            + list(self.phase_gainsNu) + list(self.amp_gainsNu)
                            + list(self.phase_gainsT) + list(self.amp_gainsT))

        # Check phase center:

        if not isinstance(self.phase_center, str):
            self._printError("'phase_center' must be a string!")
        else:
            if len(self.phase_center) == 0:
                self.phrefset = False
            else:
                try:
                    self.phrefset = True
                    csys = cs.newcoordsys(direction=True)
                    dirstr = self.phase_center.split()
                    if len(dirstr) == 2:
                        csys.setdirection(refcode="J2000", refval=self.phase_center)
                    else:
                        csys.setdirection(refcode=dirstr[0], refval=" ".join(dirstr[1:]))
                    csys.convertdirection("J2000")
                    self.refpos = np.copy(csys.torecord()['direction0']['crval'])
                except Exception:
                    self._printError("'phase_center' is not a CASA-formatted sky coordinate!")

        # Did the user forget how does this task work? :D
        for param in self.var:
            if not isinstance(param, str):
                self._printError("'var' must be a list of strings!")
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
            paridx = zip([m.start() + 1 for m in re.finditer(r'\[', self.var[i])],
                         [m.start() for m in re.finditer(r'\]', self.var[i])])
            maxpar = max([maxpar] + list(map(float, [self.var[i][ss[0]:ss[1]] for ss in paridx])))
            if component not in self.implemented:
                msg = "Model component '" + str(component) + "' is not known!\n"
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
        paridx = zip([m.start() + 1 for m in re.finditer(r'\[', self.scalefix)],
                     [m.start() for m in re.finditer(r'\]', self.scalefix)])
        maxpar = max([maxpar] + list(map(float, [self.scalefix[ss[0]:ss[1]] for ss in paridx])))

        # Get the model parameters for the gains:
        #  self.maxGainTerm = 0
        #  self.NgainAnts = {}
        for key in self.phase_gainsNuT:
            term = self.phase_gainsNuT[key]
            paridx = zip([m.start() + 1 for m in re.finditer(r'\[', term)],
                         [m.start() for m in re.finditer(r'\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.amp_gainsNuT:
            term = self.amp_gainsNuT[key]
            paridx = zip([m.start() + 1 for m in re.finditer(r'\[', term)],
                         [m.start() for m in re.finditer(r'\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.phase_gainsNu:
            term = self.phase_gainsNu[key]
            paridx = zip([m.start() + 1 for m in re.finditer(r'\[', term)],
                         [m.start() for m in re.finditer(r'\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.amp_gainsNu:
            term = self.amp_gainsNu[key]
            paridx = zip([m.start() + 1 for m in re.finditer(r'\[', term)],
                         [m.start() for m in re.finditer(r'\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.phase_gainsT:
            term = self.phase_gainsT[key]
            paridx = zip([m.start() + 1 for m in re.finditer(r'\[', term)],
                         [m.start() for m in re.finditer(r'\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        for key in self.amp_gainsT:
            term = self.amp_gainsT[key]
            paridx = zip([m.start() + 1 for m in re.finditer(r'\[', term)],
                         [m.start() for m in re.finditer(r'\]', term)])
            maxpar = max([maxpar] + list(map(float, [term[ss[0]:ss[1]] for ss in paridx])))

        # The fixed model should contain CONSTANT variables:
        for fcomp in self.fixedvar:
            for param in fcomp.split(','):
                if not isinstance(param, float):
                    self._printError("Fixed variables must be a list of floats (or strings representing floats)!")
                    return False

        # There should be as many p_inis as parameters!
        if len(self.p_ini) != maxpar + 1:
            self._printError("'p_ini' is of length %i, but there are %i parameters used!" %
                             (len(self.p_ini), maxpar + 1))
            return False

        # Check for consistency in the fixed model:
        if len(self.fixed) > 0:
            self.takeModel = 'model_column' in self.fixed
        else:
            self.takeModel = False

        if self.takeModel:
            self._printInfo("MODEL COLUMN will be taken as fixed model.\n"
                            + "skipping the 'fixedvar' column and all other fixed components")
            self.fixed = ['model_column']
            self.fixedvar = []
        else:
            for i, component in enumerate(self.fixed):
                checkpars = self.fixedvar[i].split(',')
                if component not in self.implemented:
                    msg = "Model component '" + str(component) + "' is not known!\n"
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
            self._printInfo("setting 'outfile' to its default (i.e., 'modelfit.dat')")
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
        """Reads all the inputs parameters and performs sanity checks.

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

        self._printDebug("uvmultifit::checkInputs")
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
                self._printInfo("SPW is a LIST! User BEWARE! Any fit in *spectral mode*\n"
                                + "will fit the model to each MS separately!")
        elif len(self.spw) > 0:
            self.spw = list([self.spw[0]])
        else:
            self._printError("Bad formatted spw!\n")
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
                    self._printError("'scans' should be a list of integers (or a list of lists of integers, "
                                     + "if there are several measurement sets).")
                    return False
            if len(self.scans) != len(self.vis):
                self._printError("List of (lists of) scans does not have the same length "
                                 + "as the list of measurement sets!")
                return False

            for si, sc in enumerate(self.scans):
                if isinstance(sc, str):
                    self.scans[si] = list(map(int, sc.split(',')))
        except Exception:
            self._printError("'scans' should be a list of integers (or a list of lists of integers, "
                             + "if there are several measurement sets).")
            return False

        # Check dimensions of vis, spw, model, etc.:
        if len(self.vis) != len(self.spw) and len(self.spw) > 1:
            self._printError("The length of 'spw' is not equal to the length of 'vis'!")
            return False

        if not isinstance(self.stokes, str):
            self._printError("'stokes' must be a string!")
            return False

        if not isinstance(self.only_flux, bool):
            self._printError("'only_flux' must be a boolean!")

        if self.only_flux:
            if len(self.p_ini) != len(self.model):
                self._printError("If only_flux=True, number of parameters must be equal to "
                                 + "number of model components!")

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
                self._printError("Measurement set %s does not exist!" % (visi))
                return False

        # Can the required column exist?
        if self.column not in ['data', 'corrected_data', 'corrected']:
            self._printError("'column' can only take values 'data' or 'corrected'!")
            return False
        if self.column == 'corrected':
            self.column = 'corrected_data'

        self.pointing = []
        self.sourscans = []
        phasedirs = [{} for v in self.vis]
        self.field_id = []

        # Get number of antennas:
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
                if isinstance(self.field, (int, str)):
                    fitest = int(self.field)
                else:
                    fitest = int(self.field[vi])
                phasedirs[vi][fitest] = ms.range('phase_dir')['phase_dir']['direction'][:, fitest]
                self.field_id[-1].append(fitest)
            except Exception:
                if isinstance(self.field, str):
                    aux = str(self.field)
                    self.field = [aux for v in self.vis]

                for f, field in enumerate(allfields):
                    if self.field[vi] in field:
                        self.field_id[-1].append(f)
                        phasedirs[vi][f] = ms.range('phase_dir')['phase_dir']['direction'][:, f]
                if len(self.field_id[-1]) == 0:
                    self._printError("field %s is not in %s" % (self.field[vi], v))
                    ms.close()
                    return False

            ms.close()
        self.phasedirs = phasedirs

        # Ref. position:
        if not self.phrefset:
            self._printInfo("setting phase center on first scan")
            self.refpos = self.phasedirs[0][min(self.phasedirs[0].keys())]
        else:
            self._printInfo("setting phase center on %s" % self.phase_center)

        # Find out all the scans where this field id is observed:
        for vi, v in enumerate(self.vis):
            success = ms.open(v)
            if not success:
                self._printError("failed to open measurement set '+v+'!")
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
                        # fi = self.field_id[vi].index(myfield)
                        self.pointing[-1].append(phasedirs[vi][myfield])

            self.pointing[-1] = np.array(self.pointing[-1])

        for vi, v in enumerate(self.vis):
            if len(self.scans[vi]) > 0:
                goodscid = [x for x in self.scans[vi] if x in self.sourscans[vi]]
                if len(goodscid) != len(self.scans[vi]):
                    badscid = [x for x in self.scans[vi] if x not in self.sourscans[vi]]
                    msg = 'the following scans do NOT correspond to source %s: ' % (str(self.field))
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

            # freqdic = ms.getspectralwindowinfo()
            spwchans = ms.range(['num_chan'])['num_chan']

            aux = channeler(self.spw[j], width=self.chanwidth, maxchans=spwchans)
            if aux[0]:
                ranges = list(aux[1])
            else:
                self._printError(aux[1] + "Something seems to be wrong with the 'spw' number %i." % (vi))
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
                self._printError("Bad Stokes parameter %s" % self.stokes)
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
                except Exception:
                    self._printError("Cannot convert to '+self.stokes+'!")
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
                except Exception:
                    self._printError("Cannot convert to '+self.stokes+'")
                    return False
            else:
                self._printError("Polarization " + self.stokes + " not understood.")
                return False

            ms.close()

        # Try to compile the equations for the variables:
        #  if data_changed:
        #  try:
        #    del self.mymodel
        #  except Exception:
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

        self._printDebug("uvmultifit::_setEngineWgt")
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

        self._printDebug("uvmultifit::_setWgtEq")
        tempfloat = 0.0

        if not isinstance(self.dish_diameter, dict):
            try:
                self.dish_diameter = float(self.dish_diameter)
            except ValueError:
                self._printError("The dish diameter must be a number! (in meters)")
                return False

        tb.open(os.path.join(self.vis[0], 'ANTENNA'))
        self.antnames = tb.getcol('NAME')

        if self.pbeam:
            self._printInfo("""
You selected to apply primary-beam correction.\n
PLEASE, remember that the beam is being approximated\n
with a Gaussian, so it may not be very accuracte far\n
from the pointing direction.\n""")
            if isinstance(self.dish_diameter, float):
                if self.dish_diameter == 0.0:
                    try:
                        tempfloat = np.copy(tb.getcol('DISH_DIAMETER'))
                    except Exception:
                        self._printInfo("dish diameter column not found in antenna tables!")
                    tempfloat = np.zeros(len(self.antnames))
                else:
                    self._printInfo("an antenna diameter of %.3f m will be applied" % self.dish_diameter)
                    tempfloat = np.array([self.dish_diameter for a in self.antnames])

            elif isinstance(self.dish_diameter, dict):
                tempfloat = np.array([0.0 for a in self.antnames])
                for anam in self.dish_diameter.keys():
                    for anid in range(len(self.antnames)):
                        antids = re.search(anam, self.antnames[anid])
                        if 'start' in dir(antids):
                            tempfloat[anid] = self.dish_diameter[anam]
                self._printInfo("manual antenna-size setting")
                for anid in range(len(self.antnames)):
                    self._printInfo("antenna %s has a diameter of %.2fm" % (self.antnames[anid], tempfloat[anid]))

            else:
                self._printError("BAD dish_diameter! Should be a float or a dict!")

            if np.max(tempfloat) == 0.0:
                self._printError("The antenna diameters are not set in the ms.\n"
                                 + "Please, set it manually or turn off primary-beam correction.\n")
                return False
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
    def readData(self, del_data=False):
        """Reads the data, according to the properties ``vis, column, chanwidth``, etc.

        It then fills in the properties ``averdata, averfreqs, averweights, u, v, w``, etc.

        Each one of these properties is a list with the data (one list item per spectral window/scan).

        A previous successful run of function ``checkInputs()`` is assumed.

        .. note:: Instead of re-reading data from scratch, using the same uvmultifit instance, a
        better approach may be to restart CASA and create a fresh uvmultifit instance with the
        new data, avoiding some memory leakage related to potential hidden references to the
        data in the IPython's *recall* prompt."""

        self._printDebug("uvmultifit::readData")
        tic = time.time()

        # if del_data:  # data_changed:
        #     # self.deleteData()
        #     #    self.clearPointers(0)
        #     #    self.clearPointers(1)
        #     pass

        self.success = False
        # self._printDebug("inside readData")

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

        # maxDist = 0.0

        # Read data for each spectral window:
        self.iscan = {}
        for vi in self.vis:
            self.iscan[vi] = {}

        for si in nsprang:
            self._printInfo("spectral index #%d (%d of %d)" % (si, si+1, len(nsprang)))
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
            for vis in [x for x in self.spwlist if x[3] <= si]:
                msname = self.vis[vis[1]]

                self._printInfo("opening measurement set '" + msname + "'")
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
                            self._printWarning("spw %i has more than one Data Description ID!" % sp)

                        maskspw = np.zeros(np.shape(crosscorr), dtype=np.bool)

                        for ddi in DDs:
                            maskspw = np.logical_or(maskspw, SPW == ddi)

                        maskspw *= crosscorr

                        # For the first ms in the list, read the frequencies of the spw.
                        # All the other mss will be assumed to have the same frequencies:
                        # if True:
                        tb.open(os.path.join(msname, 'SPECTRAL_WINDOW'))
                        origfreqs = tb.getcol('CHAN_FREQ')
                        self.averfreqs[si] = np.array([np.average(origfreqs[r]) for r in rang])
                        nfreq = len(rang)
                        tb.close()
                        self._printInfo("reading scans for spw %i" % sp)

                        # Read all scans for this field id:
                        for sc, scan in enumerate(self.sourscans[vis[1]]):
                            tb.open(msname)

                            masksc = maskspw * (tb.getcol('SCAN_NUMBER') == int(scan))
                            fieldids = list(np.sort(np.unique(tb.getcol('FIELD_ID')[masksc])))

                            tb.close()

                            for fieldid in fieldids:
                                self._printInfo("reading scan #%i (%i of %i), field: %i" %
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
                                        self._printInfo("offset: %.2e RA (tsec) %.2e Dec (asec)" %
                                                        (phshift[0] / 15., phshift[1]))

                                    # Average spectral channels:
                                    _, ndata = np.shape(datamask)
                                    # ntimes = len(times)
                                    # ntav = int(max([1, round(float(ntimes) / self.timewidth)]))
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
                                    # GaussWidth = 2. * (self.uvtaper / 1.17741)**2.
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
                self.avermod[si] = np.require(np.concatenate(modelscanAv, axis=0),
                                              requirements=['C', 'A'])

            self.averweights[si] = np.require(np.concatenate(weightscan, axis=0),
                                              dtype=np.float64, requirements=['C', 'A'])
            self.u[si] = np.require(np.concatenate(uscan, axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.v[si] = np.require(np.concatenate(vscan, axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.w[si] = np.require(np.concatenate(wscan, axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.t[si] = np.require(np.concatenate(tscan, axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.tArr[si] = np.require(np.concatenate(tArray, axis=0),
                                       dtype=np.float64, requirements=['C', 'A'])
            self.tIdx[si] = np.require(np.concatenate(tIndex, axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])
            self.RAshift[si] = np.require(np.concatenate(RAscan, axis=0),
                                          dtype=np.float64, requirements=['C', 'A'])
            self.Decshift[si] = np.require(np.concatenate(Decscan, axis=0),
                                           dtype=np.float64, requirements=['C', 'A'])
            self.Stretch[si] = np.require(np.concatenate(Stretchscan, axis=0),
                                          dtype=np.float64, requirements=['C', 'A'])
            self.ant1[si] = np.require(np.concatenate(ant1scan, axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])
            self.ant2[si] = np.require(np.concatenate(ant2scan, axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])

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

            # try:
            #     del GaussFact
            # except Exception:
            #     pass
            gc.collect()

        # Initial and final times of observations (time reference for proper motions):
        self.t0 = np.min([np.min(ti) for ti in self.t])
        self.t1 = np.max([np.max(ti) for ti in self.t])

        tac = time.time()
        self._printInfo("reading took %.2f seconds" % (tac - tic))

        self.success = True

        NIFs = len(self.averdata)
        self.Nspw = NIFs

        for i in range(self.Nspw):
            self._printInfo("there are %i integrations (%i visibs.) in spw %i" %
                            (len(self.tArr[i]), len(self.t[i]), i))

        #  import pickle as pk
        #  kk = open('DEBUG.dat', 'w')
        #  pk.dump([self.tIdx, self.tArr, self.t, self.ant1, self.ant2, self.u, self.v, self.w, self.averdata], kk)
        #  kk.close()
        #  raw_input('INPUT')

        #  self.clearPointers(0)

        # Set pointers to data, model, etc.:
        # self.initData(del_data=del_data)
        # self._printDebug("leaving readData")

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

        self._printDebug("uvmultifit::computeFixedModel")

        if self.takeModel:
            self._printInfo("fixed model was taken from model column\n"
                            + "Notice that if you are RE-FITTING, you'll need to RELOAD the model column!")
        else:
            self._printInfo("going to compute fixed model (may need quite a bit of time)")
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
    def initData(self):
        """ Initiates the data pointers of the 'modeler' instance.

        The 'modeler' instance stores pointers to the data and metadata, the compiled model (and
        fixedmodel), the parameter values, and all the methods to compute residuals, Chi Squared, etc."""

        self._printDebug("uvmultifit::initData")
        # Reset pointers for the modeler:
        # self.mymodel.deleteData()

        # Set number of spectral windows:
        gooduvm = uvmod.setNspw(int(self.Nspw))

        if gooduvm != self.Nspw:
            self._printError("Error in the C++ extension!\n")
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
                except Exception:
                    self._printError("You already used the model column! Should read it again!")

            ########
            # Array of booleans, to determine if a datum enters the fit:
            self.mymodel.fittable[-1] = self.mymodel.wgt[-1][:, 0] > 0.0
            # self.mymodel.wgtcorr[-1] = np.require(np.copy(np.tile(PBFactor,(len(self.model), 1))),
            #                                       requirements=['C', 'A'])
            self.mymodel.wgtcorr[-1] = np.require(np.copy(PBFactor), requirements=['C', 'A'])

            del PBFactor

            self.mymodel.fittablebool[-1] = np.require(np.copy(self.mymodel.fittable[-1]
                                                               ).astype(np.bool), requirements=['C', 'A'])

            # print("calling setData")
            # print(f"spidx: {spidx}")
            # print(f"self.mymodel.uv[-1][0]: {self.mymodel.uv[-1][0]}")
            # print(f"self.mymodel.uv[-1][1]: {self.mymodel.uv[-1][1]}")
            # print(f"self.mymodel.uv[-1][2]: {self.mymodel.uv[-1][2]}")
            # print(f"self.mymodel.wgt[-1]: {self.mymodel.wgt[-1]}")
            # print(f"self.mymodel.data[-1]: {self.mymodel.data[-1]}")
            # print(f"self.mymodel.output[-1]: {self.mymodel.output[-1]}")
            # print(f"self.mymodel.freqs[-1]: {self.mymodel.freqs[-1]}")
            # print(f"self.mymodel.fittable[-1]: {self.mymodel.fittable[-1]}")
            # print(f"self.mymodel.wgtcorr[-1]: {self.mymodel.wgtcorr[-1]}")
            # print(f"self.mymodel.dt[-1]: {self.mymodel.dt[-1]}")
            # print(f"self.mymodel.dtArr[-1]: {self.mymodel.dtArr[-1]}")
            # print(f"self.mymodel.dtIdx[-1]: {self.mymodel.dtIdx[-1]}")
            # print(f"self.mymodel.offset[-1][0]: {self.mymodel.offset[-1][0]}")
            # print(f"self.mymodel.offset[-1][1]: {self.mymodel.offset[-1][1]}")
            # print(f"self.mymodel.offset[-1][2]: {self.mymodel.offset[-1][2]}")
            # print(f"self.mymodel.ants[-1][0]: {self.mymodel.ants[-1][0]}")
            # print(f"self.mymodel.ants[-1][1]: {self.mymodel.ants[-1][1]}")
            # print(f"self.mymodel.isGain[-1]: {self.mymodel.isGain[-1]}")
            # print(f"self.mymodel.Nants: {self.mymodel.Nants}")
            gooduvm = uvmod.setData(spidx, self.mymodel.uv[-1][0], self.mymodel.uv[-1][1], self.mymodel.uv[-1][2],
                                    self.mymodel.wgt[-1], self.mymodel.data[-1],
                                    self.mymodel.output[-1], self.mymodel.freqs[-1],
                                    self.mymodel.fittable[-1], self.mymodel.wgtcorr[-1], self.mymodel.dt[-1],
                                    self.mymodel.dtArr[-1], self.mymodel.dtIdx[-1],
                                    self.mymodel.offset[-1][0], self.mymodel.offset[-1][1], self.mymodel.offset[-1][2],
                                    self.mymodel.ants[-1][0], self.mymodel.ants[-1][1],
                                    self.mymodel.isGain[-1], self.Nants)

            # if not gooduvm:
            #     self._printError("Error in the C++ extension!\n")
            #     return False
            #
            # try:
            #     for spidx in range(self.Nspw - 1, -1, -1):
            #         del self.avermod[spidx]
            # except Exception:
            #     pass
            # gc.collect()
        # self._printDebug("leaving initData")
        return True

    ############################################
    #
    #  SET MODEL ARRAYS AND FILL MODEL DATA TO BE SENT TO MODELER
    #
    def initModel(self):
        """ Allocates memory for the modeler data, which will be used by the C++ extension.

        Also compiles the model variables. It is a good idea to run this method everytime that
        the user reads new data to be fitted (i.e., after a call to readData())

        It is indeed MANDATORY to run this function if the model to be fitted has changed."""

        self._printDebug("uvmultifit::initModel")
        if not self.mymodel.initiated:
            # self.clearPointers(1)
            # self.mymodel.deleteModel()
            # else:
            self._printInfo("UVMultiFit compiled model(s) do not seem to exist yet.")

        # Set number of threads:
        ncpu = int(self.NCPU)
        gooduvm = uvmod.setNCPU(ncpu)
        if gooduvm != ncpu:
            self._printError("Error in the C++ extension!")
            return False

        # Initial setup modeler:
        self.mymodel.setup(self.model, self.var, self.fixed, self.fixedvar, self.scalefix,
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
        self._printDebug("model initiated")
        self._printDebug("length of model: %s" % (str(len(self.model))))
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
        self._printInfo("going to compile models")
        self.mymodel.compileAllModels()

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
        self._printInfo("going to run setModel")
        gooduvm = uvmod.setModel(self.mymodel.imod, self.mymodel.Hessian, self.mymodel.Gradient,
                                 self.mymodel.varbuffer, self.mymodel.varfixed, self.mymodel.dpar,
                                 self.mymodel.propRA, self.mymodel.propDec, self.refpos,
                                 self.mymodel.parDependence, self.mymodel.GainBuffer, compFixed, self.isMixed)

        # Delta to estimate numerical derivatives:
        self.mymodel.minnum = self.minDp

        if not gooduvm:
            self._printError("Error in the C++ extension!")
            return False

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
    def fit(self, redo_fixed=True, reinit_model=False, save_file=True, reset_flags=False):
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

        **reset_flags** : `bool`
          Default is False. This is used to clear out the status of *bad data* that may have been
          set by spetial routines of **UVMultiFit** (e.g., the Quinn Fringe Fitter). Default is
          to NOT reset flags (this option should be OK most of the time)."""

        self._printDebug("uvmultifit::fit")
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
                self._printInfo("selecting data by Modified Julian Date")
                timeWindow = np.logical_and(self.t[si] >= self.MJDrange[0],
                                            self.t[si] <= self.MJDrange[1])
                self.mymodel.fittablebool[si][:] = np.logical_and(timeWindow, self.mymodel.fittablebool[si])
                del timeWindow

            self.mymodel.fittable[si][:] = self.mymodel.fittablebool[si][:].astype(np.int8)

        #    else:
        #
        #      self.mymodel.fittable[si][:] = 1
        #      self.mymodel.fittablebool[si][:] = 1

        # Check if there is data available:
            unflagged = np.sum(self.mymodel.wgt[si][self.mymodel.fittablebool[si], :] != 0.0, axis=0)
            # ntot = np.sum(self.mymodel.fittablebool[si])
            if self.OneFitPerChannel:
                if np.sum(unflagged == 0.0) > 0:
                    self._printError("not enough data for this time range! Channels: "
                                     + str(list(np.where(unflagged == 0.0))))
                    self.allflagged = True
                    notfit[si] = list(np.where(unflagged == 0.0)[0])
                    #  return False
            else:
                datatot += np.sum(unflagged)

        #   self.mymodel.output[si][:] = 0.0
        if datatot == 0 and not self.OneFitPerChannel:
            self._printError("not enough data for this time range!\n")
            self.allflagged = True
            return False

        #  for si in range(nspwtot):
        #    self.mymodel.wgtcorr[si][:] = -self.mymodel.KfacWgt
        self._printInfo("now fitting model")

        # Initialize model:
        if reinit_model:
            if self.mymodel.initiated:
                goodinit = self.initModel()
            if self.mymodel.failed or goodinit is False:
                self._printError("bad model (re)initialization! %s" % self.mymodel.resultstring)
                return False

        for i in range(len(self.p_ini) + 1):
            self.mymodel.varbuffer[i][:] = 0.0
            self.mymodel.varfixed[i][:] = 0.0

        # FIT!!!
        ndata = 0.0

        ##################
        # CASE OF SPECTRAL-MODE FIT:
        self._printDebug("one fit per channel: %s" % (str(self.OneFitPerChannel)))
        if self.OneFitPerChannel:
            fitparams = [[] for j in range(nspwtot)]
            fiterrors = [[] for j in range(nspwtot)]
            covariance = [[] for j in range(nspwtot)]
            ChiSq = [[] for j in range(nspwtot)]
            Nvis = [[] for j in range(nspwtot)]

            # Fit channel-wise for each spw:
            for si in range(nspwtot):
                rang = np.shape(self.mymodel.wgt[si])[1]

                for nuidx in range(rang):
                    self._printInfo("fitting channel " + str(nuidx+1) + " of " + str(rang) + " in spw " + str(si))

                    self.mymodel.currspw = si
                    self.mymodel.currchan = nuidx

                    # Compute fixed model:
                    if redo_fixed and len(self.mymodel.fixed) > 0 and not self.takeModel:
                        self.computeFixedModel()

                    # Fit with simplex (if asked for):

                    if nuidx not in notfit[si]:

                        if self.method == 'simplex':
                            fitsimp = _mod_simplex(self.mymodel.ChiSquare, self.p_ini,
                                                   args=(self.bounds, self.p_ini),
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
                                _ = self.mymodel.residuals(fit[0], mode=-3)
                            elif self.write_model == 2:
                                _ = self.mymodel.residuals(fit[0], mode=-4)
                            elif self.write_model == 3:
                                _ = self.mymodel.residuals(fit[0], mode=-5)

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
            self._printInfo("fitting to all frequencies at once")

            # This will tell the modeller to solve in continuum mode:
            self.mymodel.currspw = -1
            self.mymodel.currchan = -1

            # Compute fixed model:
            if redo_fixed and len(self.mymodel.fixed) > 0 and not self.takeModel:
                self._printDebug("Generating fixed model. May take some time")
                self.computeFixedModel()
                self._printDebug("Done!")

            # Pre-fit with simplex:
            if self.method == 'simplex':
                fitsimp = _mod_simplex(self.mymodel.ChiSquare, self.p_ini,
                                       args=(self.bounds, self.p_ini), relxtol=self.SMPtune[0],
                                       relstep=self.SMPtune[1], maxiter=self.SMPtune[2] * len(self.p_ini))
                fit = [fitsimp[0], np.zeros((len(self.p_ini), len(self.p_ini))), fitsimp[1]]
            else:
                # Bound least-squares fitting:
                fit = self.mymodel.LMMin(self.p_ini)
                if not fit:
                    return False

            # self._printInfo("fit results")
            # for x in fit:
            #     print(x)

            # Estimate the parameter uncertainties and save the model in the output array:
            if self.savemodel:
                if self.write_model == 1:
                    _ = self.mymodel.residuals(fit[0], mode=-3)
                elif self.write_model == 2:
                    _ = self.mymodel.residuals(fit[0], mode=-4)
                elif self.write_model == 3:
                    _ = self.mymodel.residuals(fit[0], mode=-5)

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

        self._printInfo("the reduced Chi Squared will be set to 1 by re-scaling the visibility weights.")
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

        self._printInfo("done fitting, now saving to file '%s'" % (self.outfile))

        outf = open(self.outfile, 'w')
        outf.write("# MODELFITTING RESULTS FOR MS: " + ','.join(self.vis) + "\n")
        outf.write("# NOTICE THAT THE REDUCED CHI SQUARE SHOWN HERE IS THE VALUE\n")
        outf.write("# *BEFORE* THE RE-SCALING OF THE VISIBILITY WEIGHTS.\n")
        outf.write("# HOWEVER, THE PARAMETER UNCERTAINTIES ARE *ALREADY* GIVEN\n")
        outf.write("# FOR A REDUCED CHI SQUARE OF 1.\n")
        if isinstance(Nvis, list):
            DOG = int(np.average(Nvis))
        else:
            DOG = Nvis
        outf.write("# AVG. NUMBER OF DEGREES OF FREEDOM: %i\n" % DOG)
        if self.pbeam:
            if isinstance(self.dish_diameter, float):
                outf.write("# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. USING A DISH DIAMETER OF: %.3f METERS\n" %
                           self.dish_diameter)
            else:
                outf.write("# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. USING: \n")
                for iant in range(len(self.antnames)):
                    outf.write("#    FOR ANTENNA %s A DIAMETER OF %.2fm \n" %
                               (self.antnames[iant], self.userDiameters[iant]))

        else:
            outf.write("# PRIMARY-BEAM CORRECTION HAS NOT BEEN APPLIED.\n")

        outf.write("###########################################\n")
        outf.write("#\n# MODEL CONSISTS OF:\n")
        for m, mod in enumerate(self.model):
            outf.write("# '" + mod + "' with variables: " + self.var[m] + "\n")
        if len(self.fixed) > 0:
            outf.write("#\n#\n# FIXED MODEL CONSISTS OF:\n")
            if self.takeModel:
                outf.write("# THE MODEL COLUMN FOUND IN THE DATA\n")
            else:
                for m, mod in enumerate(self.fixed):
                    var2print = list(map(float, self.fixedvar[m].split(',')))
                    outf.write(
                        "# '" + mod + "' with variables: " + ' '.join(['%.3e'] * len(var2print)) % tuple(var2print))
                outf.write("#\n#  - AND SCALING FACTOR: %s\n" % self.scalefix)

        outf.write("#\n# INITIAL PARAMETER VALUES:\n")
        for p0i, p0 in enumerate(self.p_ini):
            if self.bounds is not None:
                outf.write("#  p[%i] = %.5e with bounds: %s\n" % (p0i, p0, str(self.bounds[p0i])))
            else:
                outf.write("#  p[%i] = %.5e with no bounds\n" % (p0i, p0))

        outf.write("#\n##########################################\n")
        parshead = []

        for pp in range(len(self.p_ini)):
            parshead.append(pp)
            parshead.append(pp)

        headstr = (
            "# Frequency (Hz)     " + "  p[%i]       err[%i]   " * len(self.p_ini) + "   red. ChiSq\n") % tuple(parshead)
        outf.write(headstr)

        formatting = "%.12e   " + "%.4e " * (2 * npars) + "   %.4e \n"
        if not self.OneFitPerChannel:
            toprint = tuple([np.average(self.averfreqs)] + prtpars + [ChiSq])
            outf.write(formatting % toprint)

        else:
            for spwvi in self.spwlist:
                for r, rr in enumerate(spwvi[2]):
                    k = spwvi[3] + r
                    rang = rr[1]
                    for nu, freq in enumerate(self.averfreqs[k]):
                        toprint = tuple([freq] + prtpars[k][nu] + [ChiSq[k][nu]])
                        outf.write(formatting % toprint)
        outf.close()
        tac = time.time()

        self._printInfo("fit took %.2f seconds" % (tac - tic))
        #  uvmod.unsetWork()

        return True


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
        # totflux = np.sum(cleans[:, 2])

        # Put shifts in arcseconds w.r.t. image center:
        cleans[:, 0] = (cleans[:, 0] - refpix[0]) * deltaxy[0] * 180. / np.pi * 3600.
        cleans[:, 1] = (cleans[:, 1] - refpix[1]) * deltaxy[1] * 180. / np.pi * 3600.

        # Set model as a set of deltas:
        model = ['delta' for cl in cleans]
        var = ['%.3e, %.3e, %.3e' % tuple(cl) for cl in cleans]

        # Send results:
        ia.close()
        return [True, model, var, [cleanlocs[0], cleanlocs[1]]]
    except Exception:
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
                 return_all=False):
    """
    Minimization of scalar function of one or more variables using the
    Nelder-Mead algorithm.

    Options for the Nelder-Mead algorithm are:
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
                        # print('simj', sim[j])
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
    # warnflag = 0

    return [x, fval]
