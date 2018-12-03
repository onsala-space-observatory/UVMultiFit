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
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#

import os, sys
global goodclib, ms

__version__ = "3.0.0"
date = 'MAR 2018'

################
# Import all necessary modules.

mypath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(mypath)
########
# Execute twice, to avoid the silly (and harmless)
# error regarding the different API versions of
# numpy between the local system and CASA:
try:
    import _uvmultimodel as uvmod
    goodclib = True
    print '\nC++ shared library loaded successfully\n'
except:
    goodclib = False
    print '\n There has been an error related to the numpy'
    print ' API version used by CASA. This is related to uvmultifit'
    print ' (which uses the API version of your system) and should'
    print ' be *harmless*.\n'

if not goodclib:
    try:
        import _uvmultimodel as uvmod
        goodclib = True
        print '\nC++ shared library loaded successfully\n'
    except:
        goodclib = False
#############

import sys

import time
import numpy as np
import scipy as sp
import scipy.optimize as spopt
import scipy.special as spec
import os
import gc
from taskinit import gentools
import re
import shutil
from clearcal_cli import clearcal_cli as clearcal
ms = gentools(['ms'])[0]
ia = gentools(['ia'])[0]
cs = gentools(['cs'])[0]

greetings = '\n ##########################################################################\n'
greetings += ' # UVMULTIFIT --  ' + date + '. EUROPEAN ALMA REGIONAL CENTER (NORDIC NODE).  #\n'
greetings += ' #       Please, add the UVMULTIFIT reference to your publications:       #\n'
greetings += ' #                                                                        #\n'
greetings += ' #      Marti-Vidal, Vlemmings, Muller, & Casey 2014, A&A, 563, A136      #\n'
greetings += ' #                                                                        #\n'
greetings += ' ##########################################################################\n\n'


class uvmultifit(object):
    """ .

      ALL DEFAULTS:


     uvmultifit(self,vis='', spw='0', column = 'data', field = 0, scans = [],
             uniform=False, chanwidth = 1, timewidth = 1, stokes = 'I',
             write='', MJDrange=[-1.0,-1.0], model=['delta'], 
             var=['p[0],p[1],p[2]'], p_ini=[0.0,0.0,1.0], phase_center = '',
             fixed=[], fixedvar=[], scalefix='1.0', outfile = 'modelfit.dat',
             NCPU=4, pbeam=False, ldfac = 1.22, dish_diameter=0.0, gridpix=0, 
             OneFitPerChannel=True, bounds=None, cov_return=False,
             finetune=False, uvtaper=0.0, method='levenberg', wgt_power=1.0, 
             LMtune=[1.e-3,10.,1.e-5,200,1.e-3], SMPtune=[1.e-4,1.e-1,200], 
             only_flux=False, proper_motion = 0.0, HankelOrder = 80, 
             amp_gains = {}, phase_gains = {})


      Fits a multi-component model to the visibilities in a measurement set. 
      The fit is performed either to each spectral channel (or chunck of 
      spectral channels) or to all channels together. In case of fit to each 
      frequency bin, no correlation between channels (e.g., no spectral 
      index) is assumed, so the fit to each channel is independent of the 
      fits to the other channels. In case of fit to all channels, the 
      spectral index can be set as either a fixed value or a fitting parameter.

      The uncertainties are estimated from the Jacobian matrix and 
      scaled so that the reduced Chi squared equals unity. Null 
      (or too small) uncertaintiesm may indicate an unsuccessful fit.

      The best-fit results are stored in an external file and also returned
      to the user as a dictionary. The elements of the dictionary are:
      'Frequency','Parameters','Uncertainties', and 'Reduced Chi squared'.

      The array named 'Reduced Chi squared' is the reduced post-fit Chi 
      squared, computed with the uncertainties taken from the data. That means
      it is computed *before* scaling the data uncertainties (and hence the 
      best-fit parameter uncertainties) to scale the reduced Chi squared to 
      its expected value (i.e., unity). Quite large values are likely 
      indicative of too low data weights. This should be harmless to the fit, 
      since uvmultifit scales the weight a-posteriori to make the reduced Chi 
      square equal to its expected value. Too high values of Chi squared could 
      also be indicative of a bad fit. In such a case, the a-posteriori 
      parameter uncertainties would also be large. The user should check these 
      uncertainties as an assessment of the quality of the fit.

      The user can save a "uvmultifit instance" with call to 'uvmultifit' like
      this one:

         myfit = uvmultifit(vis=... )

      Then, the user can look at the contents of the 'myfit.result' dictionary 
      to recover the fitting results. For instance, if the user wants to plot 
      the first fitted parameter (i.e., p[0], see below) against frequency for 
      the third spectral window, the command would be:

         plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3][:,0])

      To plot the second fitted parameter (i.e., p[1]), just execute:

         plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3][:,1])
     
      and so on. If there is only 1 fitting parameter:

         plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3])

      NOTE: The fitted parameters are ordered by spectral window *in the order 
        given by the "spw" parameter*. Hence, if spw='2,3', then the first 
        element of the list of parameters (i.e., myfit.result['Parameters'][0]) 
        will be the parameters fitted for the spw number 2.

      NOTE 2: For fits to all channels together (see below), 
        myfit.result['Parameters'] will only have one entry (i.e., the 
        parameter values fitted for all channels at once).

      NOTE 3: 'uvmultifit' will return an instance of a class with many 
        properties. One of these properties is a set of arrays with all the 
        data (averaged in time and/or frequency, according to the "chanwidth" 
        and "timewidth" parameters given by the user; see below). Hence, a 
        uvmultifit instance can occupy quite a bit of physical memory. It may 
        be a good idea to remove the instance from memory (i.e., using 
        "del myfit" and "gc.collect()") or restart CASA once the user has got 
        the results of the fit (this is the preferred approach, since it's 
        difficult to remove all references to an object in an
        Interactive-Python session, due to the "interactivity" of the prompt).


      The input parameters of uvmultifit() are:

      vis = Name of the measurement set. It can also be a list of MS names.

      spw = String or list of strings. These are the spectral window(s) and 
            channel selection(s) to be fitted. For each spw, the user can 
            select many channel ranges in the usual CASA way. If one string is 
            given, all ms are supposed to have the same spectral configuration. 
            If a list of strings is given, one per ms, this restriction 
            doesn't apply.
            
            Example 1: vis = ['Cband.ms',  'Lband.ms']
                       spw = ['0',  '0:10~100']
          
            In this example, the spw number 0 in the 'Cband.ms' measurement set,
            together with the channels 10 to 100 of spw 0 in 'Lband.ms' will be
            fitted. 'Cband.ms' and 'Lband.ms' don't have to have the same 
            frequency configuration.

            Example 2: vis = ['Cband1.ms',  'Cband2.ms']
                       spw = '0:1~50;2'

            In this example, the channels 1 to 50 of spw 0 and the whole spw 2
            will be fitted. 'Cband1.ms' and 'Cband2.ms' must both have the same 
            frequency configuration (i.e., same spw frequencies and channels)


      column = The data column. It can be either 'data' or 'corrected'

      field = The id number (or name) of the target source in the ms.

      pbeam = Boolean. If false, the primary-beam correction is not applied. 
              This is *important* for fitting mosaic data.

      dish_diameter = In case that the antenna diameters cannot be read from 
              the datasets, the user must provide it (in meters).

      ldfac = Float. Proportionality constant between lambda/Diameter and 
               the FWHM of the primary beam (assumed to be Gaussian!). I.e.:
           
               FWHM = ldfac*lambda/Diameter

               Normally, ldfac = 1.22 should be fine, although 1.00 works
               better with data from simobserve.

      gridpix = Integer. Default 0. If >0, the visibilities are gridded using
             the CASA task *clean*. The primary-beam corrections are then 
             applied by clean, so the "pbeam" and "dish_diameter" parameters of 
             uvmultifit are overriden by clean. The number of pixels in Fourier 
             space is set by the value of gridpix. THIS OPTION MAY BE USEFUL 
             FOR VERY LARGE SOURCES, WITH SIZES OF THE ORDER OF (OR LARGER 
             THAN) THE PRIMARY BEAM.

      scans = List of integers; default []. The id numbers of the scans to load.
              Default means to load all the scans of the selected source. If 
              multiple measurement sets are selected, there should be a list of 
              lists (i.e., one list per measurement set). 
              For instance, if vis=["ms1.ms","ms2.ms"], then scans=[[1,2,3],[]] 
              would select the scans 1, 2, and 3 from the "ms1.ms" dataset and 
              all the scans from the "ms2.ms" dataset.

      chanwidth = number of spectral channels to average in each chunk of 
                  data to be fitted.

      timewidth = number of time channels (i.e., integration times) to 
                  average in each chunk of data. The averaging is always 
                  cut at the scan boundaries (so a very large timewidth 
                  means to average down to one visibility per scan).

      MJDrange = list of two floats. These are the initial and final Modified 
                 Julian Dates of the data to be used in the fitting. All the 
                 data are loaded a-priori, and the MJDrange condition is 
                 applied afterwards. This way, variability studies can be 
                 performed efficiently. Default (i.e., <= 0.0) means not to 
                 select data based on JD time range.

      stokes = polarization product. Can be any of 'PI','I','Q','U','V'. It 
               also accepts individual correlation products: 'XX', 'YY', 'XY',
               'YX', 'RR', 'LL','LR', or 'LR'. 
               Default is 'I'. If 'PI' is given, the program will compute 'I' 
               whenever possible and use 'XX', 'YY', 'RR' or 'LL' otherwise. 
               This way, the amount of data used in a polarization-independent 
               fits is maximized.

      model = list of model components to fit. Each component is given as a 
              string. 
          Possible models are: 'delta', 'disc', 'Gaussian', 'ring', 'sphere'
                               'bubble', 'expo', 'power-2', 'power-3'
                               'GaussianRing'

      var = list of strings. These are the variables for each model. 
             The variables can be set to *any* algebraic expression 
             involving the fitting parameters (being the ith parameter 
             represented by 'p[i]') and the observing frequency in Hz 
             (represented by 'nu'). Any numpy function can also be called,
             using the prefix 'np' (e.g., 'p[0]*np.power(nu,p[1])'). 
             See some examples below.

      p_ini = list of the initial values of the parameters. This is expected
              to be a list of floats.

      phase_center = string. The sky position where all components are
                     referenced to. If an empty string is given, the phase 
                     center will be that of the first field id whose is being
                     read (i.e., if a mosaic is being read, the first pointing
                     will be the phase center).
                     If the string is not empty, the program expects to find
                     a coordinate in CASA format (i.e., 'J2000 RA Dec'), where 
                     RA is in format '00h00m00.0s' and Dec is in format 
                     '00d00m00.0s'.

      fixed = Like 'model', but defines a model with completely fixed
              variables (i.e., defined only by numbers; not parameters).
              This model will be computed only once, hence making the 
              code execution faster.
 
      fixedvar = Like 'var', but refers to the 'fixed' model. Hence,
                  it is expected to be either a list of numbers or a list 
                  of strings representing numbers.
              
      scalefix = string representing a function that sets a scaling 
                 factor for the fixed model. It can indeed be a function of 
                 the fitting parameters (e.g., scalefix='p[0]' will multiply 
                 the overall flux density of the fixed model by p[0]) and 
                 the observing frequency 'nu'.
                  

      OneFitPerChannel = boolean. If True, fits are performed to the different
                  frequency channels, one by one. If False, one common fit is 
                  performed to all data. In this case, the user may also want 
                  to fit the spectral indices of the components.


      outfile = Name of the file to store results.

      bounds = list of boundaries (i.e., minimum and maximum allowed values) 
               for the fitting parameters. 'None' means that no bound is 
               defined. If the list is empty, no bounds are assumed for any 
               parameter.
 
      cov_return = If True, the covariance matrix for each fit is added to 
                the returning dictionary (with key name 'covariance')

      uvtaper = float. Default is 0.0. If not 0.0, the weights of the 
                visibilities are multiplied by a Gaussian in Fourier space,
                whose HWHM is the value of \'uvtaper\', *in meters*.

      uniform = Boolean. Default is False. If True, the weights of all data 
                are made equal. Notice that a uvtaper can still be applied 
                (i.e., \'uvtaper\' parameter).

      finetune = Boolean. Default is False. If set to true, the fit is not
                 performed, but only a 'uvmultifit' instance is created. The
                 user can then run the different methods of the instance 
                 by him/herself (see the help text for each method). This can
                 be useful, for instanse, if the user wants to try many 
                 different models (and/or subtract many different fixed 
                 models), without having to reload the data every time after 
                 each step. 

      wgt_power = Float. Default is 1. Power index of the visibility weights 
                 in the computation of the Chi square. wgt_power = 1 would be 
                 the "statistically justified" value, although other values 
                 may be tested if the user suspects that the visibility 
                 uncertainties are not well estimated in his/her dataset.

      method =  String. Default is 'simplex'. Possible values are 'simplex' and 
                'levenberg'. Sometimes, the least-squares minimization may not 
                converge well with the Levenberg-Marquardt method (if the 
                model is complicated and/or the uv coverage is quite sparse). 
                Levenberg-Marquardt also requires a lot of memory. In all 
                these cases, a Chi square minimization using the SIMPLEX 
                algorithm may work beter. However, SIMPLEX may require more 
                function evaluations.

           ##############
            WARNING: The SIMPLEX method does not currently provide parameter
                       uncertainties. It just returns the reduced Chi squared!
           ##############


      write = String. Default is ''. If it is set to 'model', the best-fit 
                   model is saved in the 'model' column of the measurement set.
                   If it is set to 'residuals', the fit residuals are saved 
                   in the 'corrected' column of the measurement set.

                   Currently, this only works if stokes is set to either 
                   'PI', 'I', or an individual correlation product (like 
                   'XX' or 'XY') AND if both timewidth and chanwidth are 
                   set to 1.

      NCPU = Integer. Default is 4. Number of threads allowed to run.

      SMPtune = list of 2 floats and one integer. Used to fine-tune the Simplex 
                algorithm. Change only if you really know what you are doing.

              -SMPtune[0] -> Default is 1.e-4. Maximum allowed error in the 
                  parameters, from the search of the Chi Square minimum. 
              
              -SMPtune[1] -> Default is 1.e-1. Relative size of the first step 
                 in the parameter search.

              -SMPtune[2] -> Maximum number of iterations allowed per fitting 
                 parameter. Default: 200

      LMtune = list of 3 floats and one integer. Used to fine-tune the 
               Levenberg-Marquardt algorithm. Change only if you really know 
               what you are doing.

             -LMtune[0] -> The Lambda factor to weight up the Newton component 
               (i.e., H + Lambda*H_diag = Res), where H is the Hessian and 
               Res the vector of residuals. Default: 1.e-3

             -LMtune[1] -> The factor to multiply (divide) Lambda if the 
                iterator worsened (improved) the Chi Squared. Default: 10.

             -LMtune[2] -> The maximum relative error allowed for the ChiSq. 
                           Default: 1.e-5

             -LMtune[3] -> Maximum number of iterations allowed per fitting 
                parameter. Default: 200

             -LMtune[4] -> (optional) Maximum relative error allowed for 
                           parameters. Default: 1.e-3. 
                           If not provided, LMtune[4] = LMtune[2]


      only_flux: Boolean. Default False. If True, the program assumes that 
                 only the flux densities of all components are going to be 
                 fitted. Furthermore, the fitting parameters shall be just 
                 equal to the flux densities being fitted, and should be given 
                 in the same order as the list of model components. In these 
                 cases, using only_flux = True can accelerate the fit quite a 
                 bit, especially if there are many components to be fitted.

      proper_motion: list of 2-element float lists: Each element (i.e., 2-float 
                     list) is the proper motion, in RA and Dec, of each model 
                     component. The units are arc-seconds per year. Proper 
                     motions cannot be fitted yet, but may be fittable in 
                     future versions of uvmultifit.
                     Default is a float = 0.0, meaning all proper motions are 
                     null. If the proper motions are not null, the position 
                     used as reference for the fit will correspond to that of 
                     the first scan of the observed field.

      HankelOrder: Only used for models without an analytic expression
                   (GaussianRing). In these cases, UVMULTIFIT performs the 
                   Hankel transform by using the series expansion of the 
                   Bessel function J0. The order of this expansion is set
                   by this keyword. The larger the distance in Fourier space
                   (i.e., the larger the source, in beam units), the higher
                   HankelOrder should be, to keep the numerical precision.
                   Default is 80. For cases of sources with sizes similar 
                   to the synthesized beam, HankelOrder=80 should suffice.
                   The use of too high values may cause Overflow errors!
                   


      amp_gains: Dictionary (default empty). Controls whether to solve for
                 antenna amplitude gains and source-model parameters 
                 *simultaneously*. The keys of the dictionary are the 
                 indices of the antennas to self-calibrate. The values of the 
                 dictionary are strings, whose elements are comma-separated
                 functions of the fitting parameters. Each element of this
                 list of functions will be the coefficient of a time 
                 polynomial used to model the gains. For instance, if we want
                 to self-calibrate antennas number 1 and 2, using a linear 
                 time dependence of the amplitude gains, then:

          amp_gains = {1:'p[0],p[1]',2:'p[2],p[3]'}

          This will solve for the gain of antenna 1 as: GAIN = p[0] + p[1]*t
          and the gain for antenna 2 as: GAIN = p[2] + p[3]*t

          where "t" is the time elapsed since the first observed visibility,
          measured in seconds. The fitting parameters here are treated in 
          exactly the same way as the source parameters, so bounds can also 
          be set.

      phase_gains: Same as for amp_gains

      EXAMPLE: If the want to solve for a constant amplitude scale of 
               antenna #2, just set:

          amp_gains = {2: 'p[0]'}

        This way, p[0] will be the amplitude gain of antenna #2. The source
        parameters can be defined using other parameters, p[i], which will be
        solved simultaneously with p[0].

      NOTE: There is a new function, "pieceWise(a0, a1, t0, t1)" that can be 
       used when defining the "amp_gains" and "phase_gains" functions. It 
       returns a linear interpolation (between times t0 and t1), using a0 as 
       the value of the function in t0 and a1 as the value in t1. The value 
       for times earlier than t0 are extrapolated as a0, and the times after 
       t1 are extrapolated as a1. This function is useful to approximate the 
       antenna gains as piece-wise linear functions. For instance, a 
       piece-wise function with two nodes would be:

       pieceWise(a0,a1,t0,t1) + piecewise(0,a2-a1,t1,t2)

       The total function would cross the points: (a0,t0), (a1,t1), (a2,t2)

     




      -----------------------------------------

      MODEL DETAILS: Each model depends on a list of variables that can 
      be written as *any* algebraic combination of fitting parameters. Some 
      explanatory examples are given in the MODEL EXAMPLES section below. 
      Here is the list of currently supported models, and their variables:

      -'delta'    -> Variables: RA, Dec, Flux
      -'disc'     -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle
      -'ring'     -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle
      -'Gaussian' -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle
      -'sphere'   -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle
      -'bubble'   -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle
      -'expo'     -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle
      -'power-2'  -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle
      -'power-3'  -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle
      -'GaussianRing'  -> Variables: RA, Dec, Flux, Major, Ratio, PositionAngle, Sigma


        'sphere' stands for an optically-thin uniform filled sphere.
        'bubble' stands for a uniform spherical surface.
        'expo' stands for an exponential radial flux decay.
        'power-2' stands for a decay as 1/(r^2 + r0^2) (notice that in 
                this case, the flux is the integral from r=0 to r=r0)
        'power-3' stands for a decay as 1/(1 + (2^2/3 - 1)(r/r0)^2)^(3/2).

        'GaussianRing' stands for a radial profile following the expression:

           F*exp(-(Q - R)**2./(2*Sigma**2.))

         NOTE: The GaussianRing is built by a series expansion of the Bessel
         function in the computation of the Hankel transform. The order of the
         expansion (default is 10) can be set in the new HankelOrder keyword.
         The bigger the source (in beam units) the higher HankelOrder should
         be, to keep the numerical precision. For cases of source sizes similar
         to the beam, an order of 5 should suffice (FT precision within 1%). 

      * RA and Dec are the shifts w.r.t. the phase center (in arcsec)
      * Flux is the total flux density of the component (in Jy)
      * Major is the diameter along the major axis
      * Ratio is the size ratio between the reference axis and the other axes 
              (i.e., it is set to 1.0 for circularly-symmetric sources)
      * PositionAngle is the angle of the reference axis, from North to East 
        (in deg.)

      * Sigma (whenever it is used) is an auxiliary variable for models that 
        need more than one 'size-like' parameter to be defined (e.g., the 
        GaussianRing, where we have the size of the ring and its width).   

      ----------------
      MODEL EXAMPLES:
      ----------------

      0.- One delta, with position and flux set as free parameters:

         model = ['delta']
         var = ['p[0],p[1],p[2]']

         In this case, the RA shift w.r.t. the image center is p[0] (in 
         arcsec), the Dec shift is p[1] (also in arcsec), and the flux 
         density is p[2] (in Jy).

         If we know that the RA and Dec shifts of the delta w.r.t. the
         phase center cannot be larger than 1 arcsec, we could set:

         bounds = [[-1,1], [-1,1], None]

         If we also want to force the flux density to be positive:

         bounds = [[-1,1], [-1,1], [0,None]]


      1.- Two deltas, being the position of the second one fixed w.r.t. 
          the position of the first one. Let's say that the second delta
          is shifted at 0.5 and 0.6 arcsec, in RA and Dec (respectively),
          from the first delta, and we want the position of the first
          delta to be free in our fit. Then:
     
         model = ['delta','delta']
         var = ['p[0],p[1],p[2]','p[0]+0.5, p[1]+0.6, p[3]']
      
         In this case, p[0] and p[1] are the RA and Dec of the first 
         delta; p[2] is the flux density of the first delta; and p[3] 
         is the flux density of the second delta.
         (notice that the RA and Dec position of the second delta is 
         equal to that of the first delta plus a fixed shift).



      2.- A ring plus a delta at its center. The absolute position of 
          the compound source is free, and the ring is circular.

         model = ['ring','delta']
         var = ['p[0],p[1],p[2],p[3],1.0,0.0','p[0],p[1],p[4]']

         In this case, p[0] and p[1] are the RA and Dec of both 
         components (i.e., the ring and the delta); p[2] is the total 
         flux density of the ring and p[3] is its diameter; p[4] is the 
         flux density of the delta. Notice that the axes Ratio of the
         ring is set constant (and unity) and the position angle is
         also set constant (although it's meaningless for Ratio=1).

         For extended models, it is a good idea to bound the size to
         have positive values. Hence, in this case:

         bounds = [None,None,None,[0,None]]

         In case we also want to bound the fitted fluxes to be 
         positive, we would have:

         bounds = [None,None,[0,None],[0,None],[0,None]]


      3.- Like example 1, but fixing also the flux-density of the 
          second delta to be 2.5 times that of the first delta:

         model = ['delta','delta']
         var = ['p[0],p[1],p[2]','p[0]+0.5,p[3]+0.6,p[2]*2.5']


      4.- A circularly-symmetric disc with a hole (i.e., with its inner 
          half subtracted):

       model = ['disc','disc']
       var = 
        ['p[0],p[1],4/3*p[2],p[3],1,0','p[0],p[1],-p[2]/3,p[3]/2,1,0']

       In this case, p[0] and p[1] are the RA and Dec shifts, 
       respectively; p[2] is the flux density of the disc with the hole
       (i.e., with its inner half subtracted), and p[3] is the disc 
       diameter. Notice that the hole in the disc has been created by 
       adding a *negative* disc of size equals to 1/2 of the size 
       of the larger disc, and flux density equals to -1/4 of that
       of the larger disc. The overall effect of both discs is that of 
       one single disc with a hole (i.e., with no emission at radii 
       < 1/2 of the outer radius). 


       5.- A delta component with a spectral index, fitted to the 
           whole dataset (i.e., setting OneFitPerChannel=False):

             model = ['delta']
             var = ['p[0],p[1],p[2]*(nu/1.e9)**p[3]']
   
          In this case, the fitted spectral index will be p[3] and p[2]
          will be the flux density at 1GHz. Notice that p_ini (the list 
          of initial values for the parameters) must also have the a 
          priori value of the spectral index. 
          For this example, p_ini could be:

             p_ini = [0.0,0.0,2.3,-0.7]

          for a source close to the center, with an approximate flux 
          of 2.3 Jy at 1GHz and an a priori spectral index of -0.7.

          If the spectral index is well known, it can of course be fixed 
          in the fit. In such a case:

             model = ['delta']
             var = ['p[0],p[1],p[2]*(nu/1.e9)**(-0.7)']
             p_ini = [0.0,0.0,2.3]

          NOTE: Fitting sizes *and* spectral indices at the same time
          may produce crazy results, since both things are quite coupled
          in Fourier space. Hence, some care should be taken when fitting
          to all frequencies together. Notice also that, unless your bandwidth
          is quite wide and/or your SNR is quite high, any fit of the 
          spectral index may not be reliable.


       6.- A filled sphere with its inner half (i.e., the core) removed. 
           The sphere is fixed at the image center:

          model = ['sphere','sphere']
          var = ['0,0,9/8*p[0],p[1],1,0','0,0,-p[0]/8,p[1]/2,1,0']

          In this case, p[0] is the total flux density and p[1] is the outer
          diameter. This example is similar to that of the disc with a hole, 
          but in this case the model is a sphere where we have removed all 
          the emission at radii < 0.5 times the outer radius.


       7.- A disc with a hole of variable size. In this case, the smaller 
          disc with negative flux density, that we use to generate the hole, 
          must *always* have the same surface intensity (in absolute value) 
          that the larger positive disc. Hence, and if we fix the disc 
          position at the image center (for simplicity), we have:

          model = ['disc','disc']
          var = ['0,0,p[0]*p[1]**2./(p[1]**2.-1),p[2],1,0',
                        '0,0,-p[0]/(p[1]**2.-1),p[2]/p[1],1,0']

          In this case, the flux density of the disc with the hole is p[0],
          the outer dimaeter is p[2], and p[1] is the ratio of the outer 
          size to the inner size. 


 """

    global sys, goodlib, ms, uvmod, goodclib, time, np, sp, spopt, spec, os, gc, gentools, re, clearcal, ms, cs

    ############################################
    #
    #  FREE MEMORY
    #

    def deleteData(self):

        for i in range(len(self.averdata) - 1, -1, -1):
            del self.averdata[i]
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

        for mdi in self.iscancoords:
            Npar = len(mdi)
            for mdp in range(Npar - 1, -1, -1):
                del mdi[mdp]
            del mdi
        del self.iscancoords

        del self.averdata, self.averweights, self.averfreqs, self.v, self.u, self.w, self.t
        del self.ant1, self.ant2, self.RAshift, self.Decshift, self.Stretch

    def __del__(self):

        self.deleteData()
        del self.mymodel

############################################
#
#  CREATE INSTANCE
#

    def __init__(self,
                 vis='',
                 spw='0',
                 column='data',
                 field=0,
                 scans=[],
                 uniform=False,
                 chanwidth=1,
                 timewidth=1,
                 stokes='I',
                 write='',
                 MJDrange=[-1.0, -1.0],
                 ldfac=1.22,
                 model=['delta'],
                 var=['p[0],p[1],p[2]'],
                 p_ini=[0.0, 0.0, 1.0],
                 phase_center='',
                 fixed=[],
                 fixedvar=[],
                 scalefix='1.0',
                 outfile='modelfit.dat',
                 NCPU=4,
                 pbeam=False,
                 dish_diameter=0.0,
                 gridpix=0,
                 OneFitPerChannel=True,
                 bounds=None,
                 cov_return=False,
                 finetune=False,
                 uvtaper=0.0,
                 method='levenberg',
                 wgt_power=1.0,
                 LMtune=[1.e-3, 10., 1.e-5, 200, 1.e-3],
                 SMPtune=[1.e-4, 1.e-1, 200],
                 only_flux=False,
                 proper_motion=0.0,
                 HankelOrder=80,
                 phase_gains={},
                 amp_gains={}):
        """ Just the constructor method, for class instantiation."""

        self.first_time = True
        self.implemented = [
            'delta', 'disc', 'ring', 'Gaussian', 'sphere', 'bubble', 'expo',
            'power-2', 'power-3', 'GaussianRing'
        ]
        numvars = [3, 6, 6, 6, 6, 6, 6, 6, 6,
                   7]  # Number of variables for the models
        self.isNumerical = ['GaussianRing'
                            ]  # Models that need a series expansion for J0.
        self.uniform = uniform
        self.vis = vis
        self.column = column
        self.averdata = []
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
        self.timewidth = timewidth
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
        self.gridpix = 0
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
            self.write_model = {'': 0, 'model': 1, 'residuals': 2}[write]
        except:
            self.printError(
                'ERROR: keyword \'write\' should be set to either \'\', \'model\', or \'residuals\' '
            )
            return False

#  try:
#   test = self.scipylm
#  except:
#   self.scipylm = False

# Start instance:
        self.startUp()

############################################
#
#  STARTERS
#

###########################################
# This method will be overriden in the GUI, to avoid execution of CheckInputs() and the fit:

    def startUp(self):

        if not goodclib:
            self.printError(
                "ERROR: C++ library cannot be loaded! Please, contact the Nordic ARC node."
            )
            return False

        self.printInfo(greetings)

        # Check parameters and read the data in:
        goodread = self.checkInputs(data_changed=True)
        if goodread:
            goodread = self.readData(data_changed=False)
            self.initData()
        else:
            self.printError('\n\n ABORTING UVMULTIFIT! CHECK INPUTS!\n')

# Fit if finetune==False:
        if goodread:
            if self.finetune == False:
                goodfit = self.fit()
                if not goodfit:
                    self.printError('\n FAILED FIT!\n')
                elif self.write_model in [1, 2]:
                    if self.timewidth == 1 and self.stokes not in [
                            'Q', 'U', 'V'
                    ]:
                        self.printInfo('\nWriting into measurement set(s)\n')
                        self.writeModel()
                    else:
                        msg = '\nCANNOT write into measurement set!\n'
                        msg += '\nIf you want to fill-in the model (or corrected) column:\n'
                        msg += '    1.- \'timewidth\' and \'chanwidth\' should be set to 1\n'
                        msg += '    2.- \'stokes\' should be set to either \'I\', \'PI\', or a corr. product.'
                        self.printError(msg)

                if goodfit:
                    self.printInfo(
                        '\n\n\nFit done!! And UVMULTIFIT class instantiated successfully!\n'
                    )
            else:
                self.printInfo(
                    '\n\n\nUVMULTIFIT class instantiated successfully!\n')

############################################
#
#  PRINT MESSAGES AND ERRORS
#

###################
# Functions overriden in GUI mode:

    def printError(self, message):
        sys.stdout.write(message)
        sys.stdout.flush()
        raise ValueError(message)

    def printInfo(self, message):
        sys.stdout.write(message)
        sys.stdout.flush()
###################
#################################################

############################################
#
#  WRITE MODEL (OR RESIDUALS)
#

    def writeModel(self):
        """ Writes the predictions of the compiled model(s) into the \'model\' 
      column of the measurement set(s). This only executes if \'stokes\' is 
      set to 'PI','I', or an individual correlation product (like, e.g., 
      'XX' or 'XY'). This function should not be called directly by the user.
  """

        self.printInfo(
            '\n\n WARNING: WRITING TO MOSAICS IS EXPERIMENTAL AND MAY NOT WORK!\n'
        )

        for vi, v in enumerate(self.vis):

            # Get the columns of parallel-hand correlations:
            success = ms.open(v, nomodify=False)
            if not success:
                self.printError('ERROR! %s cannot be openned in write mode' %
                                (v))
                return False

            spws = map(int, self.iscan[v].keys())
            ms.selectinit(datadescid=spws[0])
            polprods = [
                x[0] for x in list(ms.range(['corr_names'])['corr_names'])
            ]
            if self.stokes in polprods:
                polii = [polprods.index(self.stokes)]
            elif self.stokes in ['PI', 'I']:
                if 'XX' in polprods:
                    polii = [polprods.index('XX'), polprods.index('YY')]
                elif 'RR' in polprods:
                    polii = [polprods.index('RR'), polprods.index('LL')]
                else:
                    self.printError(
                        'Stokes not understood for %s! Will not update the model column!'
                        % (v))
                    ms.close()
                    return False
            else:
                self.printError(
                    'Stokes not understood for %s! Will not update the model column!'
                    % (v))
                ms.close()
                return False

            if self.write_model == 1:
                column = 'model_data'
            elif self.write_model == 2:
                column = 'corrected_data'

            for sp in spws:
                for scan in map(int, self.iscan[v][sp].keys()):
                    sys.stdout.write('\r Doing spw %i, scan_id %i, in %s' %
                                     (sp, scan, v))
                    sys.stdout.flush()
                    ms.selectinit(datadescid=sp)
                    success = ms.select({'scan_number': scan})

                    moddata = ms.getdata([column], ifraxis=True)
                    tmod = moddata[column]

                    if len(tmod) > 0:
                        modcol = True
                        tmod[:] = 0.0
                    else:
                        msg = '\n There\'s no requested column in %s\nfor scan %i of spw %i' % (
                            v, scan, sp)
                        msg += ' You should run \'clearcal\' with \'addmodel\' = True\'\n'
                        self.printError(msg)
                        modcol = False

                    if modcol:
                        info = self.iscan[v][sp][scan]
                        crosscorr = info[-1]
                        nnu = np.shape(self.mymodel.output[info[0]])[1]
                        nb, nt = info[2:4]

                        re = np.reshape(self.mymodel.output[info[0]][info[
                            1]:info[1] + nb * nt, :], (nt, nb, nnu))
                        # im = np.reshape(self.mymodel.output[info[0]][1][info[1]:info[1]+nb*nt,:],(nt,nb,nnu))

                        for poli in polii:
                            for nui, r in enumerate(info[4]):
                                moddata[column][
                                    poli, r,
                                    crosscorr, :] = re[:, :,
                                                       nui]  #+1.j*im[:,:,nui]) #[np.newaxis,:,:]
                #       moddata[column][poli,r,:nt,:nb] = (re[:,:,nui]+1.j*im[:,:,nui])[np.newaxis,:,:]

                        ms.putdata(moddata)

            ms.close()
            self.printInfo("\n %s written successfully!\n" % column)

############################################
#
#  SANITY CHECKS AND DEFINING BASIC PARAMETERS
#

    def checkOrdinaryInputs(self):
        """ Performs some sanity checks on the input parameters.
      This function should not be called directly by the user."""

        if type(self.model) is str:
            self.model = list([self.model])
        else:
            self.model = list(self.model)

        if type(self.var) is str:
            self.var = list([self.var])
        else:
            self.var = list(self.var)

        if type(self.fixedvar) in [str, float]:
            self.fixedvar = list([self.fixedvar])
        else:
            self.fixedvar = list(self.fixedvar)

        if type(self.fixed) is str:
            self.fixed = list([self.fixed])
        else:
            self.fixed = list(self.fixed)

        if type(self.p_ini) is [float, str]:
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

        if type(self.phase_gains) is not dict:
            self.printError('\n \'phase_gains\' must be a dictionary!')
        else:
            for key in self.phase_gains.keys():
                if type(key) is not int:
                    self.printError(
                        '\n The keys of \'phase_gains\' must be integers!')

        if type(self.amp_gains) is not dict:
            self.printError('\n \'amp_gains\' must be a dictionary!')
            for key in self.amp_gains.keys():
                if type(key) is not int:
                    self.printError(
                        '\n The keys of \'amp_gains\' must be integers!')

#  if len(self.phase_gains.keys())==0 and len(self.amp_gains.keys())==0:
#     self.useGains = False
#  else:
#     self.useGains = True
        self.useGains = list(
            set(self.phase_gains.keys() + self.amp_gains.keys()))

        # Check phase center:

        if type(self.phase_center) is not str:
            self.printError('\n \'phase_center\' must be a string!')
        else:
            if len(self.phase_center) == 0:
                self.phrefset = False
            else:
                # if True:
                try:
                    self.phrefset = True
                    csys = cs.newcoordsys(direction=True)
                    dirstr = self.phase_center.split()
                    if len(dirstr) == 2:
                        csys.setdirection(
                            refcode="J2000", refval=self.phase_center)
                    else:
                        csys.setdirection(
                            refcode=dirstr[0], refval=" ".join(dirstr[1:]))
                    csys.convertdirection("J2000")
                    self.refpos = np.copy(
                        csys.torecord()['direction0']['crval'])
            #  else:
                except:
                    self.printError(
                        '\n \'phase_center\' is not a CASA-formatted sky coordinate!'
                    )

# Did the user forget how does this task work? :D
        for param in self.var:
            if type(param) is not str:
                self.printError('\n \'var\' must be a list of strings!')
                return False

# Get the number of parameters and check the model consistency:
        maxpar = 0
        lf = [
            'GaussLine', 'GaussLine', 'LorentzLine', 'LorentzLine', 'power',
            'maximum', 'minimum'
        ]  # Typical function that can be used.
        for i, component in enumerate(self.model):
            vic = self.var[i].count
            checkpars = self.var[i].split(',')
            nfuncs = sum(
                map(vic, lf)
            )  # 2*(self.var[i].count('GaussLine')+self.var[i].count('LorentzLine'))
            paridx = zip(
                [m.start() + 1 for m in re.finditer('\[', self.var[i])],
                [m.start() for m in re.finditer('\]', self.var[i])])
            maxpar = max([maxpar] + map(
                float, [self.var[i][ss[0]:ss[1]] for ss in paridx]))
            if not component in self.implemented:
                msg = '\nModel component \'' + str(
                    component) + '\' is not known!\n'
                fmt = 'Supported models are:' + ' \'%s\' ' * len(
                    self.implemented)
                msg += fmt % tuple(self.implemented)
                self.printError(msg)
                return False
            if (component in self.implemented and self.numvar[component] !=
                (len(checkpars) - nfuncs) != 6):
                self.printError(
                    'Wrong number of variables (%i) in \'%s\' model.' %
                    (len(checkpars) - nfuncs, component))
                return False

# Scalefix must be a string representing a function:
        if not (type(self.scalefix) is str):
            self.printError('\'scalefix\' should be a string!')
            return False

        #print self.scalefix
# Get the overal number of parameters (i.e., from the variables AND the scalefix string):
        paridx = zip([m.start() + 1 for m in re.finditer('\[', self.scalefix)],
                     [m.start() for m in re.finditer('\]', self.scalefix)])
        maxpar = max([maxpar] +
                     map(float, [self.scalefix[ss[0]:ss[1]] for ss in paridx]))

        # Get the model parameters for the gains:
        #  self.maxGainTerm = 0
        #  self.NgainAnts = {}
        for key in self.phase_gains:
            term = self.phase_gains[key]
            paridx = zip([m.start() + 1 for m in re.finditer('\[', term)],
                         [m.start() for m in re.finditer('\]', term)])
            maxpar = max(
                [maxpar] + map(float, [term[ss[0]:ss[1]] for ss in paridx]))

        for key in self.amp_gains:
            term = self.amp_gains[key]
            paridx = zip([m.start() + 1 for m in re.finditer('\[', term)],
                         [m.start() for m in re.finditer('\]', term)])
            maxpar = max(
                [maxpar] + map(float, [term[ss[0]:ss[1]] for ss in paridx]))

# The fixed model should contain CONSTANT variables:
        for fcomp in self.fixedvar:
            for pi, param in enumerate(fcomp.split(',')):
                try:
                    temp = float(param)
                except:
                    self.printError(
                        '\n Fixed variables must be a list of floats (or strings\n representing floats)!'
                    )
                    return False

    # print self.fixedvar

# There should be as many p_inis as parameters!
        if len(self.p_ini) != maxpar + 1:
            self.printError(
                '\n \'p_ini\' is of length %i, but there are %i parameters used!\n ABORT!'
                % (len(self.p_ini), maxpar + 1))
            return False

# Check for consistency in the fixed model:
        for i, component in enumerate(self.fixed):
            checkpars = self.fixedvar[i].split(',')
            if not component in self.implemented:
                msg = '\nModel component \'' + str(
                    component) + '\' is not known!\n'
                fmt = 'Supported components are ' + ' \'%s\' ' * len(
                    self.implemented)
                msg += fmt % tuple(self.implemented)
                self.printError(msg)
                return False
            if (component in self.implemented
                    and self.numvar[component] != len(checkpars) != 6):
                self.printError(
                    'Wrong number of variables (%i) in \'%s\' fixed model.' %
                    (len(checkpars), component))
                return False

# Set outfile name:
        if self.outfile == '':
            self.printInfo(
                '\nSetting \'outfile\' to its default (i.e., \'modelfit.dat\').\n'
            )
            self.outfile = 'modelfit.dat'

# Set (and check) the bounds in the fit:
        if (self.bounds is None) or (len(self.bounds) == 0):
            self.bounds = None
        else:
            for b, bound in enumerate(self.bounds):
                if bound == None:
                    self.bounds[b] = [None, None]
                if self.bounds[b][0] != None and self.p_ini[b] <= self.bounds[b][0]:
                    self.printError(
                        'Initial value (%.2e) of parameter %i is lower (or equal) than its lower boundary (%.2e)!'
                        % (self.p_ini[b], b, self.bounds[b][0]))
                    return False
                if self.bounds[b][1] != None and self.p_ini[b] >= self.bounds[b][1]:
                    self.printError(
                        'Initial value (%.2e) of parameter %i is larger (or equal) than its upper boundary (%.2e)!'
                        % (self.p_ini[b], b, self.bounds[b][1]))
                    return False

        if (self.bounds is not None) and (len(self.bounds) != len(self.p_ini)):
            self.printError(
                'Length of \'bounds\' list (%i) is not equal to number of parameters (%i)!'
                % (len(self.bounds), len(self.p_ini)))
            return False

        return True

############################################
#
#  MORE SANITY CHECKS
#

    def checkInputs(self, data_changed=False):
        """ Reads all the inputs related to data selection, checks for 
      self-consistency, and sets a number of internal variables
      for the correct data reading, averaging, and formatting.
      It's good to run this before fitting, if the user has changed 
      some parameters to redo a fit. """

        import re

        # Some preliminary (self-consistency) checks of the parameters:

        self.savemodel = True

        # As the function says, check ordinary inputs :)
        success = self.checkOrdinaryInputs()
        if not success:
            return False

# We always work with lists here:
        if type(self.vis) is str:
            self.vis = list([self.vis])
        else:
            self.vis = list([str(ss) for ss in self.vis])

# Set the right format for the spectral window(s):
        if type(self.spw) is str:
            self.spw = list([self.spw])
        elif len(self.spw) > 1:
            self.spw = list([str(ss) for ss in self.spw])
            if self.OneFitPerChannel:
                self.printInfo(
                    '\n\n SPW is a LIST! User BEWARE! Any fit in *spectral mode* \n WILL fit the model to each MS separately!\n\n'
                )
        elif len(self.spw) > 0:
            self.spw = list([self.spw[0]])
        else:
            self.printError('\nBad formatted spw!\n')
            return False

# Set list of scans:
        try:
            if type(self.scans) is list and len(self.scans) == 0:
                self.scans = [[] for v in self.vis]

# We work with lists of nvis elements:
            self.scans = list(self.scans)
            if type(self.scans[0]) is int:
                if len(self.vis) == 1:
                    self.scans = [self.scans]
                else:
                    self.printError(
                        '\n\n \'scans\' should be a list of integers (or a list of lists of integers,\nif there are several measurement sets).\n ABORTING!'
                    )
                    return False
            if len(self.scans) != len(self.vis):
                self.printError(
                    '\n\n List of (lists of) scans does not have the same length as the list of measurement sets!\n ABORTING!'
                )
                return False

            for si, sc in enumerate(self.scans):
                if type(sc) is str:
                    self.scans[si] = map(int, sc.split(','))

        except:

            # print self.scans
            self.printError(
                '\n\n \'scans\' should be a list of integers (or a list of lists of integers,\nif there are several measurement sets).\n ABORTING!'
            )
            return False

# Check dimensions of vis, spw, model, etc.:

        if len(self.vis) != len(self.spw) and len(self.spw) > 1:
            self.printError(
                '\n\nThe length of \'spw\' is not equal to the length of \'vis\'!\n Aborting!\n'
            )
            return False

        if not (type(self.stokes) is str):
            self.printError('\'stokes\' must be a string!')
            return False

        if type(self.only_flux) is not bool:
            self.printError('\'only_flux\' must be a boolean!')

        if self.only_flux:
            if len(self.p_ini) != len(self.model):
                self.printError(
                    'If only_flux=True, number of parameters must be equal to number of model components!\n Aborting!\n'
                )

        if self.proper_motion == 0.0:
            self.proper_motion = [[0., 0.] for i in self.model]
        elif type(self.proper_motion) is not list:
            self.printError('\'proper_motion\' must be a list!')
        else:
            if len(self.proper_motion) != len(self.model):
                self.printError(
                    'The length of \'proper_motion\' must be equal to the number of model components!'
                )
            for pi in self.proper_motion:
                if type(pi) is not list:
                    self.printError(
                        'The elements of \'proper_motion\' must be lists of two floats!'
                    )
                elif len(pi) != 2:
                    self.printError(
                        'The elements of \'proper_motion\' must be lists of two floats!'
                    )
                elif type(pi[0]) is not float or type(pi[1]) is not float:
                    self.printError(
                        'The elements of \'proper_motion\' must be lists of two floats!'
                    )

# Do the mss exist?
# print self.vis, len(self.vis)
        for visi in self.vis:
            if not os.path.exists(visi):
                self.printError("\nMeasurement set %s does not exist!" %
                                (visi))
                return False

# Does the required column exist?
        if self.column not in ['data', 'corrected_data', 'corrected']:
            self.printError(
                '\n\'column\' can only take values \'data\' or \'corrected\'!')
            return False
        if self.column == 'corrected':
            self.column = 'corrected_data'

        self.pointing = []
        self.sourscans = []
        phasedirs = [{} for v in self.vis]
        self.field_id = []

        # Open MS and look for the selected data:
        for vi, v in enumerate(self.vis):
            #   phasedirs.append([])
            self.field_id.append([])
            success = ms.open(v)
            #   phasedirs[vi] = {}
            if not success:
                self.printError('Failed to open measurement set ' + v + '!')
                return False

            allfields = list(ms.range('fields')['fields'])

            try:
                if type(self.field) in [str, int]:
                    fitest = int(self.field)
                else:
                    fitest = int(self.field[vi])
                phasedirs[vi][fitest] = ms.range('phase_dir')['phase_dir'][
                    'direction'][:, fitest]
                self.field_id[-1].append(fitest)
            except:

                if type(self.field) is str:
                    aux = str(self.field)
                    self.field = [aux for v in self.vis]

                for f, field in enumerate(allfields):
                    if self.field[vi] in field:
                        self.field_id[-1].append(f)
                        phasedirs[vi][f] = ms.range('phase_dir')['phase_dir'][
                            'direction'][:, f]
                if len(self.field_id[-1]) == 0:
                    self.printError('\nERROR! Field %s is not in %s' %
                                    (self.field[vi], v))
                    ms.close()
                    return False

            ms.close()

    # print phasedirs
        self.phasedirs = phasedirs

        # Ref. position:
        if not self.phrefset:
            self.printInfo('\nSetting phase center on first scan\n')
            self.refpos = self.phasedirs[0][min(self.phasedirs[0].keys())]
        else:
            self.printInfo(
                '\nSetting phase center on %s\n' % self.phase_center)

# Find out all the scans where this field id is observed:
        for vi, v in enumerate(self.vis):
            success = ms.open(v)
            if not success:
                self.printError('Failed to open measurement set ' + v + '!')
                return False

            info = ms.getscansummary()
            self.sourscans.append([])
            self.pointing.append([])
            for key in info.keys():
                if info[key]['0']['FieldId'] in self.field_id[vi]:
                    self.sourscans[-1].append(int(key))
                for fieldid in info[key].keys():
                    myfield = info[key][fieldid]['FieldId']
                    if myfield in self.field_id[vi]:
                        fi = self.field_id[vi].index(myfield)
                        #    self.sourscans[-1].append(int(key))
                        #   print fi
                        self.pointing[-1].append(phasedirs[vi][myfield])

            self.pointing[-1] = np.array(self.pointing[-1])

        for vi, v in enumerate(self.vis):
            if len(self.scans[vi]) > 0:
                goodscid = filter(lambda x: x in self.sourscans[vi],
                                  self.scans[vi])
                if len(goodscid) != len(self.scans[vi]):
                    badscid = filter(lambda x: x not in self.sourscans[vi],
                                     self.scans[vi])
                    msg = '\n ERROR! The following scans do NOT correspond to source %s:\n' % (
                        str(self.field))
                    msg += str(badscid)
                    self.printError(msg)
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
                self.printError('Failed to open measurement set ' + v + '!')
                return False

            freqdic = ms.getspectralwindowinfo()
            spwchans = ms.range(['num_chan'])['num_chan']

            aux = channeler(
                self.spw[j], width=self.chanwidth, maxchans=spwchans)
            if aux[0]:
                ranges = list(aux[1])
            else:
                self.printError(
                    aux[1] +
                    '\nSomething seems to be wrong with the \'spw\' number %i. \n ABORTING!'
                    % (vi))
                return False

            nspws = range(len(spwchans))

            # These are the spws with selected channels for the fit:
            selspws = list(
                [list([i, ranges[i]]) for i in nspws if len(ranges[i]) > 0])
            self.spwlist.append(list(
                [j, vi, selspws,
                 spwi]))  # spwlist[vi][2] es una lista con [spwid,chanranges]
            # spwlist[vi][3]+index(spwlist[vi][2]) es la spw "local"
            spwi += {True: 0, False: len(selspws)}[len(self.spw) == 1]

            #################
            # Deal with polarization:

            ms.selectinit(datadescid=selspws[0][0])
            polprods = [
                x[0] for x in list(ms.range(['corr_names'])['corr_names'])
            ]
            self.pol2aver.append(np.zeros(len(polprods)))
            self.polmod.append(
                0
            )  # 0: normal,   1: multiply by i,   2: pol. independent, 3: just one product

            if self.stokes not in ['I', 'Q', 'U', 'V'] + polprods:
                self.printError('\n Bad Stokes parameter %s' % self.stokes)
                return False

        # print '\n\nSELECTED: ',self.stokes, polprods

            if self.stokes in polprods:
                self.pol2aver[-1][polprods.index(self.stokes)] = 1.0
                self.polii.append([polprods.index(self.stokes)])
                self.polmod[-1] = 3
        #   print '\n\nSELECTED: ',polprods.index(self.stokes)
            elif 'RR' in polprods:
                try:
                    if self.stokes == 'I':
                        self.polii.append(
                            [polprods.index('RR'),
                             polprods.index('LL')])
                        self.pol2aver[-1][polprods.index('RR')] = 0.5
                        self.pol2aver[-1][polprods.index('LL')] = 0.5
                    if self.stokes == 'PI':
                        self.polii.append(
                            [polprods.index('RR'),
                             polprods.index('LL')])
                        self.pol2aver[-1][polprods.index('RR')] = 0.5
                        self.pol2aver[-1][polprods.index('LL')] = 0.5
                        self.polmod[-1] = 2
                    if self.stokes == 'Q':
                        self.polii.append(
                            [polprods.index('RL'),
                             polprods.index('LR')])
                        self.pol2aver[-1][polprods.index('RL')] = 0.5
                        self.pol2aver[-1][polprods.index('LR')] = 0.5
                    if self.stokes == 'U':
                        self.polii.append(
                            [polprods.index('RL'),
                             polprods.index('LR')])
                        self.pol2aver[-1][polprods.index('RL')] = 0.5
                        self.pol2aver[-1][polprods.index('LR')] = 0.5
                        self.polmod[-1] = 1
                    if self.stokes == 'V':
                        self.polii.append(
                            [polprods.index('RR'),
                             polprods.index('LL')])
                        self.pol2aver[-1][polprods.index('RR')] = 0.5
                        self.pol2aver[-1][polprods.index('LL')] = -0.5
                except:
                    self.printError('Cannot convert to ' + self.stokes + '!!')
                    return False
            elif 'XX' in polprods:
                #   print '\n\nHOLADOLA'
                try:
                    if self.stokes == 'I':
                        self.polii.append(
                            [polprods.index('XX'),
                             polprods.index('YY')])
                        self.pol2aver[-1][polprods.index('XX')] = 0.5
                        self.pol2aver[-1][polprods.index('YY')] = 0.5
                    if self.stokes == 'PI':
                        self.polii.append(
                            [polprods.index('XX'),
                             polprods.index('YY')])
                        self.pol2aver[-1][polprods.index('XX')] = 0.5
                        self.pol2aver[-1][polprods.index('YY')] = 0.5
                        self.polmod[-1] = 2
                    if self.stokes == 'Q':
                        self.polii.append(
                            [polprods.index('XX'),
                             polprods.index('YY')])
                        self.pol2aver[-1][polprods.index('XX')] = 0.5
                        self.pol2aver[-1][polprods.index('YY')] = -0.5
                    if self.stokes == 'U':
                        self.polii.append(
                            [polprods.index('XY'),
                             polprods.index('YX')])
                        self.pol2aver[-1][polprods.index('XY')] = 0.5
                        self.pol2aver[-1][polprods.index('YX')] = 0.5
                    if self.stokes == 'V':
                        self.polii.append(
                            [polprods.index('YX'),
                             polprods.index('XY')])
                        self.pol2aver[-1][polprods.index('YX')] = 0.5
                        self.pol2aver[-1][polprods.index('XY')] = -0.5
                        self.polmod[-1] = 1
                except:
                    self.printError('Cannot convert to ' + self.stokes + '!!')
                    return False
            else:
                self.printError('Polarization ' + self.stokes +
                                ' not understood.\n ABORTING!')
                return False
###################

#  print self.polii, self.pol2aver, self.polmod
            ms.close()

# Try to compile the equations for the variables:
#  if data_changed:
        try:
            del self.mymodel
        except:
            self.printInfo('UVMULTIFIT model does not seem to exist yet.')
            pass

        self.mymodel = modeler(
            self.model, self.var, self.fixed, self.fixedvar, self.scalefix,
            self.NCPU, self.only_flux, self.applyHankel, self.isNumerical,
            self.useGains, [self.phase_gains, self.amp_gains])
        #else:
        #  self.mymodel.var = self.var
        #  self.mymodel.model = self.model
        #  self.mymodel.fixed = self.fixed
        #  self.mymodel.fixedvar = self.fixedvar
        #  self.mymodel.scalefix = self.scalefix
        #  self.mymodel.calls = 0
        #  self.mymodel.removeFixed=False

        if self.mymodel.failed:
            self.printError(self.mymodel.resultstring)
            return False

        if not data_changed:
            self.initData()

        self.setEngineWgt()
        return True

############################################
#
#  PREPARE THE UVMOD LIBRARY AND SET THE WEIGHTING
#

    def setEngineWgt(self):
        """Prepare the PB correction and the model. 
     Not to be called directly by the user."""

        # Load the C++ library:
        if self.mymodel.Ccompmodel is None:
            self.mymodel.Ccompmodel = uvmod.modelcomp
            self.QuinnFF = uvmod.QuinnFF
        self.mymodel.NCPU = self.NCPU
        self.success = self.setWgtEq()

        return self.success

############################################
#
#  SET THE WEIGHTING AND PB-CORRECTION
#

    def setWgtEq(self):
        """Depends on setEngineWgt. Not to be called by the user."""

        # Ref. position:
        #  if not self.phrefset:
        #    self.printInfo('\nSetting phase center on first scan\n')
        #    self.refpos = self.phasedirs[0][min(self.phasedirs[0].keys())]
        #  else:
        #    self.printInfo('\nSetting phase center on %s\n'%self.phase_center)

        tempfloat = 0.0

        try:
            self.dish_diameter = float(self.dish_diameter)
        except:
            self.printError(
                "\n The dish diameter must be a number! (in meters)\n")
            return False

        if self.pbeam:
            self.printInfo(
                "\n\n You selected to apply primary-beam correction.\nPLEASE, remember that the beam is being approximated \nwith a Gaussian, so it may not be very accuracte far \nfrom the pointing direction.\n\n"
            )
            if self.dish_diameter == 0.0:
                try:
                    tb.open(os.path.append(self.vis[0], 'ANTENNA'))
                    tempfloat = np.average(tb.getcol('DISH_DIAMETER'))
                    tb.close()
                except:
                    pass
            else:
                self.printInfo(
                    "\n\n An antenna diameter of %.3f m will be applied\n" %
                    self.dish_diameter)
                tempfloat = self.dish_diameter

            if tempfloat == 0.0:
                self.printError(
                    '\nThe antenna diameters are not set in the ms. \n Please, set it manually or turn off primary-beam correction.\n'
                )
                return False
            else:
                FWHM = self.ldfac / tempfloat * (2.99e8)
                sigma = FWHM / 2.35482 * (180. / np.pi) * 3600.
                self.mymodel.KfacWgt = 1. / (
                    2. * sigma**2.)  #  (0.5*(tempfloat/1.17741)**2.)
        else:
            self.mymodel.KfacWgt = 0.0

# May refine this function in future releases:
#  self.mymodel.wgtEquation = lambda D,Kf: -D*Kf

        return True

############################################
#
#  READ THE DATA. ARRANGE ALL ARRAYS
#

    def readData(self, data_changed=True):
        """ Reads the data, according to the properties \'vis\', \'column\',
      \'chanwidth\', etc. It then fills in the properties \'averdata\',
      \'averfreqs\', \'averweights\', \'u\', \'v\', \'w\'... Each one of 
      these properties is a list with the data (one list item per 
      spectral window/scan). A previous successful run of function 
      \'checkInputs()\' is assumed.
      Instead of re-reading data from scratch using the same uvmultifit
      instance, a better approach may be to restart CASA and create a 
      fresh uvmultifit instance with the new data, avoiding some 
      memory leakage related to potential hidden references to the data 
      in the IPython \'recall\' prompt."""

        tic = time.time()

        if data_changed:
            self.deleteData()

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

            for vidx, vis in enumerate(
                    filter(lambda x: x[3] <= si,
                           self.spwlist)):  # BEWARE! was sp
                for spidx, spi in enumerate(vis[2]):
                    if vis[3] + spidx == si:

                        sp = spi[0]
                        rang = spi[1]
                        msname = self.vis[vis[1]]
                        self.iscan[msname][sp] = {}

                        self.printInfo(
                            '\n\n Opening measurement set ' + msname + '.\n')
                        success = ms.open(msname)
                        if not success:
                            self.printError('\nFailed reading data!!')
                            return False

# For the first ms in the list, read the frequencies of the spw.
# All the other mss will be assumed to have the same frequencies:
                        if True:
                            ms.selectinit(datadescid=int(sp))
                            origfreqs = ms.range('chan_freq')['chan_freq'][:,
                                                                           0]
                            self.averfreqs[si] = np.array(
                                [np.average(origfreqs[r]) for r in rang])
                            nfreq = len(rang)

                        self.printInfo('\n\n Reading scans for spw %i \n' % sp)

                        # Read all scans for this field id:

                        for sc, scan in enumerate(self.sourscans[vis[1]]):

                            self.printInfo('\r Reading scan #%i (%i of %i).' %
                                           (scan, sc + 1,
                                            len(self.sourscans[vis[1]])))

                            ms.selectinit(datadescid=int(sp))
                            success = ms.select({'scan_number': int(scan)})
                            if success:
                                fieldids = list(
                                    set(ms.getdata(['field_id'])['field_id']))
                                for fieldid in fieldids:
                                    ms.selectinit(datadescid=int(sp))
                                    success = ms.select({
                                        'scan_number':
                                        int(scan),
                                        'field_id':
                                        int(fieldid)
                                    })

                                    uvscan = ms.getdata(
                                        [
                                            'u', 'v', 'w', 'antenna1',
                                            'antenna2', 'time'
                                        ],
                                        ifraxis=True)
                                    #   baselines = [b.split('-') for b in uvscan['axis_info']['ifr_axis']['ifr_name']]

                                    # This doesn't seem to work with simobs data:
                                    #         crosscorr = np.array([b[0]!=b[1] for b in baselines])

                                    # This should be harmless (the fitter will flag autocorrs anyway:
                                    crosscorr = uvscan['antenna1'] != uvscan['antenna2']  #np.array([True for b in baselines])

                                    datascan = ms.getdata(
                                        [self.column, 'weight', 'flag'],
                                        ifraxis=True)

                                    # NOTE: There is a bug in np.ma.array that casts complex to float under certain operations (e.g., np.ma.average).
                                    # That's why we average real and imag separately.

                                    ########################################
                                    # Compute the polarization product:
                                    copyweight = np.repeat(
                                        datascan['weight'][:,
                                                           np.newaxis, :, :],
                                        np.shape(datascan['flag'])[1],
                                        axis=1)
                                    totalmask = np.logical_or(
                                        datascan['flag'][:, :, crosscorr, :],
                                        copyweight[:, :, crosscorr, :] < 0.0)

                                    origmasked = np.ma.array(
                                        datascan[self.column][:, :,
                                                              crosscorr, :],
                                        mask=totalmask,
                                        dtype=np.complex128)
                                    # The weights re weighting the RESIDUALS, and not the ChiSq terms. Hence, we divide wgt_power by 2.:
                                    origweight = np.ma.power(
                                        np.ma.array(
                                            copyweight[:, :, crosscorr, :],
                                            mask=totalmask),
                                        self.wgt_power / 2.)
                                    datamask = 0.0
                                    weightmask = 0.0
                                    flagmask = 0.0

                                    polavg = [
                                        pol != 0.0
                                        for pol in self.pol2aver[vis[1]]
                                    ]

                                    if self.polmod[vis[1]] == 2:
                                        flagmask = np.ma.logical_and(
                                            totalmask[self.polii[vis[1]][
                                                0], :, :, :],
                                            totalmask[self.polii[vis[1]][
                                                1], :, :, :])
                                        datamask = np.ma.average(
                                            origmasked[polavg, :].real, axis=0
                                        ) + 1.j * np.ma.average(
                                            origmasked[polavg, :].imag, axis=0)
                                        weightmask = np.ma.sum(
                                            origweight[polavg, :], axis=0)
                                    else:
                                        if self.polmod[vis[1]] == 3:
                                            flagmask = np.ma.logical_or(
                                                False, totalmask[self.polii[
                                                    vis[1]][0], :, :, :])
                                        else:
                                            flagmask = np.ma.logical_or(
                                                totalmask[self.polii[vis[1]][
                                                    0], :, :, :],
                                                totalmask[self.polii[vis[1]][
                                                    1], :, :, :])
                                        for pol in self.polii[vis[1]]:
                                            datamask += origmasked[
                                                pol, :] * self.pol2aver[vis[
                                                    1]][pol]
                                            weightmask += origweight[pol, :]
                                        if self.polmod[vis[1]] == 1:
                                            datamask *= 1.j

                                    weightmask[flagmask] = 0.0

                                    #########################################

                                    # Free some memory:
                                    del datascan, origmasked, origweight
                                    del copyweight, totalmask, flagmask

                                    phshift = 3600. * 180. / np.pi * (
                                        self.phasedirs[vis[1]][fieldid] -
                                        self.refpos)
                                    strcos = np.cos(
                                        self.phasedirs[vis[1]][fieldid][1])

                                    if phshift[0] != 0.0 or phshift[1] != 0.0:
                                        self.printInfo(
                                            '\r\t\t\t\t Offset: %.2e RA (tsec) %.2e Dec (asec).'
                                            % (phshift[0] / 15., phshift[1]))

# Average spectral channels:
                                    nchannel, nbaselines, ntimes = np.shape(
                                        datamask)
                                    ntav = int(
                                        max([
                                            1,
                                            round(
                                                float(ntimes) / self.timewidth)
                                        ]))
                                    datatemp = np.ma.zeros(
                                        (nfreq, nbaselines, ntimes),
                                        dtype=np.complex128)
                                    weighttemp = np.ma.zeros(
                                        (nfreq, nbaselines, ntimes))

                                    if self.chanwidth == 1:  # and False:
                                        concRan = [c[0] for c in rang]
                                        datatemp[:, :] = datamask[concRan, :]
                                        weighttemp[:, :] = weightmask[
                                            concRan, :]
                                    else:
                                        for nu in range(nfreq):
                                            datatemp[nu, :] = np.ma.average(
                                                datamask[rang[nu], :].real,
                                                axis=0) + 1.j * np.ma.average(
                                                    datamask[rang[nu], :].imag,
                                                    axis=0)
                                            weighttemp[nu, :] = np.ma.sum(
                                                weightmask[rang[nu], :],
                                                axis=0)

# Average in time and apply uvtaper:
                                    GaussWidth = 2. * (
                                        self.uvtaper / 1.17741)**2.
                                    nbas = np.sum(crosscorr)
                                    ui = np.zeros((nbas, ntav))
                                    vi = np.zeros((nbas, ntav))
                                    wi = np.zeros((nbas, ntav))
                                    ant1i = np.zeros(
                                        (nbas, ntav), dtype=np.int32)
                                    ant2i = np.zeros(
                                        (nbas, ntav), dtype=np.int32)
                                    timei = np.zeros((nbas, ntav))
                                    tArrayi = np.zeros(ntav)
                                    tIndexi = np.zeros(
                                        (nbas, ntav), dtype=np.int32)
                                    RAoffi = np.zeros((nbas, ntav))
                                    Decoffi = np.zeros((nbas, ntav))
                                    Stretchi = np.zeros((nbas, ntav))
                                    avercompl = np.ma.zeros(
                                        (nfreq, nbas, ntav),
                                        fill_value=0,
                                        dtype=np.complex128)
                                    averwgt = np.ma.zeros(
                                        (nfreq, nbas, ntav), fill_value=0)

                                    if self.timewidth == 1:  # and False:

                                        ui[:, :] = uvscan['u'][crosscorr, :]
                                        vi[:, :] = uvscan['v'][crosscorr, :]
                                        wi[:, :] = uvscan['w'][crosscorr, :]
                                        ant1i[:, :] = uvscan['antenna1'][
                                            crosscorr][:, np.newaxis]
                                        ant2i[:, :] = uvscan['antenna2'][
                                            crosscorr][:, np.newaxis]
                                        timei[:, :] = uvscan['time'][
                                            np.newaxis, :]
                                        tArrayi[:] = uvscan['time']
                                        tIndexi[:, :] = np.arange(
                                            ntav,
                                            dtype=np.int32)[np.newaxis, :]
                                        RAoffi[:, :] = float(phshift[0])
                                        Decoffi[:, :] = float(phshift[1])
                                        Stretchi[:, :] = float(strcos)

                                        avercompl[:, :, :] = datatemp
                                        if self.uvtaper > 0:
                                            averwgt[:, :, :] = weighttemp * np.exp(
                                                -(ui * ui + vi * vi) /
                                                GaussWidth)
                                        else:
                                            averwgt[:, :, :] = weighttemp

                                    else:

                                        ant1s = uvscan['antenna1'][crosscorr]
                                        ant2s = uvscan['antenna2'][crosscorr]
                                        for nt in range(ntav):
                                            t0 = nt * self.timewidth
                                            t1 = min([
                                                ntimes,
                                                (nt + 1) * self.timewidth
                                            ])
                                            uu = uvscan['u'][crosscorr, t0:t1]
                                            vv = uvscan['v'][crosscorr, t0:t1]
                                            ww = uvscan['w'][crosscorr, t0:t1]
                                            ui[:, nt] = np.average(
                                                uu,
                                                axis=1)  # Baseline dimension
                                            vi[:, nt] = np.average(
                                                vv,
                                                axis=1)  # Baseline dimension
                                            wi[:, nt] = np.average(
                                                ww,
                                                axis=1)  # Baseline dimension
                                            ant1i[:,
                                                  nt] = ant1s  # Baseline dimension
                                            ant2i[:,
                                                  nt] = ant2s  # Baseline dimension
                                            timei[:, nt] = np.average(
                                                uvscan['time'][t0:t1]) / 86400.
                                            tArrayi[nt] = time[
                                                0,
                                                nt]  #np.average(uvscan['time'][t0:t1])/86400.
                                            tIndexi[:, nt] = nt
                                            RAoffi[:, nt] = float(phshift[0])
                                            Decoffi[:, nt] = float(phshift[1])
                                            Stretchi[:, nt] = float(strcos)

                                            if self.uvtaper > 0.0:
                                                #   uvdist = uu*uu + vv*vv
                                                GaussFact = np.exp(
                                                    -(uu * uu + vv * vv) /
                                                    GaussWidth)
                                            else:
                                                GaussFact = np.ones(
                                                    np.shape(uu))

                                            broadwgt = weighttemp[:, :, t0:t1]
                                            avercompl[:, :, nt] = np.ma.average(
                                                datatemp[:, :, t0:t1].real,
                                                axis=2,
                                                weights=broadwgt
                                            ) + 1.j * np.ma.average(
                                                datatemp[:, :, t0:t1].imag,
                                                axis=2,
                                                weights=broadwgt)
                                            if self.uniform:
                                                averwgt[:, :, nt] = np.ma.sum(
                                                    np.ones(
                                                        np.shape(broadwgt)) *
                                                    GaussFact[
                                                        np.newaxis, :, :],
                                                    axis=2)
                                            else:
                                                averwgt[:, :, nt] = np.ma.sum(
                                                    broadwgt * GaussFact[
                                                        np.newaxis, :, :],
                                                    axis=2)

                                    BtT = nbas * ntav

                                    # Flatten the arrays in the baseline-time dimensions and transpose.
                                    # ms.getdata gives us (pol,freq,baseline,time). We averaged in pol, and we want (baseline*time,freq).
                                    # The reason of the later, is that having freq at the end will easy the concatenation of all the arrays
                                    # of each scan (since rows go first!):
                                    # UV coordinates, visibilities, and weights:

                                    ant1scan.append(ant1i.reshape(BtT))
                                    ant2scan.append(ant2i.reshape(BtT))
                                    uscan.append(ui.reshape(BtT))
                                    vscan.append(vi.reshape(BtT))
                                    wscan.append(wi.reshape(BtT))
                                    tscan.append(timei.reshape(BtT))
                                    tArray.append(tArrayi)

                                    if len(tIndex) > 1:
                                        #     print sc, tIndex, np.shape(tIndex)
                                        tIndexi += tIndex[-1][-1] + 1
                                    tIndex.append(tIndexi.reshape(BtT))

                                    RAscan.append(RAoffi.reshape(BtT))
                                    Decscan.append(Decoffi.reshape(BtT))
                                    Stretchscan.append(Stretchi.reshape(BtT))
                                    avercomplflat = avercompl.reshape(
                                        (nfreq, nbas * ntav))
                                    weightflat = averwgt.reshape((nfreq,
                                                                  nbas * ntav))
                                    datascanAv.append(
                                        np.transpose(avercomplflat))
                                    #      datascanim.append(np.transpose(avercomplflat.imag))
                                    weightscan.append(
                                        np.transpose(
                                            weightflat.reshape((nfreq,
                                                                nbas * ntav))))

                                    # Release memory:
                                    del ant1i, ant2i, ui, vi, wi, timei, RAoffi, Decoffi, Stretchi, averwgt, avercomplflat, tIndexi, tArrayi
                                    del avercompl, weightflat
                                    if self.timewidth > 1:
                                        del broadwgt, uu, vv, ww

# Useful info for function writeModel() and for pointing correction:
                                    self.iscan[msname][sp][scan] = [
                                        int(si),
                                        int(i0scan),
                                        int(ntav),
                                        int(nbas),
                                        list(rang),
                                        np.copy(crosscorr)
                                    ]
                                    self.iscancoords[si].append([
                                        i0scan, i0scan + ntav * nbas,
                                        phshift[0], phshift[1]
                                    ])

                                    i0scan += nbas * ntav
                    #    print '\n\n',sp,scan,self.iscan[msname][sp][scan], i0scan

                        ms.close()

# Concatenate all the scans in one single array. Notice that we separate real and imag and save them
# as floats. This is because ctypes doesn' t handle complex128.

#   print si
            self.averdata[si] = np.concatenate(
                datascanAv, axis=0)  #,np.concatenate(datascanim,axis=0)]
            self.averweights[si] = np.concatenate(weightscan, axis=0)
            self.u[si] = np.concatenate(uscan, axis=0)
            self.v[si] = np.concatenate(vscan, axis=0)
            self.w[si] = np.concatenate(wscan, axis=0)
            self.t[si] = np.concatenate(tscan, axis=0)
            self.tArr[si] = np.concatenate(tArray, axis=0)
            self.tIdx[si] = np.concatenate(tIndex, axis=0)

            self.RAshift[si] = np.concatenate(RAscan, axis=0)
            self.Decshift[si] = np.concatenate(Decscan, axis=0)
            self.Stretch[si] = np.concatenate(Stretchscan, axis=0)
            self.ant1[si] = np.concatenate(ant1scan, axis=0)
            self.ant2[si] = np.concatenate(ant2scan, axis=0)

            # Free some memory:
            #   del uu,vv,ww,avercomplflat,weightflat
            del datatemp, weighttemp, uvscan  #, avercompl,averwgt
            for dda in datascanAv:
                del dda
#   for dda in datascanim:
#     del dda
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

            del datascanAv  #, datascanim
            del weightscan, tscan, uscan, vscan, wscan, tArray, tIndex
            del RAscan, Decscan, Stretchscan, ant1scan, ant2scan

            try:
                del GaussFact
            except:
                pass

            gc.collect()

# Initial time of observations (reference for proper motions):
        self.t0 = np.min([np.min(ti) for ti in self.t])

        self.Nants = 1 + np.max(
            [np.max(self.ant1[IF]) for IF in range(len(self.ant1))] +
            [np.max(self.ant2[IF]) for IF in range(len(self.ant2))])

        self.printInfo('\n\n Done reading\n')
        tac = time.time()
        self.printInfo('\nReading took %.2f seconds.\n' % (tac - tic))

        self.success = True

        # if self.clib and self.mymodel.Ccompmodel is not None:
        #   self.printInfo('Freeing pointers\n')
        #   uvmod.del_pointers()

        if data_changed:
            self.initData()

        return True

############################################
#
#  COMPUTE MODELS TO BE DIRECTLY SUBTRACTED FROM THE DATA
#

    def computeFixedModel(self):
        """ Computes the value of the fixed model on all the u-v data points. It saves
       the results in the \'mymodel.output\' property, which should have been zeroed
       previously. Notice that the fixed model is recomputed for each spectral channel 
       (i.e., if OneFitPerChannel=True), and is computed only once if OneFitPerChannel==False."""
        self.printInfo(
            '\nGoing to compute fixed model (may need quite a bit of time)\n')
        self.mymodel.residuals([0], mode=0)

############################################
#
#  COMPUTE CHI SQUARE FOR A PARTICULAR REALIZATION OF THE FIT
#

    def chiSquare(self, p):
        """ Returns a list with 2 items: The Chi Square value, computed at point \'p\' in the
       parameter space, and the number of degrees of freedom."""
        return self.mymodel.residuals(p, mode=-2)

############################################
#
#  SET SOME DATA (AND MODEL) ARRAYS
#

    def initData(self, isUVfit=True):
        """ Initiates a \'modeler\' instance, and stores it in the property \'mymodel\'.
       The \'modeler\' class stores the compiled model (and fixedmodel), the model values, 
       and all the methods to compute the residuals, the Chi Squared, etc."""

        # Maximum number of frequency channels (i.e., the maximum from all the selected spws):
        self.maxnfreq = 0

        NIFs = len(self.averdata)
        self.Nspw = NIFs

        gooduvm = uvmod.setNspw(int(NIFs))

        if gooduvm != 0:
            self.printError('\nError in the C++ extension!\n')
            return False

# Fill-in proper motions in as/day
        for i in range(len(self.model)):
            self.mymodel.propRA[i] = self.proper_motion[i][0] / 365.
            self.mymodel.propDec[i] = self.proper_motion[i][1] / 365.

        for spidx in range(NIFs):

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

            self.mymodel.freqs.append(self.averfreqs[spidx])
            self.maxnfreq = max(self.maxnfreq, len(self.averfreqs[spidx]))

            # Only have to multiply by freq, to convert these into lambda units:
            ulambda = self.FouFac * self.u[spidx] / self.LtSpeed
            vlambda = self.FouFac * self.v[spidx] / self.LtSpeed
            wlambda = self.FouFac * self.w[spidx] / self.LtSpeed

            # Data, uv coordinates, and weights:

            # Data are saved in two float arrays (i.e., for real and imag):
            self.mymodel.data[-1] = self.averdata[
                spidx]  #[0],self.averdata[spidx][1]]
            self.mymodel.wgt[-1] = self.averweights[spidx]
            self.mymodel.uv[-1] = list([ulambda, vlambda, wlambda])
            self.mymodel.offset[-1] = list([
                self.RAshift[spidx], self.Decshift[spidx], self.Stretch[spidx]
            ])
            self.mymodel.ants[-1] = list([self.ant1[spidx], self.ant2[spidx]])

            self.mymodel.Nants = self.Nants
            self.mymodel.isGain[-1] = np.zeros(
                len(self.t[spidx]), dtype=np.int8)

            for i in self.useGains:
                mask0 = self.mymodel.ants[-1][0] == i
                mask1 = self.mymodel.ants[-1][1] == i
                self.mymodel.isGain[-1][np.logical_or(mask0, mask1)] = True

# Release memory:
            del ulambda, vlambda, wlambda

            # Time spent on observation:
            self.mymodel.dt[-1] = self.t[spidx] - self.t0
            self.mymodel.t0 = self.t0
            self.mymodel.dtArr[-1] = self.tArr[spidx] - self.t0
            self.mymodel.dtIdx[-1] = self.tIdx[spidx]

            # Array to save the residuals (or model, or any output from the C library):
            self.mymodel.output[-1] = np.zeros(
                np.shape(self.averdata[spidx]), dtype=np.complex128)  #,
            # np.zeros(np.shape(self.averdata[spidx][0]))]

            #    self.mymodel.fixedmodel[-1] = [np.zeros(np.shape(self.averdata[spidx][0])),
            #                                   np.zeros(np.shape(self.averdata[spidx][0]))]

            ########
            # Array of booleans, to determine if a datum enters the fit:
            self.mymodel.fittable[-1] = np.ones(
                np.shape(self.mymodel.uv[-1][0]), dtype=np.int8)
            self.mymodel.wgtcorr[-1] = np.zeros(
                tuple([len(self.model)]) + np.shape(self.t[spidx]))

            self.mymodel.fittablebool[-1] = np.ones(
                np.shape(self.mymodel.uv[-1][0]), dtype=np.bool)

            # for i in range(NIFs):
            #  dimsarr = np.array(np.shape(self.mymodel.data[-1]),dtype=np.int32)
            #    print 'NoAutoCorr: ',np.where(self.mymodel.uv[-1][0]!=0.0)[0][:10]

            gooduvm = uvmod.setData(
                spidx, self.mymodel.uv[-1][0], self.mymodel.uv[-1][1],
                self.mymodel.uv[-1][2], self.mymodel.wgt[-1],
                self.mymodel.data[-1], self.mymodel.output[-1],
                self.mymodel.freqs[-1], self.mymodel.fittable[-1],
                self.mymodel.wgtcorr[-1], self.mymodel.dt[-1],
                self.mymodel.dtArr[-1], self.mymodel.dtIdx[-1],
                self.mymodel.offset[-1][0], self.mymodel.offset[-1][1],
                self.mymodel.offset[-1][2], self.mymodel.ants[-1][0],
                self.mymodel.ants[-1][1], self.mymodel.isGain[-1], self.Nants)

            #   print np.shape(self.mymodel.ampGainBuffer[-1]),np.shape(self.mymodel.dt[-1])

            if gooduvm != 10:
                self.printError('\nError in the C++ extension!\n')
                return False

############################################
#
#  SET MODEL ARRAYS TO BE SENT TO UVMOD
#

    def initModel(self):
        """ Allocates memory for the modeler data, which will be used by the C++ extension.
       Also compiles the model variables. Not to be called directly by the user."""
        #####
        # Array to save the variables of the model and the 'scalefix' value,
        # all of them as a function of frequency. The C library will read the variables from here
        # at each iteration and for each model component:
        #   self.mymodel.Nterms = np.zeros(self.Nants,dtype=np.int32)

        #   for ni in self.NgainAnts.keys():
        #     self.mymodel.Nterms[ni] = self.NgainAnts[ni]

        self.mymodel.varbuffer = [
            np.zeros((len(self.model), self.maxNvar + self.applyHankel,
                      self.maxnfreq)) for i in range(len(self.p_ini) + 1)
        ]
        #  print 'VARBUFFER 1 ',type(self.mymodel.varbuffer),type(self.mymodel.varbuffer[0])
        self.mymodel.varfixed = [
            np.zeros(self.maxnfreq) for i in range(len(self.p_ini) + 1)
        ]
        self.mymodel.dpar = np.zeros(len(self.p_ini), dtype=np.float64)
        self.mymodel.par2 = np.zeros((3, len(self.p_ini)), dtype=np.float64)
        #####

        #  self.printInfo("Going to run setNCPU\n")
        #  gooduvm = uvmod.setNCPU(int(self.NCPU))
        #  if gooduvm != 0:
        #    self.printError('\nError in the C++ extension!\n')
        #    return False

        self.mymodel.Hessian = np.zeros(
            (len(self.p_ini), len(self.p_ini)), dtype=np.float64)
        self.mymodel.Gradient = np.zeros(len(self.p_ini), dtype=np.float64)
        self.mymodel.imod = np.zeros(len(self.model), dtype=np.int32)
        #  self.mymodel.modGrid = np.zeros((len(self.model),self.HankelOrder),dtype=np.int32)

        self.printInfo("Going to compile models\n")
        self.mymodel.compileAllModels()

        #  NIFs = len(self.averdata); NANT = len(self.mymodel.parDependence)

        #  print 'NANT: ',NANT, self.Nants, np.max(self.ant1),np.max(self.ant2)

        self.mymodel.GainBuffer = [[[] for AI in range(self.Nants)]
                                   for spidx in range(self.Nspw)]
        #  self.mymodel.phaseGainBuffer = [[[] for AI in range(self.Nants)] for spidx in range(self.Nspw)]

        for spidx in range(self.Nspw):
            for AI in range(self.Nants):
                self.mymodel.GainBuffer[spidx][AI] = [
                    np.ones(
                        (len(self.tArr[spidx]), len(self.averfreqs[spidx])),
                        dtype=np.complex128)
                    for i in range(len(self.mymodel.parDependence[AI]))
                ]
        #  self.mymodel.phaseGainBuffer[spidx][AI] =  [np.zeros((len(self.tArr[spidx]),len(self.averfreqs[spidx]))) for i in range(len(self.mymodel.parDependence[AI]))]

        self.printInfo("Going to run setModel\n")
        gooduvm = uvmod.setModel(
            self.mymodel.imod, self.mymodel.Hessian, self.mymodel.Gradient,
            self.mymodel.varbuffer, self.mymodel.varfixed, self.mymodel.dpar,
            self.mymodel.propRA, self.mymodel.propDec, self.refpos,
            self.mymodel.parDependence, self.mymodel.GainBuffer
        )  #, self.mymodel.phaseGainBuffer ) #, self.mymodel.modGrid)
        self.mymodel.minnum = self.minDp
        self.printInfo("Success!\n")

        #  print 'VARBUFFER 2 ',type(self.mymodel.varbuffer),type(self.mymodel.varbuffer[0])

        self.printInfo("Success!\n")
        if gooduvm != 10:
            self.printError('\nError in the C++ extension!\n')
            return False

        return True

############################################
#
#  PERFORM (AND SAVE) THE FIT
#

    def fit(self, redo_fixed=True, reinit_model=True, save_file=True,
            nuidx=-1):
        """ Fits the data, using the models compiled in the \'mymodel\' property (which is an instance
      of the \'modeler\' class). If \'reinit_model\'==False (it is True by default), the models
      used are those already compiled in the \'mymodel\' property. If not, the models are 
      recompiled, according to the contents of the \'model\', \'var\', \'fixed\', and 
      \'fixedvar\' properties, and all the references to the data arrays are refreshed.

      If \'redo_fixed\'==False (it is True by default), the fixed model will not be recomputed 
      throughout the fit (THIS CAN BE DANGEROUS IF YOU ARE FITTING IN SPECTRAL-LINE MODE AND HAVE 
      DIFFERENT CHANNELS WITH QUITE DIFFERENT FREQUENCIES, SINCE THE U-V POINTS WILL *NOT* BE 
      REPROJECTED).

      Finally, if \'save_file\'==False (it is True by default) the external file with fitting 
      results will not be created.

      \'nuidx\' is a helper parameter for internal use. The user should not redefine it.

"""

        # TEST CODE:
        #  self.scipylm=False

        tic = time.time()

        self.mymodel.bounds = self.bounds
        self.mymodel.LMtune = self.LMtune

        npars = len(self.p_ini)
        nspwtot = self.spwlist[-1][3] + len(self.spwlist[-1][2])

        notfit = [[] for si in range(nspwtot)]

        # Select data according to time range:
        datatot = 0
        self.allflagged = False
        for si in range(nspwtot):

            if self.MJDrange[0] > 0.0 and self.MJDrange[1] > 0.0:
                self.printInfo('\nSelecting data by Modified Julian Date\n')
                self.mymodel.fittablebool[si][:] = np.logical_and(
                    self.t[si] >= self.MJDrange[0],
                    self.t[si] <= self.MJDrange[1])
                self.mymodel.fittable[si][:] = self.mymodel.fittablebool[
                    si][:].astype(np.int8)

            else:

                self.mymodel.fittable[si][:] = 1
                self.mymodel.fittablebool[si][:] = 1

        # Check if there is data available:
            unflagged = np.sum(
                self.mymodel.wgt[si][self.mymodel.fittablebool[si], :] != 0.0,
                axis=0)
            ntot = np.sum(self.mymodel.fittablebool[si])
            if self.OneFitPerChannel:
                if np.sum(unflagged == 0.0) > 0:
                    self.printInfo(
                        'ERROR: NOT ENOUGH DATA FOR THIS TIME RANGE! \n CHANNELS: '
                        + str(list(np.where(unflagged == 0.0))))
                    self.allflagged = True
                    notfit[si] = list(np.where(unflagged == 0.0)[0])
                #  return False
            else:
                datatot += np.sum(unflagged)

            self.mymodel.output[si][:] = 0.0

    #   self.mymodel.output[si][:] = 0.0

        if datatot == 0 and not self.OneFitPerChannel:
            self.printInfo('ERROR: NOT ENOUGH DATA FOR THIS TIME RANGE!\n')
            self.allflagged = True
            return False

        if self.first_time:
            self.first_time = False
            gooduvm = uvmod.setNCPU(int(self.NCPU))
            if gooduvm != 0:
                self.printError('\nError in the C++ extension!\n')
                return False

        for si in range(nspwtot):
            self.mymodel.wgtcorr[si][:] = -self.mymodel.KfacWgt

        self.printInfo('\nNow, fitting model \n')

        # Initialize model:
        if reinit_model:
            goodinit = self.initModel()
            #  print goodinit, self.mymodel.failed
            if self.mymodel.failed or goodinit == False:
                self.printError('ERROR: Bad model initialization! \n%s' %
                                self.mymodel.resultstring)

                return False

# Allocate memory for the threads:
# if self.first_time:
#   self.first_time = False
# else:
#   uvmod.unsetWork()

        goodinit = uvmod.setWork()
        if goodinit != 10:
            self.printError('Memory allocation error!')

        for i in range(len(self.p_ini) + 1):
            #   print '\n',i, type(self.mymodel.varbuffer)
            #   print type(self.mymodel.varbuffer[i]),'\n'
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
                print ''

                for nuidx in range(rang):
                    self.printInfo('\r Fitting channel ' + str(nuidx + 1) +
                                   ' of ' + str(rang) + ' in spw ' + str(si))

                    self.mymodel.currspw = si
                    self.mymodel.currchan = nuidx

                    # Compute fixed model:
                    if redo_fixed and len(self.mymodel.fixed) > 0:
                        self.computeFixedModel()

                    # Fit with simplex (if asked for):

                    if nuidx not in notfit[si]:

                        if self.method == 'simplex':
                            fitsimp = _mod_simplex(
                                self.mymodel.ChiSquare,
                                self.p_ini,
                                args=(self.bounds, self.p_ini),
                                disp=False,
                                relxtol=self.SMPtune[0],
                                relstep=self.SMPtune[1],
                                maxiter=self.SMPtune[2] * len(self.p_ini))
                            fit = [
                                fitsimp[0],
                                np.zeros((len(self.p_ini), len(self.p_ini))),
                                fitsimp[1]
                            ]

                        else:
                            fit = self.mymodel.LMMin(self.p_ini)
                            if not fit:
                                return False
                    #   print '\n ',si,fit
# Estimate the parameter uncertainties and save the model in the output array:
                        if self.savemodel:
                            if self.write_model == 1:
                                Chi2t = self.mymodel.residuals(fit[0], mode=-3)
                            elif self.write_model == 2:
                                Chi2t = self.mymodel.residuals(fit[0], mode=-4)

# DON'T FIT THIS CHANNEL (ALL DATA FLAGGED)
                    else:
                        fit = [[0.0 for pi in self.p_ini],
                               np.zeros((len(self.p_ini), len(self.p_ini))), 0]

                    fitparams[si].append([float(f) for f in fit[0]])

                    ndata = float(
                        np.sum(self.mymodel.wgt[si][:, nuidx] > 0.0)
                    )  #  Only add the unflagged data to compute the DoF

                    if ndata > 0.0:
                        ChiSq[si].append(fit[2] / ndata)  # Reduced ChiSquared
                    else:
                        ChiSq[si].append(
                            float(fit[2])
                        )  # There are 0 'really-free' parameters?! Watch out!

                    Nvis[si].append(ndata)

                    fiterrors[si].append([
                        np.sqrt(fit[1][i, i] * ChiSq[si][nuidx])
                        for i in range(npars)
                    ])
                    covariance[si].append(fit[1] * ChiSq[si][nuidx])

        #   print fitparams[si]
                self.fitpars = fitparams[si]

##################

##################
# CASE OF CONTINUUM-MODE FIT:
        else:

            self.printInfo('\nFitting to all frequencies at once.\n')

            # This will tell the modeller to solve in continuum mode:
            self.mymodel.currspw = -1
            self.mymodel.currchan = -1

            # Compute fixed model:
            if redo_fixed and len(self.mymodel.fixed) > 0:
                print 'Generating fixed model. May take some time'
                self.computeFixedModel()
                print 'Done!'

            # Pre-fit with simplex:
            if self.method == 'simplex':
                fitsimp = _mod_simplex(
                    self.mymodel.ChiSquare,
                    self.p_ini,
                    args=(self.bounds, self.p_ini),
                    disp=False,
                    relxtol=self.SMPtune[0],
                    relstep=self.SMPtune[1],
                    maxiter=self.SMPtune[2] * len(self.p_ini))
                fit = [
                    fitsimp[0],
                    np.zeros((len(self.p_ini), len(self.p_ini))), fitsimp[1]
                ]

            else:

                # Bound least-squares fitting:
                fit = self.mymodel.LMMin(self.p_ini)
                if not fit:
                    return False
######
# Estimate the parameter uncertainties and save the model in the output array:

            if self.savemodel:
                if self.write_model == 1:
                    Chi2t = self.mymodel.residuals(fit[0], mode=-3)
                elif self.write_model == 2:
                    Chi2t = self.mymodel.residuals(fit[0], mode=-4)

            fitparams = fit[0]

            for si in range(nspwtot):
                ndata += float(
                    np.sum(self.mymodel.wgt[si] > 0.0)
                )  #  Only add the unflagged data to compute the DoF

            if fit[2] > 0.0:
                ChiSq = fit[2] / ndata  # Reduced chi squared.
            else:
                ChiSq = fit[
                    2]  # There are 0 'really-free' parameters?! Watch out!!
            Nvis = ndata

            fiterrors = [np.sqrt(fit[1][i, i] * ChiSq) for i in range(npars)]
            covariance = fit[1] * ChiSq

######
##################

        self.printInfo(
            '\n The reduced Chi Squared will be set to 1 by re-scaling the visibility weights.\n'
        )

        # Free some memory:
        gc.collect()

        #####
        # Set the 'result' property:

        if not self.OneFitPerChannel:
            to_return = {
                'Frequency': self.averfreqs[0][0],
                'Parameters': np.array(fitparams),
                'Uncertainties': np.array(fiterrors),
                'Reduced Chi squared': ChiSq,
                'Fit': fit,
                'Degrees of Freedom': Nvis
            }
            prtpars = []
            for pp, ppar in enumerate(fitparams):
                prtpars.append(ppar)
                prtpars.append(fiterrors[pp])

        else:
            Freq = []  #{}
            Par = []  #{}
            Err = []  #{}
            Chi2 = []  #{}
            prtpars = [[[] for spi in range(len(fitparams[sp]))]
                       for sp in range(nspwtot)]

            for sp in range(nspwtot):
                Freq.append(self.averfreqs[sp])
                Par.append(np.array(fitparams[sp]))
                Err.append(np.array(fiterrors[sp]))
                Chi2.append(np.array(ChiSq[sp]))
                for pp, ppar in enumerate(fitparams[sp]):
                    for ppi, ppari in enumerate(ppar):
                        prtpars[sp][pp].append(ppari)
                        prtpars[sp][pp].append(fiterrors[sp][pp][ppi])

            to_return = {
                'Frequency': Freq,
                'Parameters': Par,
                'Uncertainties': Err,
                'Reduced Chi squared': Chi2,
                'Fit': fit,
                'Degrees of Freedom': Nvis
            }

        if self.cov_return:
            to_return['covariance'] = covariance

        self.result = to_return

        #####
        # Save results in external file:
        if not save_file:
            return True

        self.printInfo('\n DONE FITTING. Now, saving to file.\n')

        outf = open(self.outfile, 'w')
        outf.write('# MODELFITTING RESULTS FOR MS: ' + ','.join(self.vis))
        outf.write(
            '\n\n# NOTICE THAT THE REDUCED CHI SQUARE SHOWN HERE IS THE VALUE\n'
        )
        outf.write('# *BEFORE* THE RE-SCALING OF THE VISIBILITY WEIGHTS.\n')
        outf.write(
            '# HOWEVER, THE PARAMETER UNCERTAINTIES ARE *ALREADY* GIVEN\n')
        outf.write('# FOR A REDUCED CHI SQUARE OF 1.\n')
        if type(Nvis) is list:
            DOG = int(np.average(Nvis))
        else:
            DOG = Nvis
        outf.write('# AVG. NUMBER OF DEGREES OF FREEDOM: %i' % DOG)
        if self.pbeam:
            outf.write(
                '# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. USING A DISH DIAMETER OF: %.3f METERS'
                % self.dish_diameter)
        else:
            outf.write('# PRIMARY-BEAM CORRECTION HAS NOT BEEN APPLIED.')

        outf.write('\n\n###########################################\n')
        outf.write('#\n# MODEL CONSISTS OF:\n#')
        for m, mod in enumerate(self.model):
            outf.write('\n# \'' + mod + '\' with variables: ' + self.var[m])
        if len(self.fixed) > 0:
            outf.write('\n#\n#\n# FIXED MODEL CONSISTS OF:\n#')
            for m, mod in enumerate(self.fixed):
                var2print = map(float, self.fixedvar[m].split(','))
                outf.write('\n# \'' + mod + '\' with variables: ' + ' '.
                           join(['% .3e'] * len(var2print)) % tuple(var2print))
            outf.write('\n#\n#  - AND SCALING FACTOR: %s' % self.scalefix)

        outf.write('\n#\n#\n# INITIAL PARAMETER VALUES:\n#')
        for p0i, p0 in enumerate(self.p_ini):
            if self.bounds is not None:
                outf.write('\n#  p[%i] = %.5e with bounds: %s ' %
                           (p0i, p0, str(self.bounds[p0i])))
            else:
                outf.write('\n#  p[%i] = %.5e with no bounds' % (p0i, p0))

        outf.write('\n#\n##########################################\n\n')
        parshead = []

        for pp in range(len(self.p_ini)):
            parshead.append(pp)
            parshead.append(pp)

        headstr = (
            '# Frequency (Hz)   ' + 'p[%i]  error(p[%i])   ' * len(self.p_ini)
            + 'Red. Chi Sq.\n') % tuple(parshead)
        outf.write(headstr)

        if not self.OneFitPerChannel:
            formatting = "%.12e   " + "% .4e " * (2 * npars) + "   %.4e \n"
            toprint = tuple([np.average(self.averfreqs)] + prtpars + [ChiSq])
            outf.write(formatting % toprint)

        else:
            formatting = "%.12e    " + "% .4e " * (2 * npars) + "  %.4e \n"

            for spwvi in self.spwlist:
                for r, rr in enumerate(spwvi[2]):
                    k = spwvi[3] + r
                    rang = rr[1]
                    for nu, freq in enumerate(self.averfreqs[k]):
                        toprint = tuple(
                            [freq] + prtpars[k][nu] + [ChiSq[k][nu]])
                        outf.write(formatting % toprint)

        outf.close()

        tac = time.time()

        self.printInfo('\n Fit took %.2f seconds.\n\n END!\n' % (tac - tic))
        #  uvmod.unsetWork()

        return True


#####

############################################
#
# MODELER CLASS (CALLED BY UVMULTIFIT)
#


class modeler(object):
    """ Class to convert strings representing models and parameters
      into an actual equation to be used in a ChiSq visibility fitting. 
      It also interacts with the C++ extension. 
      Not to be instantiated by the user. """

    ############################################
    #
    #  FREE MEMORY JUST FOR THE MODEL-RELATED DATA:
    #

    def cleanModels(self):

        for mdspw in self.GainBuffer:
            NA = len(mdspw)
            for a in range(NA):
                NP = len(mdspw[a])
                for mdp in range(NP - 1, -1, -1):
                    del mdspw[a][mdp]
                del mdspw[a]
            del mdspw
        del self.GainBuffer

        #    for mdspw in self.phaseGainBuffer:
        #      NA = len(mdspw)
        #      for a in range(NA):
        #        NP = len(mdspw[a])
        #        for mdp in range(NP-1,-1,-1):
        #          del mdspw[a][mdp]
        #        del mdspw[a]
        #      del mdspw
        #    del self.phaseGainBuffer

        for mdi in self.varbuffer:
            del mdi
        del self.varbuffer

        for mdi in self.varfixed:
            del mdi
        del self.varfixed

        for mdi in self.dpar:
            del mdi
        del self.dpar

        for mdi in self.par2:
            del mdi
        del self.par2

        del self.Hessian, self.Gradient, self.imod

############################################
#
#  FREE MEMORY
#

    def __del__(self):

        for mdi in self.data:
            del mdi
        del self.data

        for mdi in self.wgt:
            del mdi
        del self.wgt

        for mdi in self.uv:
            del mdi[2], mdi[1], mdi[0]
            del mdi
        del self.uv

        for mdi in self.offset:
            del mdi[2], mdi[1], mdi[0]
            del mdi
        del self.offset

        for mdi in self.ants:
            del mdi[1], mdi[0]
            del mdi
        del self.ants

        del self.Nants

        for mdspw in self.GainBuffer:
            NA = len(mdspw)
            for a in range(NA):
                NP = len(mdspw[a])
                for mdp in range(NP - 1, -1, -1):
                    del mdspw[a][mdp]
                del mdspw[a]
            del mdspw
        del self.GainBuffer

        #    for mdspw in self.phaseGainBuffer:
        #      NA = len(mdspw)
        #      for a in range(NA):
        #        NP = len(mdspw[a])
        #        for mdp in range(NP-1,-1,-1):
        #          del mdspw[a][mdp]
        #        del mdspw[a]
        #      del mdspw
        #    del self.phaseGainBuffer

        for mdi in self.dt:
            del mdi
        del self.dt

        for mdi in self.dtArr:
            del mdi
        del self.dtArr

        for mdi in self.dtIdx:
            del mdi
        del self.dtIdx

        for mdi in self.output:
            del mdi
        del self.output

        for mdi in self.fittable:
            del mdi
        del self.fittable

        for mdi in self.wgtcorr:
            del mdi
        del self.wgtcorr

        for mdi in self.fittablebool:
            del mdi
        del self.fittablebool

        for mdi in self.isGain:
            del mdi
        del self.isGain

        for mdi in self.varbuffer:
            del mdi
        del self.varbuffer

        for mdi in self.varfixed:
            del mdi
        del self.varfixed

        for mdi in self.dpar:
            del mdi
        del self.dpar

        for mdi in self.par2:
            del mdi
        del self.par2

        del self.Hessian, self.Gradient, self.imod

        for mdi in self.iscancoords:
            Npar = len(mdi)
            for mdp in range(Npar - 1, -1, -1):
                del mdi[mdp]
            del mdi
        del self.iscancoords

############################################
#
#  CREATE INSTANCE
#

    def __init__(self, model, parameters, fixedmodel, fixedparameters,
                 scalefix, NCPU, only_flux, HankelOrder, isNumerical, useGains,
                 gainFunction):
        """ Just the constructor."""

        import numpy as np

        # self.minnum = np.finfo(np.float64).eps  # Step for Jacobian computation
        self.propRA = np.zeros(len(model), dtype=np.float64)
        self.propDec = np.zeros(len(model), dtype=np.float64)
        self.addfixed = False
        self.expka = 2. * np.log(2.)
        self.pow2ka = 1. / (np.pi * np.log(2.))
        self.pow3ka = np.sqrt(2.**(2. / 3.) - 1.)
        self.FouFac = (2. * np.pi) * np.pi / 180. / 3600.
        self.LtSpeed = 2.99792458e+8
        self.deg2rad = np.pi / 180.
        self.failed = False
        self.calls = 0  # Total number of calls during the fit
        self.removeFixed = False  # Tells us if the fixedmodel array exists and should be used.
        self.fixed = fixedmodel
        self.fixedvar = fixedparameters
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
        self.model = model
        self.var = parameters  # variables of the model components.
        self.scalefix = scalefix
        self.HankelOrder = HankelOrder
        self.isNumerical = isNumerical
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
        # Lists of compiled functions (one per model component).
        # These will take p and nu and return the variables of the model:
        self.varfunc = [0.0 for component in model]
        self.fixedvarfunc = [0.0 for component in fixedmodel]
        # Buffer arrays to save the values of the variables:
        self.varbuffer = []
        # Model indices (to let the C++ library know which model is which component):
        self.imod = []
        self.ifixmod = []
        # spw and channel to fit (-1 means fit to the continuum):
        self.currspw = 0
        self.currchan = 0
        self.NCPU = NCPU
        self.t0 = 0.0
        # C++ library:
        self.Ccompmodel = None
        # Parameters computed in unbound space:
        self.par2 = []
        self.bounds = []
        self.only_flux = only_flux
        self.strucvar = []
        # Levenberg-Marquardt parameters:
        self.LMtune = []
        # Gains:
        self.GainBuffer = []
        #   self.phaseGainBuffer = []
        self.parDependence = []
        self.isGain = []
        self.useGains = useGains

        KGaus = np.sqrt(1. / (4. * np.log(2.)))
        # Some useful functions:
        self.LorentzLine = lambda nu,nu0,P,G: P*0.25*G*G/(np.power(nu-nu0,2.)+(0.5*G)**2.)
        self.GaussLine = lambda nu,nu0,P,G: P*np.exp(-np.power((nu-nu0)/(G*KGaus),2.))

        self.pieceWise = lambda t,p0,p1,t0,t1: np.clip(p0+(p1-p0)*(t-t0)/(t1-t0),p0,p1)

        self.wgtEquation = lambda D, Kf: -D * Kf
        self.KfacWgt = 1.0
        self.gainFunction = gainFunction

        # List of currently-supported model components:
        self.allowedmod = [
            'delta', 'Gaussian', 'disc', 'ring', 'sphere', 'bubble', 'expo',
            'power-2', 'power-3', 'GaussianRing'
        ]

        self.resultstring = ''

############################################
#
#  COMPILE AT RUNTIME
#

# Compile the functions to make p,nu -> var:

    def compileAllModels(self):
        """ Compile all models (fixed, variable, and fixed-scale. 
        Not to be called directly by the user. """
        self.resultstring = ''
        self.failed = False
        self.compileModel()
        self.compileFixedModel()
        self.compileScaleFixed()
        self.compileGains()

    def compileGains(self):

        self.phaseAntsFunc = [lambda t, nu, p: 0. for i in range(self.Nants)]
        self.ampAntsFunc = [lambda t, nu, p: 1. for i in range(self.Nants)]

        self.parDependence = [[0] for i in range(self.Nants)]

        for ni in self.gainFunction[0].keys():
            tempstr = self.gainFunction[0][ni].replace(
                'pieceWise(', 'self.pieceWise(t,').replace(
                    't', 't[:,__import__(\'numpy\').newaxis]').replace(
                        'nu0', '%.12f' % self.freqs[0][0])
            modstr = 'self.phaseAntsFunc[' + str(
                ni) + '] = lambda t,nu,p: ' + tempstr

            parpos = [x.start() for x in re.finditer('p\[', tempstr)]
            for p0 in parpos:
                p1 = tempstr[p0:].find(']') + p0
                self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

            try:
                exec modstr in locals()
            except:
                self.failed = True
                self.resultstring = 'Syntax error in phase gain of antenna %i' % (
                    ni)
                return

        for ni in self.gainFunction[1].keys():
            tempstr = self.gainFunction[1][ni].replace(
                'pieceWise(', 'self.pieceWise(t,').replace(
                    't', 't[:,__import__(\'numpy\').newaxis]').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            modstr = 'self.ampAntsFunc[' + str(
                ni) + '] = lambda t,nu,p: ' + tempstr

            parpos = [x.start() for x in re.finditer('p\[', tempstr)]
            for p0 in parpos:
                p1 = tempstr[p0:].find(']') + p0
                self.parDependence[ni].append(1 + int(tempstr[p0 + 2:p1]))

            try:
                exec modstr in locals()
            except:
                self.failed = True
                self.resultstring = 'Syntax error in amp gain of antenna %i' % (
                    ni)
                return

        for ni in range(self.Nants):
            self.parDependence[ni] = np.unique(self.parDependence[ni]).astype(
                np.int32)

    def compileModel(self):
        """ Compiles the variable model, according to the contents of the 
        \'model\' and \'parameters\' lists."""

        # Define variable model:
        import numpy as np

        self.varfunc = [0.0 for component in self.model]

        for ii, component in enumerate(self.model):
            tempstr = self.var[ii].replace(
                'LorentzLine(', 'self.LorentLine(nu,').replace(
                    'GaussLine(', 'self.GaussLine(nu,').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            if self.only_flux:
                try:
                    tempstr2 = tempstr.split(',')
                    if len(tempstr2) > 3:
                        self.strucvar.append(
                            map(float, tempstr2[:2] + tempstr2[3:]))
                    else:
                        self.strucvar.append(map(float, tempstr2[:2]))
                except:
                    print tempstr.split(',')
                    self.resultstring = '\n If only_flux=True, all variables but the flux must be constants! Aborting!'
                    self.failed = True

            modstr = 'self.varfunc[' + str(
                ii) + '] = lambda p,nu: [' + tempstr + ']'
            try:
                exec modstr in locals()
            except:
                self.failed = True
                self.resultstring = 'Syntax error in component number %i of the variable model' % (
                    ii)
                return

            if component not in self.allowedmod:
                self.resultstring = '\n Component \'' + component + '\' is unknown. Aborting!'
                self.failed = True
                return
            else:
                self.imod[ii] = self.allowedmod.index(component)

    def compileFixedModel(self):
        """ Compiles the fixed model, according to the contents of the 
        \'fixed\' and \'fixedpars\' lists."""

        # Define fixed model:
        import numpy as np
        self.fixedvarfunc = [0.0 for component in self.fixed]

        self.ifixmod = np.zeros(len(self.fixed), dtype=np.int32)  #[]

        for ii, component in enumerate(self.fixed):
            tempstr = self.fixedvar[ii].replace(
                'LorentzLine(', 'self.LorentzLine(nu,').replace(
                    'GaussLine(', 'self.GaussLine(nu,').replace(
                        'nu0', '%.12f' % self.freqs[0][0])

            modstr = 'self.fixedvarfunc[' + str(
                ii) + '] = lambda p,nu: [' + tempstr + ']'
            #   print '\n\n MODSTR: ',modstr
            try:
                #  if True:
                exec modstr in locals()
            except:
                self.resultstring = 'Syntax error in component number %i of the fixed model' % (
                    ii)
                self.failed = True
                return

            if component not in self.allowedmod:
                self.resultstring = '\n Component \'' + component + '\' is unknown. Aborting!'
                self.failed = True
                return
            else:
                self.ifixmod[ii] = self.allowedmod.index(component)

    def compileScaleFixed(self):
        """ Compiles the scaling factor for the fixed model """

        import numpy as np
        tempstr = self.scalefix.replace(
            'LorentzLine(', 'self.LorentzLine(nu,').replace(
                'GaussLine(', 'self.GaussLine(nu,').replace(
                    'nu0', '%.12f' % self.freqs[0][0])

        scalefixedstr = 'self.compiledScaleFixed = lambda p,nu: ' + tempstr + ' + 0.0'
        try:
            exec scalefixedstr in locals()
        except:
            self.resultstring = 'Syntax error in the flux-scale equation'
            self.failed = True
            return

############################################
#
#  GET UNBOUND PARAMETER SPACE FROM BOUND PARAMETERS
#

    def getPar2(self, mode=0):
        """ Function to change fitting parameters to/from the unbound 
       space from/to the bound space and comptue the gradient for 
       the Jacobian matrix. Unbound values are in index 0, 
       gradient is in index 1, and bound values are in index 2. 
       The changes of variables are taken from the MINUIT package."""

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
                if self.bounds[p2] is None or (self.bounds[p2][0] is None
                                               and self.bounds[p2][1] is None):
                    self.par2[:, p2] = [p[p2], 1.0, p[p2]]
                elif self.bounds[p2][0] is None and self.bounds[p2][1] is not None:
                    self.par2[0, p2] = np.sqrt(
                        np.power(self.bounds[p2][1] - p[p2] + 1, 2.) - 1.)
                elif self.bounds[p2][0] is not None and self.bounds[p2][1] is None:
                    self.par2[0, p2] = np.sqrt(
                        np.power(p[p2] - self.bounds[p2][0] + 1, 2.) - 1.)
                else:
                    Kbound = (self.bounds[p2][1] - self.bounds[p2][0]) / 2.
                    self.par2[0, p2] = np.arcsin(
                        (p[p2] - self.bounds[p2][0]) / Kbound - 1.0)

        else:
            p = self.par2[0, :]
            for p2 in range(len(p)):
                if self.bounds[p2] is None or (self.bounds[p2][0] is None
                                               and self.bounds[p2][1] is None):
                    self.par2[2, p2] = p[p2]
                elif self.bounds[p2][0] is None and self.bounds[p2][1] is not None:
                    self.par2[2, p2] = self.bounds[p2][1] + 1. - np.sqrt(
                        p[p2]**2. + 1.)
                elif self.bounds[p2][0] is not None and self.bounds[p2][1] is None:
                    self.par2[2, p2] = self.bounds[p2][0] - 1. + np.sqrt(
                        p[p2]**2. + 1.)
                else:
                    Kbound = (self.bounds[p2][1] - self.bounds[p2][0]) / 2.
                    self.par2[2, p2] = self.bounds[p2][0] + Kbound * (
                        np.sin(p[p2]) + 1.)

        for p2 in range(len(p)):
            if self.bounds[p2] is None or (self.bounds[p2][0] is None
                                           and self.bounds[p2][1] is None):
                self.par2[1, p2] = 1.0
            elif self.bounds[p2][0] is None and self.bounds[p2][1] is not None:
                self.par2[1, p2] = -self.par2[0, p2] / np.sqrt(
                    self.par2[0, p2]**2. + 1.)
            elif self.bounds[p2][0] is not None and self.bounds[p2][1] is None:
                self.par2[1, p2] = self.par2[0, p2] / np.sqrt(
                    self.par2[0, p2]**2. + 1.)
            else:
                Kbound = (self.bounds[p2][1] - self.bounds[p2][0]) / 2.
                self.par2[1, p2] = Kbound * np.cos(self.par2[0, p2])

############################################
#
#  LEVENBERG-MARQUARDT LEAST-SQUARE MINIMIZATION
#

    def LMMin(self, p):
        """ Implementation of the Levenberg-Marquardt algorithm. Not to be called directly 
        by the user. """

        NITER = int(self.LMtune[3] * len(p))
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

        #  print '\n par2:',self.par2

        #  print 'Allocating memory'
        Hessian2 = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        Gradient2 = np.zeros(np.shape(self.Gradient), dtype=np.float64)
        HessianDiag = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupHess = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupGrad = np.zeros(np.shape(self.Gradient), dtype=np.float64)
        Inverse = np.zeros(np.shape(self.Hessian), dtype=np.float64)
        backupP = np.copy(self.par2[0, :])

        self.Hessian[:, :] = 0.0
        self.Gradient[:] = 0.0
        #  print 'Memory allocated!'

        if self.only_flux and self.currchan == -1:
            nnu = max([len(self.freqs[sp]) for sp in range(len(self.freqs))])
            print '\nComputing structure-only parameters'

            for midx in range(len(p)):
                if len(self.strucvar[midx]) > 3:
                    tempvar = self.strucvar[midx][:2] + [
                        p[midx]
                    ] + self.strucvar[midx][2:]
                else:
                    tempvar = self.strucvar[midx] + [p[midx]]
                for i in range(len(tempvar)):
                    for j in range(len(p) + 1):
                        self.varbuffer[j][midx, i, :nnu] = tempvar[i]

    #   print 'VARBUFFER 2 ', type(self.varbuffer)
    #   print type(self.varbuffer[0]) #[j][midx,i,:nnu] = tempvar[i]

#   print 'Fitting!'

        CurrChi = self.residuals(self.par2[2, :], -1, dof=False)
        #  print 'CurrChi  ',CurrChi
        Hessian2[:, :] = self.Hessian * (self.par2[1, :])[np.newaxis, :] * (
            self.par2[1, :])[:, np.newaxis]
        Gradient2[:] = self.Gradient * self.par2[1, :]
        backupHess[:, :] = Hessian2
        backupGrad[:] = Gradient2

        #    print '\n\n'
        #    print self.par2
        #    print '\n'
        #    print self.Hessian
        #    print '\n'
        #    print self.Gradient
        #    raw_input('HOLD')

        controlIter = 0
        for i in range(NITER):
            controlIter += 1
            for n in range(len(p)):
                HessianDiag[n, n] = Hessian2[n, n]
            try:
                goodsol = True
                Inverse[:] = np.linalg.pinv(Hessian2 + Lambda * HessianDiag)
                Dpar = np.dot(Inverse, Gradient2)
                DirDer = sum([
                    Hessian2[n, n] * Dpar[n] * Dpar[n] for n in range(len(p))
                ])
                DirDer2 = np.sum(Gradient2 * Dpar)
                TheorImpr = DirDer - 2. * DirDer2
        #   print 'Did it!'
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
                #    print 'Chi2  ', Chi2
                Hessian2[:, :] = self.Hessian * (
                    self.par2[1, :])[np.newaxis, :] * (
                        self.par2[1, :])[:, np.newaxis]
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
            #   print '\nChis: ',Chi2, CurrChi
            if relchi < 1: relchi = 1. / relchi
            todivide = np.copy(backupP)
            todivide[todivide == 0.0] = 1.0
            totest = [
                np.abs(tt) for tt in self.par2[0, :] / todivide if tt != 0.0
            ]  #; totest[totest==0.0] = 1.0
            #     print todivide,self.par2[0,:]
            #     print '\n',relchi-1.,relpar-1.
            relpar = max([{True: 1. / pb, False: pb}[pb < 1] for pb in totest])
            if relchi - 1. < functol or relpar - 1. < partol:
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
            sys.stdout.write(
                "\n\n REACHED MAXIMUM NUMBER OF ITERATIONS!\nThe algorithm may not have converged!\nPlease, check if the parameter values are meaningful.\n"
            )
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

        n = self.HankelOrder - 1

        #    print tempvar[:7]

        if imod == 'GaussianRing':  # Gaussian Ring

            a = 1.0
            k = 1. / (
                2. * np.power(tempvar[6], 2.)) * tempvar[3] * tempvar[3] / 4.
            m = 2 * n + 1

            # Recurrence relation:
            merf = (1. - spec.erf(-np.sqrt(k) * a))

            m0 = np.sqrt(np.pi / k) / 2. * merf
            m1 = np.exp(-k * a * a) / 2. / k + np.sqrt(
                np.pi / k) * a / 2. * merf

            tempvar.append(np.ones(np.shape(m1)))  # term for n=0, normalized

            if m in [0, 1]:
                return tempvar  # Only order n=0.

            else:  # Higher orders.
                res = a * m1 + 1. / (2 * k) * m0  # order m=2
                resaux = np.copy(res)
                res2 = np.copy(m1)
                for mi in range(3, m + 1):
                    if type(res) is np.float64:
                        res = a * res + (mi - 1) / (2 * k) * res2
                        res2 = resaux
                        resaux = res
                    else:
                        res[:] = a * res + (mi - 1) / (2 * k) * res2
                        res2[:] = resaux
                        resaux[:] = res
                    if np.mod(mi + 1, 2) == 0:  # (i.e., orders n=1,2,3...)
                        tempvar.append(
                            res / m1 * np.power(-1., (mi - 1) / 2) / np.power(
                                np.math.factorial((mi - 1) / 2), 2.))

            return tempvar

        else:

            raise ValueError(
                '\n\nModel %i was not correctly interpreted!\n\n' % imod)

############################################
#
#  COMPUTE RESIDUALS FOR A MODEL REALIZATION
#

    def residuals(self, p, mode=-1, dof=True):
        """ Method with a wide usage:

        mode ==  0. Compute the fixed model. Fill-in the output array with it.    
        mode == -1. Compute the Hessian matrix and the Error vector. Return the Chi2.
        mode == -2. Only return the Chi2.
        mode == -3. Add the variable model to the output array.

        This method should not be called by the user directly.

    """

        #  varsize = self.maxNvar+self.HankelOrder

        if mode in [0, -3]:
            self.calls = 0
        else:
            self.calls += 1

        if self.currchan < 0:  # Continuum (i.e., all data modelled at once)
            spwrange = range(len(self.freqs))
            nui = -1
        else:  # Spectral mode (i.e., model only the current channel of the current spw)
            spwrange = [int(self.currspw)]
            nui = int(self.currchan)

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
                for i in self.useGains:
                    for j in range(len(self.parDependence[i])):
                        ptemp = [pi for pi in p]
                        if j > 0:
                            j2 = self.parDependence[i][j] - 1
                            ptemp[j2] += self.dpar[j2]
                        if nui == -1:
                            #         self.phaseGainBuffer[spw][i][j][:] = self.phaseAntsFunc[i](self.dtArr[spw],self.freqs[spw], ptemp)
                            self.GainBuffer[spw][i][j][:] = self.ampAntsFunc[
                                i](self.dtArr[spw],
                                   self.freqs[spw], ptemp) * np.exp(
                                       1.j * self.phaseAntsFunc[i](self.dtArr[
                                           spw], self.freqs[spw], ptemp))
                        else:
                            #             self.phaseGainBuffer[spw][i][j][:,nui] = np.squeeze(self.phaseAntsFunc[i](self.dtArr[spw],self.freqs[spw][nui], ptemp))
                            self.GainBuffer[spw][i][j][:, nui] = np.squeeze(
                                self.ampAntsFunc[i]
                                (self.dtArr[spw], self.freqs[spw][nui],
                                 ptemp)) * np.exp(1.j * np.squeeze(
                                     self.phaseAntsFunc[i](self.dtArr[
                                         spw], self.freqs[spw][nui], ptemp)))

        if mode == 0:  # Just compute the fixed model and return
            self.removeFixed = True
            isfixed = True
            modbackup = self.imod[0]
            for spw in spwrange:
                self.output[spw][:] = 0.0
                # self.output[spw][1][:] = 0.0
                for midx, mi in enumerate(self.ifixmod):
                    tempvar = self.fixedvarfunc[midx](p, self.freqs[spw])
                    self.imod[0] = self.ifixmod[midx]
                    if self.imod[0] in self.isNumerical:
                        tempvar = self.gridModel(self.imod[0], tempvar)
                #  print tempvar
                    for i in range(len(tempvar)):
                        nnu = len(self.freqs[spw])
                        self.varbuffer[0][0, i, :nnu] = tempvar[i]
                    self.Ccompmodel(spw, nui, 0)
            #    print 'VARBUFFER 3 ', type(self.varbuffer)
            #    print type(self.varbuffer[0]) #[j][midx,i,:nnu] = tempvar[i]

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
                scalefx = self.compiledScaleFixed(p, self.freqs[spw])
                self.varfixed[0][:nnu] = scalefx
                for j in range(len(p)):
                    ptemp = [pi for pi in p]
                    ptemp[j] += self.dpar[j]
                    scalefx = self.compiledScaleFixed(
                        ptemp,
                        self.freqs[spw])  # Variables of current component
                    self.varfixed[j + 1][:nnu] = scalefx
            else:
                scalefx = self.compiledScaleFixed(p, self.freqs[spw][nui])
                self.varfixed[0][0] = scalefx
                for j in range(len(p)):
                    ptemp = [pi for pi in p]
                    ptemp[j] += self.dpar[j]
                    scalefx = self.compiledScaleFixed(
                        ptemp,
                        self.freqs[spw][nui])  # Variables of current component
                    self.varfixed[j + 1][0] = scalefx

            for midx, modi in enumerate(currmod):

                ptemp = [float(pi) for pi in p]

                if nui == -1:  # Continuum

                    if self.only_flux:
                        currflux = p[midx]
                        nstrucpars = len(self.strucvar[midx])
                        #  for i in range(len(p)+1):
                        self.varbuffer[0][midx, 2, :nnu] = currflux
        #       print 'VARBUFFER 4 ', type(self.varbuffer)
        #       print type(self.varbuffer[0]) #[j][midx,i,:nnu] = tempvar[i]

                    else:
                        # Variables of current component
                        tempvar = currvar[midx](p, self.freqs[spw])
                        if modi in self.isNumerical:
                            tempvar = self.gridModel(modi, tempvar)

                        for i in range(len(tempvar)):
                            self.varbuffer[0][midx, i, :nnu] = tempvar[i]

            #     print 'VARBUFFER 5 ', type(self.varbuffer)
            #     print type(self.varbuffer[0]) #[j][midx,i,:nnu] = tempvar[i]

                    if self.only_flux:
                        ptemp[midx] += self.dpar[midx]
                        self.varbuffer[midx + 1][midx, 2, :nnu] = ptemp[midx]
        #       print 'VARBUFFER 6 ', type(self.varbuffer)
        #       print type(self.varbuffer[0]) #[j][midx,i,:nnu] = tempvar[i]

                    else:

                        for j in range(len(p)):
                            ptemp = [pi for pi in p]
                            ptemp[j] += self.dpar[j]
                            # Variables of current component
                            tempvar = currvar[midx](ptemp, self.freqs[spw])
                            if modi in self.isNumerical:
                                tempvar = self.gridModel(modi, tempvar)
                            for i in range(len(tempvar)):
                                self.varbuffer[j
                                               + 1][midx, i, :nnu] = tempvar[i]

            #        print 'VARBUFFER 7 ', type(self.varbuffer)
            #        print type(self.varbuffer[0]) #[j][midx,i,:nnu] = tempvar[i]

                else:  # Spectral mode

                    tempvar = currvar[midx](
                        p,
                        self.freqs[spw][nui])  # Variables of current component
                    #   print '\n\n',len(tempvar)
                    if modi in self.isNumerical:
                        tempvar = self.gridModel(modi, tempvar)
            #   print len(tempvar),'\n\n'
                    for i in range(len(tempvar)):
                        self.varbuffer[0][midx, i, 0] = tempvar[i]

    #       print 'VARBUFFER 8 ', type(self.varbuffer)
    #       print type(self.varbuffer[0]) #[j][midx,i,:nnu] = tempvar[i]

                    for j in range(len(p)):
                        ptemp = [pi for pi in p]
                        ptemp[j] += self.dpar[j]  #self.minnum
                        tempvar = currvar[midx](ptemp, self.freqs[spw][
                            nui])  # Variables of current component
                        #     print 'ANT0: ',len(tempvar)
                        if modi in self.isNumerical:
                            #       print 'ANT: ',len(tempvar)
                            tempvar = self.gridModel(modi, tempvar)
            #       print 'DESP: ',len(tempvar)
                        for i in range(len(tempvar)):
                            self.varbuffer[j + 1][midx, i, 0] = tempvar[i]

        #       print 'VARBUFFER 9 ', type(self.varbuffer)
        #       print type(self.varbuffer[0]) #[j][midx,i,:nnu] = tempvar[i]

#    print [self.varbuffer[0][mid,:6,0] for mid in range(len(currmod))]
#    print [self.varbuffer[1][mid,:6,0] for mid in range(len(currmod))]

            ChiSq += self.Ccompmodel(spw, nui, mode)

            if mode == -2:
                if nui == -1:
                    ndata += np.sum(
                        self.wgt[spw] >
                        0.0)  # Only add the unflagged data to compute the DoF
                else:
                    ndata += np.sum(
                        self.wgt[spw][:, nui] >
                        0.0)  #  Only add the unflagged data to compute the DoF

        if mode in [-2, -1]:
            if nui < 0:
                sys.stdout.write('\r Iteration # %i. ' % (self.calls))
                sys.stdout.write('\r \t\t\t\t Achieved ChiSq:  %.8e' % (ChiSq))
                sys.stdout.flush()

        if ChiSq <= 0.0:
            raise ValueError(
                'Invalid Chi Square! Maybe the current fitted value of flux (and/or size) is negative! Please, set BOUNDS to the fit!'
            )

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

    def ChiSquare(self, p, bounds, p_ini):
        """ Just a wrapper of the \'residuals\' function, to be called
        by simplex."""
        inside = True
        if bounds != None:
            for i, bound in enumerate(bounds):
                if type(bound) is list:
                    vmin = bound[0] != None and p[i] <= bound[0]
                    vmax = bound[1] != None and p[i] >= bound[1]
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
    """ Function to convert a string with spw selection into lists
      of channels to select/average. It follows the CASA syntax."""

    if spw == '':
        spw = ','.join(map(str, range(len(maxchans))))

    entries = spw.split(',')
    output = [[] for i in maxchans]

    for entry in entries:
        check = entry.split(':')

        spws = map(int, check[0].split('~'))
        if len(spws) == 1:
            selspw = [spws[0]]
        else:
            selspw = range(spws[0], spws[1] + 1)

        for sp in selspw:

            if sp + 1 > len(maxchans):
                errstr = 'There are only %i spw in the data!\n Please, revise the \'spw\' parameter' % (
                    len(maxchans))
                return [False, errstr]

            if len(check) == 1:
                chranges = ['0~' + str(maxchans[sp] - 1)]
            else:
                chans = check[1]
                chranges = chans.split(';')
            ranges = []
            for chran in chranges:
                ch1, ch2 = map(int, chran.split('~'))
                if ch1 > ch2:
                    errstr = '%i is larger than %i. Revise channels for spw %i' % (
                        ch1, ch2, sp)
                    return [False, errstr]
                ch2 = min([ch2, maxchans[sp] - 1])
                for i in range(ch1, ch2 + 1, width):
                    ranges.append(range(i, min([(i + width), ch2 + 1])))

            output[sp] = ranges

    return [True, output]


def modelFromClean(imname, ichan=0, echan=0):
    """Reads the \'*.model\' data created by task \'clean\'
     and returns lists compatible with the \'model\' and 
     \'var\' (or with the \'fixed\' and \'fixedvar\')
     entries of uvmultifit. 
     ichan and echan are the first and last channels to
     select. The flux densities that the task will 
     return for each pixel will be the frequency average
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
        cleans = np.transpose(
            [cleanlocs[0], cleanlocs[1]] + [modarray[cleanlocs]])
        totflux = np.sum(cleans[:, 2])

        # Put shifts in arcseconds w.r.t. image center:
        cleans[:, 0] = (
            cleans[:, 0] - refpix[0]) * deltaxy[0] * 180. / np.pi * 3600.
        cleans[:, 1] = (
            cleans[:, 1] - refpix[1]) * deltaxy[1] * 180. / np.pi * 3600.

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
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#a. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#b. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the
#   distribution.
#c. Neither the name of the author nor the names of contributors may
#   be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

from numpy import asfarray


def wrap_function(function, args):
    ncalls = [0]

    def function_wrapper(x):
        ncalls[0] += 1
        return function(x, *args)

    return ncalls, function_wrapper


def _mod_simplex(func,
                 x0,
                 args=(),
                 callback=None,
                 relstep=1.e-1,
                 relxtol=1.e-4,
                 relftol=1.e-3,
                 maxiter=None,
                 maxfev=None,
                 disp=False,
                 return_all=False,
                 **unknown_options):
    """
    Minimization of scalar function of one or more variables using the
    Nelder-Mead algorithm.

    Options for the Nelder-Mead algorithm are:
        disp : bool
            Set to True to print convergence messages.
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
    x0 = asfarray(x0).flatten()
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
        sim = np.zeros((N + 1, ), dtype=x0.dtype)
    else:
        sim = np.zeros((N + 1, N), dtype=x0.dtype)
    fsim = np.zeros((N + 1, ), float)
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
    # sort so sim[0,:] has the lowest function value
    sim = np.take(sim, ind, 0)

    iterations = 1

    minx = np.sqrt(np.finfo(np.float64).eps)
    if type(relxtol) is float:
        relxtol = relxtol * np.ones(len(x0))
    else:
        relxtol = np.array(relxtol)

    while (fcalls[0] < maxfun and iterations < maxiter):
        testsim = np.copy(sim[0])
        testsim[np.abs(testsim) < minx] = minx
        testfsim = max(fsim[0], minx)
        deltas = np.abs(sim[1:] - sim[0])
        smaller = len(deltas[deltas > relxtol])

        if (smaller == 0 and np.max(np.abs(
            (fsim[0] - fsim[1:]) / testfsim)) <= relftol):
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
                        print 'simj', sim[j]
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
    """ Similar to uvmultifit, but the engine uses images created by
       task clean. 

       All keywords are equal to those of uvmultifit, with the exception
       of these:

       psf = The image (or image cube) of the Point Spread Function.
             Usually, the name is in the form \'blablabla.psf\'


       residual = The image (or image cube) of the image residuals.
                  Usually, the name is in the form \'blablabla.residual\'

            NOTE: If you want primary-beam correction, you should use the 
                  'blablabla.image', instead.


       reinvert = If immultifit has already been run and the user wants
                  to refit the same data, the Fourier inversion of the 
                  images (or image cubes) can be avoided if reinvert=False.
                  Default: True


       dBcut = Mask the uv points with a signal lower than dBCut of the 
               peak (in dB). Use this to minimize convolution-like artifacts
               in UV space. Default: -30dB => SNR = 1000.

       The user should run clean with niter=0, to ensure that no 
       deconvolution is applied to the image of residuals, so it corresponds
       exactly to the inverse Fourier transform of the visibilities.

       Notice also that keywords "vis" and "spw", "timewidth" and "chanwidth"
       are NOT used by immultifit.

   """

    def __init__(self,
                 parent=None,
                 reinvert=True,
                 start=0,
                 nchan=-1,
                 psf='',
                 residual='',
                 dBcut=-30.,
                 **kwargs):
        self.psf = str(psf)
        self.residual = str(residual)
        self.start = int(start)
        self.nchan = int(nchan)
        self.reinvert = bool(reinvert)
        #  self.scipylm = True
        self.dBcut = dBcut
        uvmultifit.__init__(self, **kwargs)

    def writeModel(self):
        self.printInfo(
            '\n\n There is no measurement set (i.e., you are running \'immultifit\').\n Will generate an image, instead.\n'
        )
        if os.path.exists(self.residual + '.immultifit'):
            shutil.rmtree(self.residual + '.immultifit')
        os.system('cp -r %s %s' % (self.residual,
                                   self.residual + '.immultifit'))

        sq2 = np.sqrt(2.)

        # Set the correct total flux of model:
        q = (self.u[0]**2. + self.v[0]**2.)**0.5
        zerosp = q < 1.0
        for nui in range(self.start, self.nchan):
            self.mymodel.output[0][:, nui - self.start][zerosp] = np.max(
                self.mymodel.output[0][:, nui - self.start])

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
            self.printInfo('\r Doing frequency channel %i of %i' %
                           (i + 1 - self.start, self.nchan - self.start))
            for j in range(imdims[0]):
                for k in range(imdims[1]):
                    cpix = imdims[1] * j + k
                    temparr[j, k] = (
                        self.mymodel.output[0][cpix, i - self.start]
                    ) * self.averweights[0][cpix, i - self.start]
            #   temparr[j,k] -= 1.j*(self.mymodel.output[0][1][cpix,i-self.start])*self.averweights[0][cpix,i-self.start]
            resim[:, :, self.stokes, i] -= np.fft.fftshift(
                np.fft.ifft2(np.fft.ifftshift(temparr))).real / sq2

        ia.putchunk(resim)
        ia.setbrightnessunit('Jy/beam')
        ia.close()

        self.printInfo('\n\n Will now save the unconvolved model image.\n')
        modname = '.'.join(
            self.residual.split('.')[:-1]) + '.fitmodel.immultifit'
        if os.path.exists(modname):
            shutil.rmtree(modname)
        os.system('cp -r %s %s' % (self.residual, modname))

        ia.open(modname)
        resim = ia.getchunk()
        imdims = np.shape(resim)
        npix = float(imdims[0] * imdims[1])
        temparr = np.zeros((imdims[0], imdims[1]), dtype=np.complex128)
        for i in range(self.start, self.nchan):
            self.printInfo('\r Doing frequency channel %i of %i' %
                           (i + 1 - self.start, self.nchan - self.start))
            for j in range(imdims[0]):
                for k in range(imdims[1]):
                    cpix = imdims[1] * j + k
                    temparr[j, k] = self.mymodel.output[0][cpix,
                                                           i - self.start]
        #    temparr[j,k] -= 1.j*(self.mymodel.output[0][1][cpix,i-self.start])
            resim[:, :, self.stokes, i] = np.fft.fftshift(
                np.fft.ifft2(np.fft.ifftshift(temparr))).real / sq2

        ia.putchunk(resim)
        ia.setbrightnessunit('Jy/pixel')
        ia.close()

############################################
#
#  SANITY CHECKS AND PARAMETER SETTINGS
#

    def checkInputs(self, data_changed=False):
        self.printInfo(
            '\n\nIn image mode, the whole image is taken.\n Will override the \'spw\', \'field\', \'MJDrange\', and \'scan\' parameters.\n\n'
        )
        self.printInfo(
            'In image mode, the channelization used is that of the image.\n Will override the \'chanwidth\' and \'timewidth\' parameters.\n\n'
        )

        self.savemodel = True
        success = self.checkOrdinaryInputs()

        if not success:
            return False

        try:
            self.stokes = abs(int(self.stokes))
        except:
            self.printInfo(
                '\n\nIn image mode, \'stokes\' must be an integer\n(i.e., the Stokes column of the image).\n'
            )
            polimag = ['I', 'Q', 'U', 'V']
            if self.stokes in polimag:
                stchan = polimag.index(self.stokes)
                self.printInfo(
                    'Since stokes is \'%s\', will try with image column %i\n' %
                    (self.stokes, stchan))
                self.stokes = stchan
            else:
                self.printError(
                    '\nPlease, set the \'stokes\' keyword to an image pol. channel\n'
                )
                return False

        try:
            ia.open(self.residual)
            imdims = np.array(ia.shape())
            imsum = ia.summary()
            ia.close()
        except:
            self.printError(
                'The residual image is not an image (or does not exist)!\n\n')
            return False

        try:
            ia.open(self.psf)
            imdims2 = np.array(ia.shape())
            ia.close()
        except:
            self.printError(
                'The PSF image is not an image (or does not exist)!\n\n')
            return False

        if max(np.abs(imdims - imdims2)) > 0:
            self.printError(
                'The dimensions of the PSF image and the image of residuals do NOT coindice!\n\n'
            )
            return False

        if imdims[2] > self.stokes:
            self.printInfo('Selecting stokes column %i' % self.stokes)
        else:
            self.printError(
                '\n\nThe images only have %i stokes columns, but you selected stokes=%i\n'
                % (imdims[2], self.stokes))
            return False

        self.MJDrange = [0., 0.]
        self.spwlist = [[0, 0, [range(imdims[3])], 0]]

        RAd = np.where(imsum['axisnames'] == 'Right Ascension')[0][0]
        Decd = np.where(imsum['axisnames'] == 'Declination')[0][0]
        self.refpos = np.array([
            imsum['refval'][RAd], imsum['refval'][Decd]
        ])  #np.copy(csys.torecord()['direction0']['crval'])

        if (self.nchan < 0): self.nchan = imdims[3]

        self.nchan += self.start

        if self.start >= imdims[3]:
            self.printError(
                'Starting channel (%i) must be smaller than number of channels (%i)'
                % (self.start, imdims[3]))
        if self.nchan > imdims[3]:
            self.printError(
                'Ending channel (%i) must be smaller than number of channels (%i)'
                % (self.nchan, imdims[3]))

# Try to compile the equations for the variables:
        if data_changed:
            try:
                del self.mymodel
            except:
                pass

            self.mymodel = modeler(
                self.model, self.var, self.fixed, self.fixedvar, self.scalefix,
                self.NCPU, self.only_flux, self.applyHankel, self.isNumerical,
                self.useGains, [self.phase_gains, self.amp_gains])
        else:
            self.mymodel.var = self.var
            self.mymodel.model = self.model
            self.mymodel.fixed = self.fixed
            self.mymodel.fixedvar = self.fixedvar
            self.mymodel.scalefix = self.scalefix
            self.mymodel.calls = 0
            self.mymodel.removeFixed = False

        if self.mymodel.failed:
            self.printError(self.mymodel.resultstring)
            return False

        if self.proper_motion == 0.0:
            self.proper_motion = [[0., 0.] for i in self.model]
        elif type(self.proper_motion) is not list:
            self.printError('\'proper_motion\' must be a list!')

        self.setEngineWgt()
        return True

############################################
#
#  NO DDE WEIGHTING
#

    def setWgtEq(self):
        self.mymodel.KfacWgt = 0.0
        return True


############################################
#
#  READ DATA (IMAGES AND PSFs)
#

    def readData(self, data_changed=True):

        if self.reinvert:

            os.system('rm -rf %s.fft.*' % self.residual)
            os.system('rm -rf %s.fft.*' % self.psf)

            self.printInfo('\nInverting images into Fourier space\n')
            ia.open(self.residual)
            ia.fft(
                real=self.residual + '.fft.real',
                imag=self.residual + '.fft.imag')
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
            freqs[i] = ia.toworld([0, 0, self.stokes,
                                   i + self.start])['numeric'][3]

        self.printInfo('\nReading gridded UV coordinates and weights\n')
        u0, v0, s0, nu0 = ia.toworld([0, 0, self.stokes,
                                      self.start])['numeric']
        u1, v1, s0, nu1 = ia.toworld([1, 1, self.stokes,
                                      self.start])['numeric']
        du = u1 - u0
        dv = v1 - v0
        for j in range(imdims[0]):
            self.printInfo('\r Reading row %i of %i' % (j, imdims[0]))
            for k in range(imdims[1]):
                cpix = imdims[1] * j + k
                ui[cpix] = (u0 + du * j)
                vi[cpix] = (v0 + dv * k)
                wgti[cpix, :] = np.copy(wgt[j, k, self.stokes, :])

        del wgt
        ia.close()
        maskwgt = wgti < np.max(wgti) * (10.**(
            self.dBcut / 10.))  # Maximum dynamic range set by "Window effect"
        self.printInfo('\n\nReading gridded visibilities\n')
        ia.open(self.residual + '.fft.real')
        datas = ia.getchunk()
        for j in range(imdims[0]):
            self.printInfo('\r Reading row %i of %i' % (j, imdims[0]))
            for k in range(imdims[1]):
                cpix = imdims[1] * j + k
                datare[cpix, :] = datas[j, k, self.stokes, self.start:
                                        self.nchan]

        del datas
        ia.close()

        self.printInfo(' ')
        ia.open(self.residual + '.fft.imag')
        datas = ia.getchunk()
        for j in range(imdims[0]):
            self.printInfo('\r Reading row %i of %i' % (j, imdims[0]))
            for k in range(imdims[1]):
                cpix = imdims[1] * j + k
                dataim[
                    cpix, :] = -datas[j, k, self.stokes, self.start:self.nchan]

        del datas
        ia.close()

        toLamb = self.LtSpeed / freqs[0]
        self.u = [ui * toLamb]
        self.v = [vi * toLamb]

        self.averdata = [(datare + 1.j * dataim) / wgti]
        bads = np.logical_or(
            np.isnan(self.averdata[0].real), np.isnan(self.averdata[0].imag))
        badmsk = np.logical_or(bads, maskwgt)
        self.averdata[0][badmsk] = 0.0
        #     self.averdata[0][1][badmsk]=0.0
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
        self.tArr = [np.linspace(0., 1., npix2)]
        self.tIdx = [np.arange(npix2, dtype=np.int32)]

        self.Nants = 2

        self.printInfo('\nDone reading.\n')

        return True
