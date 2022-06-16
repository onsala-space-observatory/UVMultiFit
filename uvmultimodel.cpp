/*
#
# UVMULTIFIT - C++ MULTI-THREADED CORE ENGINE.
#
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
*/
#include <Python.h>
#include <numpy/arrayobject.h>
#include <pthread.h>
#include <gsl/gsl_sf_bessel.h>
#include <stdio.h>
#include <sys/types.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

#if QUINN_FITTER == 0
#include "QuinnFringe.h"
#endif

#define DEBUG
void show_info(const char *var, PyObject *obj)
{
    if (obj == NULL) return;

#ifndef DEBUG
    return;
#endif

    // NPY_BOOL: 0
    // NPY_INT16: 3
    // NPY_INT32: 5
    // NPY_INT64: 7
    // NPY_UINT16: 4
    // NPY_UINT32: 6
    // NPY_UINT64: 8
    // NPY_FLOAT32: 11
    // NPY_FLOAT64: 12
    // NPY_COMPLEX64: 14
    // NPY_COMPLEX128: 15

    // std::cout << "NPY_BOOL: " << NPY_BOOL << std::endl;
    // std::cout << "NPY_INT: " << NPY_INT << std::endl;
    // std::cout << "NPY_INT16: " << NPY_INT16 << std::endl;
    // std::cout << "NPY_INT32: " << NPY_INT32 << std::endl;
    // std::cout << "NPY_INT64: " << NPY_INT64 << std::endl;
    // std::cout << "NPY_UINT16: " << NPY_UINT16 << std::endl;
    // std::cout << "NPY_UINT32: " << NPY_UINT32 << std::endl;
    // std::cout << "NPY_UINT64: " << NPY_UINT64 << std::endl;
    // std::cout << "NPY_FLOAT32: " << NPY_FLOAT32 << std::endl;
    // std::cout << "NPY_FLOAT64: " << NPY_FLOAT64 << std::endl;
    // std::cout << "NPY_COMPLEX64: " << NPY_COMPLEX64 << std::endl;
    // std::cout << "NPY_COMPLEX128: " << NPY_COMPLEX128 << std::endl;

    int ndim = PyArray_NDIM((PyArrayObject *)obj);
    int type = PyArray_TYPE((PyArrayObject *)obj);
    int size = PyArray_SIZE((PyArrayObject *)obj);

    std::cout << "'" << std::setw(10) << var
              << "' " << "is an array (" << PyArray_CheckExact(obj) << "):"
              << " of type " << std::setw(2) << type
              << " and dimension " << ndim;
    if (ndim > 1) {
        std::cout << " [";
        for (int i = 0; i < ndim; i++) {
            std::cout << " " << PyArray_DIM((PyArrayObject *)obj, i)
                      << (i == ndim-1 ? " ]" : " x");
        }
    }

    std::cout << ", size " << size << std::endl;
}

/* Docstrings */
static char module_docstring[] =
    "This module provides an interface for least-square visibility fitting.";
static char uvmultimodel_docstring[] =
    "Calculate the residuals and chi square of a multi-component model";
static char clearPointers_docstring[] =
    "Delete the data pointers.";
static char setData_docstring[] =
    "Get the data pointers.";
static char setNspw_docstring[] =
    "Set up the pointers to the data arrays.";
static char setModel_docstring[] =
    "Set up the model components and variable parameters.";
static char setNCPU_docstring[] =
    "Set up the parallelization.";
static char setWork_docstring[] =
    "Allocate memory to compute Hessian and error vector.";
static char QuinnFF_docstring[] =
    "Perform Fringe Fitting, based on the delay-rate fringe peaks, using the Quinn estimator for the peak.";

/* Available functions */
static PyObject *setNspw(PyObject *self, PyObject *args);
static PyObject *setData(PyObject *self, PyObject *args);
static PyObject *clearPointers(PyObject *self, PyObject *args);
static PyObject *setNCPU(PyObject *self, PyObject *args);
static PyObject *modelcomp(PyObject *self, PyObject *args);
static PyObject *setModel(PyObject *self, PyObject *args);
static PyObject *setWork(PyObject *self, PyObject *args);
static PyObject *QuinnFF(PyObject *self, PyObject *args);

void *writemod(void *work);

/* Module specification */
static PyMethodDef module_methods[] = {
    {"setData", setData, METH_VARARGS, setData_docstring},
    {"setNspw", setNspw, METH_VARARGS, setNspw_docstring},
    {"setModel", setModel, METH_VARARGS, setModel_docstring},
    {"setNCPU", setNCPU, METH_VARARGS, setNCPU_docstring},
    {"modelcomp", modelcomp, METH_VARARGS, uvmultimodel_docstring},
    {"setWork", setWork, METH_VARARGS, setWork_docstring},
    {"QuinnFF", QuinnFF, METH_VARARGS, QuinnFF_docstring},
    {"clearPointers", clearPointers, METH_VARARGS, clearPointers_docstring},
    {NULL, NULL, 0, NULL}
};

typedef std::complex<double> cplx64;

static int Nspw, NCPU, nui, cIF, ncomp, npar, maxnchan=0, HankelOrder=80;
static int mode, Nants;
static double cosDecRef, sinDecRef;
static bool compFixed, isModel, MixedG;
static int NparMax = 7; // Degrees of freedom for components (i.e., RA, Dec, Flux, Major, Ratio, PA, Inner).

// ALL POINTERS TO DATA AND METADATA:
struct DATA {
    std::vector<int *> ants[2];
    std::vector<int *> dtIndex;
    std::vector<double *> freqs;
    std::vector<double *> uv[3];
    std::vector<double *> wgt[2];
    std::vector<double *> dt;
    std::vector<double *> dtArray;
    std::vector<double *> RAshift;
    std::vector<double *> Decshift;
    std::vector<double *> Stretch;
    std::vector<cplx64 *> ObsVis;
    std::vector<int8_t *> fittable;
    std::vector<int8_t *> isGain;
    std::vector<int> nnu;
    std::vector<int> nt;
    double *phaseCenter;
};

// ALL POINTERS TO MODEL-RELATED STUFF:
struct MODEL {
    std::vector<cplx64 *> ModVis;
    std::vector<std::vector<std::vector<cplx64 *>>> Gain;
    std::vector<double *> vars;
    std::vector<double *> fixp;
    std::vector<double> Chi2;
    std::vector<std::vector<double>> WorkHess;
    std::vector<std::vector<double>> WorkGrad;
    std::vector<int *> parAnt;
    std::vector<int> nparAnt;
    int *models;
    double *Hessian;
    double *Gradient;
    double *dpar;
    double *muRA;
    double *muDec;
};

/* Structure to pass, as void cast, to the workers */
struct SHARED_DATA {
    std::vector<int> t0;
    std::vector<int> t1;
    int Iam;
};

static SHARED_DATA master;
std::vector<SHARED_DATA> worker;

static DATA vis;
static MODEL mod;

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "uvmultimodel",
    module_docstring,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

/* Initialize the module */
PyMODINIT_FUNC PyInit_uvmultimodel(void)
{
    // PyObject *m = Py_InitModule3("uvmultimodel", module_methods, module_docstring);
    PyObject *m = PyModule_Create(&moduledef);
    if (m == NULL) return NULL;

    // Initiate variables with dummy values:
    NCPU = 0;
    Nspw = 0;
    npar = -1;
    Nants = 0;
    ncomp = -1;
    master.t0.resize(1);
    master.t1.resize(1);

    /* Load `numpy` functionality. */
    import_array();
    return m;
}

static PyObject *clearPointers(PyObject *self, PyObject *args)
{
    int i;
    if (!PyArg_ParseTuple(args, "i", &i)) {
        printf("FAILED clearPointers!\n");
        return NULL;
    }
    printf("inside uvmod.clearPointers: %i\n", i);

    switch (i) {
      case 0:
        break;
      case 1:
        break;
      case 2:
        // NCPU = -1;
        // Nspw = -1;
        // npar = -1;
        // Nants = -1;
        // ncomp = -1;
        master.t0.resize(1);
        master.t1.resize(1);
        isModel = false;
    }

    return Py_BuildValue("i", i);
}

// Taylor expansion of the Hankel transform, for special models (GaussianRing):
void GridModel(int imod, double UVRad, double *varp, int *currow, double &tempAmp)
{
    double UVRadP = 1.0;
    double UVRadSq = UVRad*UVRad/4.;
    int i;

    if (imod == 9) {
        tempAmp = 0.0;
        for (i = 0; i < HankelOrder; i++) {
            tempAmp += varp[currow[NparMax+i]]*UVRadP;
            UVRadP *= UVRadSq;
        }
    }
}

/* Code for the workers (NOT CALLABLE FROM PYTHON). */
void *writemod(void *work)
{
    SHARED_DATA *mydata = (SHARED_DATA *)work;

    int nu0, nu1;
    if (nui == -1) {
        nu0 = 0;
        nu1 = vis.nnu[cIF];
    } else {
        nu0 = nui;
        nu1 = nui+1;
    }

    int k, i, t, j, m, p, currow[NparMax+HankelOrder];

    double phase, uu, vv, ww, UVRad, Ampli, Ampli0, DerivRe[npar], DerivIm[npar];
    double SPA, CPA, tA, tB, tempChi, deltat, Ellip, Ellip0;
    Ellip0 = 1.;
    Ellip = 1.;
    Ampli0 = 1.;
    double wgtcorrA, tempRe[npar+1], tempIm[npar+1], tempAmp;
    tempChi = 0.0;

    const double deg2rad = 0.017453292519943295;
    const double sec2rad = 3.1415926535/180./3600.;
    const double radsec2 = pow(3.1415926535/180./3600.,2.); // 2.3504430539097885e-11;
    int Iam = mydata->Iam;
    int widx = (Iam == -1) ? 0 : Iam;
    // int pmax = (mode == -1) ? npar+1 : 1;
    // int mmax = (mode == 0) ? 1 : ncomp;
    int pmax = 0, mmax = 0;
    bool write2mod = false, writeDer = false;

    switch (mode) {
      case  0:
        mmax = 1;
        write2mod = true;         // COMPUTE FIXED MODEL
        pmax = 1;
        writeDer  = false;
        break;
      case -1:
        mmax = ncomp;
        write2mod = false;        // HESS & ERRORS & CHISQ
        pmax = npar+1;
        writeDer = true;
        break;
      case -2:
        mmax = ncomp;
        write2mod = false;        // GET CHISQ
        pmax = 1;
        writeDer  = false;
        break;
      case -3:
        mmax = ncomp;
        write2mod = true;         // ADD VARMOD TO FIXED
        pmax = 1;
        writeDer  = false;
        break;
      case -4:
        mmax = ncomp;
        write2mod = true;         // ADD VARMOD TO FIXED
        pmax = 1;
        writeDer  = false;
        break;
      case -5:
        mmax = 0;
        write2mod = false;        // CALIBRATE
        pmax = 1;
        writeDer  = false;
        break;
    }

    // bool write2mod = (mode == -3 || mode == -4 || mode == -5 || mode >= 0);
    // bool writeDer = (mode==-1);
    bool EllipChanged;
    double tempD, tempR, tempI, cosphase, sinphase, cosphase0, sinphase0, PA;
    double wterm, rsh, dsh, wamp, tempres0, tempres1, ll, mm, PBcorr;
    std::vector<cplx64> totGain(pmax);

    int ant1, ant2, pdep, currTIdx, kT;
    cplx64 GCplx;
    bool calibrate, PBlimit;

    for (t = mydata->t0[cIF]; t < mydata->t1[cIF]; t++) {
        deltat = vis.dt[cIF][t];

        currTIdx = vis.dtIndex[cIF][t];
        ant1 = vis.ants[0][cIF][t];
        ant2 = vis.ants[1][cIF][t];

        if (vis.fittable[cIF][t]!=0) {
            for (i = nu0; i < nu1; i++) {
                j = (nui!=-1) ? 0 : i;
                k = vis.nnu[cIF]*t+i;

                kT = MixedG ? vis.nnu[cIF]*currTIdx+i : vis.nnu[cIF];
                tempD = vis.wgt[0][cIF][k]*vis.wgt[0][cIF][k];

                if (vis.isGain[cIF][t] != 0) {
                    for (p = 0; p < pmax; p++) {
                        calibrate = false;
                        GCplx = cplx64(1.0);

                        for (pdep = 0; pdep < mod.nparAnt[ant1]; pdep++) {
                            if (mod.parAnt[ant1][pdep] == p) {
                                if (MixedG) {
                                    GCplx *= mod.Gain[cIF][ant1][pdep][kT];
                                } else {
                                    GCplx *= mod.Gain[cIF][ant1][pdep][i]*mod.Gain[cIF][ant1][pdep][kT+currTIdx];
                                }
                                calibrate=true;
                            }
                        }

                        if (!calibrate) {
                            if (MixedG) {
                                GCplx *= mod.Gain[cIF][ant1][0][kT];
                            } else {
                                GCplx *= mod.Gain[cIF][ant1][0][i]*mod.Gain[cIF][ant1][0][kT+currTIdx];
                            }
                        }

                        calibrate = false;
                        for (pdep = 0; pdep < mod.nparAnt[ant2]; pdep++) {
                            if (mod.parAnt[ant2][pdep] == p) {
                                if (MixedG) {
                                    GCplx *= std::conj(mod.Gain[cIF][ant2][pdep][kT]);
                                } else {
                                    GCplx *= std::conj(mod.Gain[cIF][ant2][pdep][i]*mod.Gain[cIF][ant2][pdep][kT+currTIdx]);
                                }
                                calibrate=true;
                            }
                        }

                        if (!calibrate) {
                            if (MixedG) {
                                GCplx *= std::conj(mod.Gain[cIF][ant2][0][kT]);
                            } else {
                                GCplx *= std::conj(mod.Gain[cIF][ant2][0][i]*mod.Gain[cIF][ant2][0][kT+currTIdx]);
                            }
                        }

                        totGain[p] = GCplx;
                    }
                }

                uu = vis.freqs[cIF][i]*vis.uv[0][cIF][t];
                vv = vis.freqs[cIF][i]*vis.uv[1][cIF][t];
                ww = vis.freqs[cIF][i]*vis.uv[2][cIF][t];

                for (p = 0; p < pmax; p++) {
                    tempRe[p] = 0.0;
                    tempIm[p] = 0.0;
                }

                for (m = 0; m < mmax; m++) {
                    PBlimit = false;

                    for (p = 0; p < NparMax+HankelOrder; p++) {
                        currow[p] = m*maxnchan*(NparMax+HankelOrder)+j+p*maxnchan;
                    }

                    for (p = 0; p < pmax; p++) {
                        if (!PBlimit && (p == 0 || (mod.vars[p][currow[0]] != mod.vars[0][currow[0]] || mod.vars[p][currow[1]] != mod.vars[0][currow[1]]))) {

                            // Project RA:
                            rsh = (mod.vars[p][currow[0]] + mod.muRA[m]*deltat)/cosDecRef - vis.RAshift[cIF][t];

                            // Project Dec:
                            dsh = mod.vars[p][currow[1]] - vis.Decshift[cIF][t] + mod.muDec[m]*deltat;

                            // Get source-centered shifts (l,m):
                            tempR = vis.Decshift[cIF][t]*sec2rad+vis.phaseCenter[1];
                            tempI = tempR + dsh*sec2rad;
                            mm = asin(cos(tempR)*sin(tempI)-sin(tempR)*cos(tempI)*cos(rsh*sec2rad))/sec2rad;
                            ll = atan(cos(tempI)*sin(rsh*sec2rad)/(cos(tempR)*cos(tempI)*cos(rsh*sec2rad)+sin(tempR)*sin(tempI)))/sec2rad;
                            PBcorr = vis.wgt[1][cIF][t]*(mm*mm + ll*ll)*vis.freqs[cIF][i]*vis.freqs[cIF][i];
                            // PBcorr = vis.wgt[1][cIF][m*vis.nt[cIF]+t]*(mm*mm + ll*ll)*vis.freqs[cIF][i]*vis.freqs[cIF][i];

                            // if (p==0&&t==10) {printf("PBCorr: %.3e m=%i\n",PBcorr,m);}

                            wgtcorrA = exp(PBcorr);
                            PBlimit = false;
                            // ACTIVATE THIS TO AVOID OVER-COMPUTING MODEL VISIBILITIES IN VERY LARGE MOSAICS (i.e. 3-SIGMA PBEAM CUTOFF):
                            // if (PBcorr<-3.0) {wgtcorrA = exp(PBcorr); PBlimit = false;} else {wgtcorrA = 0.0; PBlimit=true;}
                            phase = ll*uu + mm*vv;
                            wamp = sqrt(1. - (ll*ll + mm*mm)*radsec2);
                            wterm = ww*(wamp - 1.)*sec2rad;
                            cosphase = cos(phase+wterm);
                            sinphase = sin(phase+wterm);
                            if (p == 0) {
                                cosphase0 = cosphase;
                                sinphase0 = sinphase;
                            }
                        } else {
                            cosphase = cosphase0;
                            sinphase = sinphase0;
                            PBlimit = false;
                        }

                        if (!PBlimit && (p == 0 || (mod.vars[p][currow[5]] != mod.vars[0][currow[5]] || mod.vars[p][currow[4]] != mod.vars[0][currow[4]])) ) {
                            EllipChanged = true;
                            if (mod.models[m] != 0) {
                                PA = mod.vars[p][currow[5]]*deg2rad;
                                SPA = sin(PA);
                                CPA = cos(PA);

                                tA = (uu*CPA - vv*SPA)*mod.vars[p][currow[4]];
                                tB = (uu*SPA + vv*CPA);
                                Ellip = sqrt(tA*tA+tB*tB);
                            } else {
                                Ellip = 1.0;
                            }
                            if (p == 0) {
                                Ellip0 = Ellip;
                            }
                        } else {
                            Ellip = Ellip0;
                            EllipChanged = false;
                        }

                        if ( !PBlimit && (p == 0 || EllipChanged || (mod.vars[p][currow[3]] != mod.vars[0][currow[3]] || mod.vars[p][currow[6]] != mod.vars[0][currow[6]]) )) {
                            if (mod.models[m] != 0) {
                                UVRad = Ellip*(mod.vars[p][currow[3]]/2.0);
                                tempAmp = 1.0;
                                if (mod.models[m] > 8) {
                                    GridModel(mod.models[m], UVRad, mod.vars[p], currow, tempAmp);
                                }

                                if (UVRad > 0.0) {
                                    switch (mod.models[m]) {
                                      case 1:
                                        Ampli = exp(-0.3606737602*UVRad*UVRad);
                                        break;
                                      case 2:
                                        Ampli = 2.0*gsl_sf_bessel_J1(UVRad)/UVRad;
                                        break;
                                      case 3:
                                        Ampli = gsl_sf_bessel_J0(UVRad);
                                        break;
                                      case 4:
                                        Ampli = 3.0*(sin(UVRad)-UVRad*cos(UVRad))/(UVRad*UVRad*UVRad);
                                        break;
                                      case 5:
                                        Ampli = sin(UVRad)/UVRad;
                                        break;
                                      case 6:
                                        Ampli = pow(1.+2.0813689810056077*UVRad*UVRad,-1.5);
                                        break;
                                      case 7:
                                        Ampli = 0.459224094*gsl_sf_bessel_K0(UVRad);
                                        break;
                                      case 8:
                                        Ampli = exp(-UVRad*1.3047660265);
                                        break;
                                    default:
                                        Ampli = tempAmp;
                                    }
                                } else {
                                    vis.wgt[0][cIF][k] = 0.0;
                                    Ampli=1.0;
                                }

                            } else {
                                Ampli = 1.0;
                            }

                            Ampli *= wgtcorrA;
                            if (p == 0) {
                                Ampli0 = Ampli;
                            }

                        } else if (!PBlimit) {
                            Ampli = Ampli0;
                        } else {
                            Ampli = 0.0;
                        }

                        if (!PBlimit) {
                            if (vis.isGain[cIF][t] != 0) {
                                tempR = mod.vars[p][currow[2]]*Ampli*cosphase;
                                tempI = mod.vars[p][currow[2]]*Ampli*sinphase;
                                tempRe[p] += totGain[p].real()*tempR - totGain[p].imag()*tempI;
                                tempIm[p] += totGain[p].imag()*tempR + totGain[p].real()*tempI;
                            } else {
                                tempRe[p] += mod.vars[p][currow[2]]*Ampli*cosphase;
                                tempIm[p] += mod.vars[p][currow[2]]*Ampli*sinphase;
                            }
                        }
                    }
                }

                if (compFixed && (mode != -5)) { // Add from model array (scaling with flux)
                    if (vis.isGain[cIF][t]!=0) {
                        for (p = 0; p < pmax; p++) {
                            if (vis.isGain[cIF][t]!=0) {
                                tempR = mod.ModVis[cIF][k].real()*mod.fixp[p][j];
                                tempI = mod.ModVis[cIF][k].imag()*mod.fixp[p][j];
                                tempRe[p] += totGain[p].real()*tempR - totGain[p].imag()*tempI;
                                tempIm[p] += totGain[p].imag()*tempR + totGain[p].real()*tempI;
                            } else {
                                tempRe[p] += mod.ModVis[cIF][k].real()*mod.fixp[p][j];
                                tempIm[p] += mod.ModVis[cIF][k].imag()*mod.fixp[p][j];
                            }
                        }
                    }
                }

                if (write2mod) { // Add to model array:
                    mod.ModVis[cIF][k] = cplx64(tempRe[0],tempIm[0]);
                }

                // Save residuals instead:
                if (mode == -4) {
                    mod.ModVis[cIF][k] = vis.ObsVis[cIF][k] - mod.ModVis[cIF][k];
                } else if (mode == -5) {
                    mod.ModVis[cIF][k] = vis.ObsVis[cIF][k]/totGain[0];
                }

                tempres0 = (vis.ObsVis[cIF][k].real()-tempRe[0])*vis.wgt[0][cIF][k];
                tempChi += tempres0*tempres0;
                tempres1 = (vis.ObsVis[cIF][k].imag()-tempIm[0])*vis.wgt[0][cIF][k];
                tempChi += tempres1*tempres1;

                if (writeDer) {
                    for (p = 0; p < npar; p++) {
                        DerivRe[p] = (tempRe[p+1]-tempRe[0])/mod.dpar[p];
                        DerivIm[p] = (tempIm[p+1]-tempIm[0])/mod.dpar[p];
                        mod.WorkGrad[widx][p] += (vis.ObsVis[cIF][k].real()-tempRe[0])*tempD*DerivRe[p];
                        mod.WorkGrad[widx][p] += (vis.ObsVis[cIF][k].imag()-tempIm[0])*tempD*DerivIm[p];
                    }

                    for (p = 0; p < npar; p++) {
                        for (m = p; m < npar; m++) {
                            mod.WorkHess[widx][npar*p+m] += tempD*(DerivRe[p]*DerivRe[m] + DerivIm[p]*DerivIm[m]);
                        }
                    }
                }
            }
        }
    }

    if (writeDer) {
        for (p = 0; p < npar; p++) {
            for (m = p; m < npar; m++) {
                mod.WorkHess[widx][npar*m+p] = mod.WorkHess[widx][npar*p+m];
            }
        }
    }

    if (Iam != -1) {
        mod.Chi2[Iam] = tempChi;
        pthread_exit((void*) 0);
    } else {
        mod.Chi2[0] = tempChi;
        return (void*) 0;
    }
}
// END OF CODE FOR THE WORKERS

// SET THE NUMBER OF IFS (IT REINITIATES ALL SHARED DATA VARIABLES!):
// USAGE FROM PYTHON: setNspw(i) where i is the number of SPW
static PyObject *setNspw(PyObject *self, PyObject *args)
{
    int i = 0;
    if (!PyArg_ParseTuple(args, "i",&i)) {
        printf("FAILED setNspw!\n");
        fflush(stdout);
        return NULL;
    }
    printf("inside uvmod.setNspw: %i\n", i);

    vis.nnu.resize(i);
    vis.nt.resize(i);
    master.t0.resize(i);
    master.t1.resize(i);

    vis.freqs.resize(i);

    vis.ants[0].resize(i);
    vis.ants[1].resize(i);

    vis.dtIndex.resize(i);
    vis.dtArray.resize(i);

    vis.uv[0].resize(i);
    vis.uv[1].resize(i);
    vis.uv[2].resize(i);
    vis.ObsVis.resize(i);
    mod.ModVis.resize(i);
    vis.fittable.resize(i);
    vis.isGain.resize(i);
    vis.wgt[1].resize(i);
    vis.wgt[0].resize(i);
    vis.dt.resize(i);
    vis.RAshift.resize(i);
    vis.Decshift.resize(i);
    vis.Stretch.resize(i);
    mod.Gain.resize(i);

    Nspw = i;
    // printf("successfully set Nspw to %i\n", Nspw);

    PyObject *ret = Py_BuildValue("i", Nspw);
    return ret;
}

// PREPARE PARALLELIZATION (MUST BE RUN AFTER setData)
// USAGE FROM PYTHON: setNCPU(i) where i is the num. of threads allowed
static PyObject *setNCPU(PyObject *self, PyObject *args)
{
    int i, j, k, k0, k1;
    if (!PyArg_ParseTuple(args, "i", &i)) {
        printf("FAILED setNCPU!\n");
        fflush(stdout);
        return NULL;
    }
    printf("inside uvmod.setNCPU: %i\n",i);

    // printf("Preparing memory for %i workers\n", i);
    // for (j = 0; j < NCPU; j++) {
    // delete[] worker[j].t0;
    // delete[] worker[j].t1;
    // }

    // if (NCPU > 0) {
    //     // delete[] worker;
    //     // delete[] mod.Chi2;
    //     // delete[] mod.WorkHess;
    //     // delete[] mod.WorkGrad;
    // }

    worker.resize(i);
    // printf("Preparing workers for %i spectral windows\n", Nspw);
    for (k = 0; k < i; k++) {
        if (Nspw <= 0) Nspw = 1;
        worker[k].t0.resize(Nspw);
        worker[k].t1.resize(Nspw);
        worker[k].Iam = k;
    }

    for (j = 0; j < Nspw; j++) {
        int nperproc = (int)((double)vis.nt[j])/((double)i);
        for (k = 0; k < i; k++) {
            k0 = nperproc*k;
            if (k == i-1) {
                k1 = vis.nt[j];
            } else {
                k1 = nperproc*(k+1);
            }
            worker[k].t0[j] = k0;
            worker[k].t1[j] = k1;
        }
    }

    mod.WorkHess.resize(i);
    mod.WorkGrad.resize(i);
    mod.Chi2.resize(i);
    // for (int j = 0; j < i; j++) {
    //     mod.WorkHess[j] = NULL;
    //     mod.WorkGrad[j] = NULL;
    // }
    NCPU = i;

    PyObject *ret = Py_BuildValue("i", NCPU);
    return ret;
}

// Fill-in the DATA arrays and master SHARED_DATA object
// USAGE FROM PYTHON: setData(IF, arrlist) where arrlist is list of data arrays.
//                    and IF is the IF number (setNspw must be run first!)
static PyObject *setData(PyObject *self, PyObject *args)
{
    // PyObject *obj;
    PyObject *pu, *pv, *pw, *pwgt, *preal, *poreal, *tArr, *tIdx;
    PyObject *pfreqs, *pfittable, *ant1l, *ant2l;
    PyObject *pwgtcorr, *dtime, *RAoffset, *Decoffset, *Stretchoff, *iG;
    int IF;

    if (!PyArg_ParseTuple(args, "iOOOOOOOOOOOOOOOOOOi", &IF, &pu, &pv, &pw, &pwgt, &preal, &poreal,
                          &pfreqs, &pfittable, &pwgtcorr, &dtime, &tArr, &tIdx, &RAoffset, &Decoffset,
                          &Stretchoff, &ant1l, &ant2l, &iG, &Nants)) {
        printf("FAILED setData!\n");
        fflush(stdout);
        return NULL;
    }

    printf("inside uvmod.setData: IF: %i, Ants: %i\n", IF, Nants);

    show_info("pu", pu);
    show_info("pv", pv);
    show_info("pw", pw);

    vis.uv[0][IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(pu, NPY_DOUBLE, NPY_IN_ARRAY));
    vis.uv[1][IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(pv, NPY_DOUBLE, NPY_IN_ARRAY));
    vis.uv[2][IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(pw, NPY_DOUBLE, NPY_IN_ARRAY));

    show_info("pwgt", pwgt);
    show_info("pwgtcorr", pwgtcorr);
    vis.wgt[0][IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(pwgt, NPY_DOUBLE, NPY_IN_ARRAY));
    vis.wgt[1][IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(pwgtcorr, NPY_DOUBLE, NPY_IN_ARRAY));

    show_info("preal", preal);
    show_info("poreal", poreal);
    vis.ObsVis[IF] = (cplx64 *)PyArray_DATA(PyArray_FROM_OTF(preal, NPY_COMPLEX128, NPY_IN_ARRAY));
    mod.ModVis[IF] = (cplx64 *)PyArray_DATA(PyArray_FROM_OTF(poreal, NPY_COMPLEX128, NPY_IN_ARRAY));

    show_info("pfreqs", pfreqs);
    vis.freqs[IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(pfreqs, NPY_DOUBLE, NPY_IN_ARRAY));

    show_info("pfittable", pfittable);
    vis.fittable[IF] = (int8_t *)PyArray_DATA(PyArray_FROM_OTF(pfittable, NPY_INT8, NPY_IN_ARRAY));

    show_info("dtime", dtime);
    vis.dt[IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(dtime, NPY_DOUBLE, NPY_IN_ARRAY));

    show_info("tArr", tArr);
    vis.dtArray[IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(tArr, NPY_DOUBLE, NPY_IN_ARRAY));

    show_info("tIdx", tIdx);
    vis.dtIndex[IF] = (int *)PyArray_DATA(PyArray_FROM_OTF(tIdx, NPY_INT, NPY_IN_ARRAY));

    show_info("RAoffset", RAoffset);
    show_info("Decoffset", Decoffset);
    show_info("Stretchoff", Stretchoff);
    vis.RAshift[IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(RAoffset, NPY_DOUBLE, NPY_IN_ARRAY));
    vis.Decshift[IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(Decoffset, NPY_DOUBLE, NPY_IN_ARRAY));
    vis.Stretch[IF] = (double *)PyArray_DATA(PyArray_FROM_OTF(Stretchoff, NPY_DOUBLE, NPY_IN_ARRAY));

    show_info("ant1l", ant1l);
    show_info("ant2l", ant2l);
    vis.ants[0][IF] = (int *)PyArray_DATA(PyArray_FROM_OTF(ant1l, NPY_INT, NPY_IN_ARRAY));
    vis.ants[1][IF] = (int *)PyArray_DATA(PyArray_FROM_OTF(ant2l, NPY_INT, NPY_IN_ARRAY));

    show_info("iG", iG);
    vis.isGain[IF] = (int8_t *)PyArray_DATA(PyArray_FROM_OTF(iG, NPY_INT8, NPY_IN_ARRAY));

    vis.nt[IF] = PyArray_DIM(preal, 0);
    vis.nnu[IF] = PyArray_DIM(preal, 1);

    if (vis.nnu[IF] > maxnchan) {
        maxnchan = vis.nnu[IF];
    }
    printf("inside uvmod.setData: vis.nnu[%d]: %i (maxchan = %i)\n", IF, vis.nnu[IF], maxnchan);

    Py_RETURN_TRUE;
}

// Fill-in the MODEL arrays
// USAGE FROM PYTHON: setModel(arrlist) where arrlist is list of data arrays.
PyObject *setModel(PyObject *self, PyObject *args)
{
    PyObject *HessArr, *GradArr, *modArr, *VarArr, *FixArr, *tempArr, *dparArr;
    PyObject *propRA, *propDec, *refpos, *parDep, *aG;

    int i, j, IF,isFixed, isMixed;

    if (!PyArg_ParseTuple(args, "OOOOOOOOOOOii", &modArr, &HessArr, &GradArr,
                          &VarArr, &FixArr, &dparArr,
                          &propRA, &propDec, &refpos,
                          &parDep, &aG,&isFixed, &isMixed)) {
        printf("FAILED setModel!\n");
        fflush(stdout);
        return NULL;
    }
    printf("\ninside uvmod.setModel: isFixed: %i, isMixed: %i\n", isFixed, isMixed);

    isModel = true;

    compFixed = isFixed == 1;
    MixedG = isMixed == 1;

    show_info("modArr", modArr);
    show_info("HessArr", HessArr);
    show_info("GradArr", GradArr);
    show_info("dparArr", dparArr);
    show_info("propRA", propRA);
    show_info("propDec", propDec);
    show_info("refpos", refpos);

    mod.models = (int *)PyArray_DATA(PyArray_FROM_OTF(modArr, NPY_INT32, NPY_IN_ARRAY));
    mod.Hessian = (double *)PyArray_DATA(PyArray_FROM_OTF(HessArr, NPY_DOUBLE, NPY_IN_ARRAY));
    mod.Gradient = (double *)PyArray_DATA(PyArray_FROM_OTF(GradArr, NPY_DOUBLE, NPY_IN_ARRAY));
    mod.dpar = (double *)PyArray_DATA(PyArray_FROM_OTF(dparArr, NPY_DOUBLE, NPY_IN_ARRAY));
    mod.muRA = (double *)PyArray_DATA(PyArray_FROM_OTF(propRA, NPY_DOUBLE, NPY_IN_ARRAY));
    mod.muDec = (double *)PyArray_DATA(PyArray_FROM_OTF(propDec, NPY_DOUBLE, NPY_IN_ARRAY));

    vis.phaseCenter = (double *)PyArray_DATA(PyArray_FROM_OTF(refpos, NPY_DOUBLE, NPY_IN_ARRAY));
    cosDecRef = cos(vis.phaseCenter[1]);
    sinDecRef = sin(vis.phaseCenter[1]);

    ncomp = PyArray_DIM(modArr, 0);
    npar = PyArray_DIM(GradArr, 0);
    Nants = (int)PyList_Size(parDep);
    printf("Nants = %d\n", Nants);

    mod.nparAnt.resize(Nants);
    mod.parAnt.resize(Nants);

    for (i = 0; i < Nants; i++) {
        mod.nparAnt[i] = PyArray_DIM(PyList_GetItem(parDep, i), 0);
        mod.parAnt[i] = (int *)PyArray_DATA(PyArray_FROM_OTF(PyList_GetItem(parDep, i), NPY_INT32, NPY_IN_ARRAY));
    }

    for (IF = 0; IF < Nspw; IF++) {
        mod.Gain[IF].resize(Nants);
        for (j = 0; j < Nants; j++) {
            mod.Gain[IF][j].resize(mod.nparAnt[j]);
            for (i = 0; i < mod.nparAnt[j]; i++) {
                tempArr = PyList_GetItem(PyList_GetItem(PyList_GetItem(aG, IF), j), i);
                // vsize = PyArray_SIZE(tempArr);
                mod.Gain[IF][j][i] = (cplx64 *)PyArray_DATA(PyArray_FROM_OTF(tempArr, NPY_COMPLEX128, NPY_IN_ARRAY));
            }
        }
    }

    HankelOrder = PyArray_DIM(PyList_GetItem(VarArr, 0), 1) - NparMax;
    printf("parameters: %d %d %d\n", HankelOrder, NparMax, npar);

    // delete[] mod.vars;
    // delete[] mod.fixp;
    mod.vars.resize(npar+1);
    mod.fixp.resize(npar+1);

    for (i = 0; i < (npar+1); i++) {
        tempArr = PyList_GetItem(VarArr, i);
        mod.vars[i] = (double *)PyArray_DATA(PyArray_FROM_OTF(tempArr, NPY_DOUBLE, NPY_IN_ARRAY));
    }

    for (i = 0; i < (npar+1); i++) {
        tempArr = PyList_GetItem(FixArr, i);
        mod.fixp[i] = (double *)PyArray_DATA(PyArray_FROM_OTF(tempArr, NPY_DOUBLE, NPY_IN_ARRAY));
    }

    Py_RETURN_TRUE;
}

// Allocate memory for the workers.
// USAGE FROM PYTHON: setWork() with no arguments.
//                    (setNCPU and setModel must be run first!)
static PyObject *setWork(PyObject *self, PyObject *args)
{
    int i;
    for (i = 0; i < NCPU; i++) {
        mod.WorkHess[i].resize(npar*npar);
        mod.WorkGrad[i].resize(npar);
    }
    //   printf("\n     setWork %i\n\n",NCPU);

    PyObject *ret = Py_BuildValue("i", 10);
    return ret;
}

// Deallocate the memory allocated with setWork.
// USAGE FROM PYTHON: unsetWork() with no arguments.
//                    (obviously, setWork must be run first!)
// static PyObject *unsetWork(PyObject *self, PyObject *args)
// {
//     // int i;
//     printf("\n UNSET WORK! %i\n", NCPU);
//     // for (i = 0; i < NCPU; i++) {
//     //     if (mod.WorkHess[i]) delete mod.WorkHess[i];
//     //     if (mod.WorkGrad[i]) delete mod.WorkGrad[i];
//     // }

//     // delete WorkHess;
//     // delete WorkGrad;
//     PyObject *ret = Py_BuildValue("i", 10);
//     return ret;
// }

/* Main Python function. It spreads the work through the workers */
// USAGE FROM PYTHON: modelcomp(IF, nui, opts) (see Python code for info).
static PyObject *modelcomp(PyObject *self, PyObject *args)
{
    void *status;
    double totChi = 0.0;
    int i,j;
    int nparsq = npar*npar;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "iii", &cIF, &nui, &mode)) {
        printf("FAILED modelcomp!\n");
        fflush(stdout);
        return NULL;
    }
    printf("\ninside uvmod.modelcomp: IF: %i nui: %i mode: %i\n", cIF, nui, mode);

    // Zero the workers memory:
    for (i = 0; i < NCPU; i++) {
        for (j = 0; j < nparsq; j++) {
            mod.WorkHess[i][j] = 0.0;
        }
        for (j = 0; j < npar; j++) {
            mod.WorkGrad[i][j] = 0.0;
        }
    }

    if (NCPU > 1) {
        /* Code for the case NCPU>1.
           Define the workers and perform the parallel task. */
        pthread_t MyThreads[NCPU];
        pthread_attr_t attr;

        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

        for (i = 0; i < NCPU; i++) {
            pthread_create(&MyThreads[i], &attr, writemod, (void *)&worker[i]);
        }
        pthread_attr_destroy(&attr);

        for (i = 0; i < NCPU; i++) {
            pthread_join(MyThreads[i], &status);
        }
    } else {
        /* Case of one single process (NCPU=1).
           Now, the master will do all the work. */
        master.t0[cIF] = 0;
        master.t1[cIF] = vis.nt[cIF];
        master.Iam = -1;

        writemod((void *)&master);
    }

    /* Add-up the Chi square, the error vector, and the  Hessian */
    for (i = 0; i < NCPU; i++) {
        totChi += mod.Chi2[i];
        for (j = 0; j < npar; j++) {
            mod.Gradient[j] += mod.WorkGrad[i][j];
        }
        for (j = 0; j < nparsq; j++) {
            mod.Hessian[j] += mod.WorkHess[i][j];
        }
    }

    /* Update references and set the return value */
    PyObject *ret = Py_BuildValue("d", totChi);

    return ret;
}

static PyObject *QuinnFF(PyObject *self, PyObject *args)
{
    int IFFit, refant, doModel, doGlobal;
    int ErrStat = 0; // To track errors in GFF (not implemented)

    if (!PyArg_ParseTuple(args, "iiii", &IFFit, &refant, &doModel, &doGlobal)) {
        printf("FAILED QuinnFringe!\n");
        fflush(stdout);
        return NULL;
    }
    PyObject *ret;

#if QUINN_FITTER == 0
    if (IFFit >= Nspw) {
        ret = Py_BuildValue("i", -1);
        printf("\n spw is too high!");
        return ret;
    }

    QuinnFringe *FringeFit = new QuinnFringe(Nants, vis.nt[IFFit], vis.nnu[IFFit],
                                             vis.ObsVis[IFFit], mod.ModVis[IFFit],
                                             vis.ants[0][IFFit], vis.ants[1][IFFit],
                                             vis.dt[IFFit],vis.fittable[IFFit],
                                             vis.freqs[IFFit],vis.wgt[0][IFFit]);
    int result = FringeFit->GFF(refant, doGlobal, doModel);
    if (result != 0) {
        ret = Py_BuildValue("i", result);
        return ret;
    }

    // Return the gains:
    double *Rates = new double[Nants];
    double *Delays = new double[Nants];
    double *Phases = new double[Nants];
    double *Bins = new double[2];

    ErrStat = FringeFit->getRates(Rates);
    if (ErrStat != 0) return NULL;
    ErrStat = FringeFit->getDelays(Delays);
    if (ErrStat != 0) return NULL;
    ErrStat = FringeFit->getPhases(Phases);
    if (ErrStat != 0) return NULL;
    ErrStat = FringeFit->getBins(Bins);
    if (ErrStat != 0) return NULL;

    PyObject *PyRate, *PyDelay, *PyPhase;

    npy_intp Dims[1];
    Dims[0] = Nants;
    printf("\nNants: %i\n",Nants);

    PyRate = PyArray_SimpleNewFromData(1, Dims, NPY_FLOAT64, (void *)Rates);
    PyDelay = PyArray_SimpleNewFromData(1, Dims, NPY_FLOAT64, (void *)Delays);
    PyPhase = PyArray_SimpleNewFromData(1, Dims, NPY_FLOAT64, (void *)Phases);

    //printf("\n Ant 1: %.4e %.4e %.4e\n", Rates[0], Delays[0], Phases[0]);
    ret = Py_BuildValue("[O,O,O,d,d]", PyDelay, PyRate, PyPhase, Bins[0], Bins[1]);
    if (FringeFit) delete FringeFit;
#else
    printf("\n QUINN FITTER NOT INSTALLED!");
    ret = Py_BuildValue("i", -1);
#endif
    return ret;
}
