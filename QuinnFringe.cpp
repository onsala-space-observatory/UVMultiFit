/*
   QUINNFRINGE - A GLOBAL FRINGE FITTER BASED ON DELAY/RATE MATRICES.
   Copyright (C) 2018  Ivan Marti-Vidal
   Nordic Node of EU ALMA Regional Center (Onsala, Sweden)
   Max-Planck-Institut fuer Radioastronomie (Bonn, Germany)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <fftw3.h>
#include <gsl/gsl_linalg.h>
#include "QuinnFringe.h"

typedef std::complex<double> cplx64;

QuinnFringe::QuinnFringe(int Na, int Nt, int Nc,
                         cplx64 *ObsV,
                         cplx64 *ModV,
                         int *A1,
                         int *A2,
                         double *tA,
                         int8_t *dofit,
                         double *freqs,
                         double *wgts) :
    Nant(Na), Ntime(Nt), Nchan(Nc), ObsVis(ObsV), ModVis(ModV), Ant1(A1), Ant2(A2),
    Freqs(freqs), fittable(dofit), Times(tA), DataWeights(wgts)

{
    NBas = Nant*Nant;

    // Allocate memory for antenna gains:
    // Memory is not freed (i.e., the arrays from previous runs already belong to Python):
    Rates = new double[Nant];
    Delays = new double[Nant];
    Phases = new double[Nant];

    Dnu = Freqs[Nchan-1]-Freqs[0];

    NData = new int[NBas];

    int i, j, k;
    k = 0;

    BasNum = new int*[Nant];
    for (i = 0; i < Nant; i++) {
        BasNum[i] = new int[Nant];
        for (j = 0; j < Nant; j++) {
            BasNum[i][j] = k;
            k += 1;
        }
    }

}

QuinnFringe::~QuinnFringe() {
    int i;
    delete[] Rates;
    delete[] Delays;
    delete[] Phases;
    delete[] NData;
    for (i = 0; i < Nant; i++) {
        delete[] BasNum[i];
    }
    delete[] BasNum;
}

int QuinnFringe::getRates(double *OutRat) {
    int i;
    for (i = 0; i < Nant; i++) {
        OutRat[i] = Rates[i];
    }
    return 0;
}

int QuinnFringe::getDelays(double *OutDel) {
    int i;
    for (i = 0; i < Nant; i++) {
        OutDel[i] = Delays[i];
    }
    return 0;
}

int QuinnFringe::getPhases(double *OutPhs) {
    int i;
    for (i = 0; i < Nant; i++) {
        OutPhs[i] = Phases[i];
    }
    return 0;
}

int QuinnFringe::getBins(double *OutBins) {
    OutBins[0] = 1./Dnu;
    OutBins[1] = 1./DtMin;
    return 0;
}

// Quinn Estimator of the FFT peak with sub-bin precision:
double QuinnFringe::Tau(double x) {
    return 0.25*log(3.*x*x + 6.*x + 1.) - sqrt(6.)/24.*log((x+1.-sqrt(2./3.))/(x+1.+sqrt(2./3.)));
}

double QuinnFringe::Estimate(cplx64 *FFTVec) {
    double Denom = FFTVec[1].real()*FFTVec[1].real() + FFTVec[1].imag()*FFTVec[1].imag();
    double AP = (FFTVec[2].real()*FFTVec[1].real() + FFTVec[2].imag()*FFTVec[1].imag())/Denom;
    double AM = (FFTVec[0].real()*FFTVec[1].real() + FFTVec[0].imag()*FFTVec[1].imag())/Denom;

    double DP = -AP/(1.-AP);
    double DM =  AM/(1.-AM);

    return (DP + DM)/2. + Tau(DP*DP) - Tau(DM*DM);
}

int QuinnFringe::GFF(int REFANT, int DOGLOBAL, int DOMODEL) {

    int i, j, k, l, t, inu, knu;

    // Will implement window search in the future.
    // For now, no effective window:
    int npix = 10000;

    // Will implement SNR cutoff in the future.
    // For now, no SNR cutoff:
    double minSNR = 0.0;

    // Get The dimensions for the matrices:
    for (i = 0; i < NBas; i++) {
        NData[i] = 0;
    }

    int i0, i1;
    switch (DOGLOBAL) {
      case 0:
        i0 = REFANT;
        i1 = REFANT+1;
        break;
      case 1:
        i0 = 0;
        i1 = Nant;
        break;
      default:
        printf("\n UNRECOGNIZED GLOBAL OPTION. WILL NOT GLOBALIZE!");
        i0 = REFANT;
        i1 = REFANT+1;
    }

    for (i = i0; i < i1; i++) {
        for (j = 0; j < Nant; j++) {
            for (t = 0; t < Ntime; t++) {
                if (Ant1[t] == i && Ant2[t] == j && i != j && fittable[t]) {
                    NData[BasNum[i][j]] += 1;
                }
            }
        }
    }

    cplx64 *aroundPeak = new cplx64[3];

    double* BLRates = new double[NBas];
    double* BLDelays = new double[NBas];
    double* BLPhases = new double[NBas];
    double* Weights = new double[NBas];

    // FFT FOR EACH BASELINE:
    int MaxDim = Nchan*NData[0];
    double BasWgt;
    int prevNvis = NData[0];
    double tmin,tmax;

    for (j = 1; j < NBas; j++) {
        if (NData[j]*Nchan>MaxDim) {
            MaxDim=NData[j]*Nchan;
        }
    }

    int *Tindex[NBas];
    for (j = 1; j < NBas; j++) {
        i = NData[j];
        i += i==0?1:0;
        Tindex[j] = new int[i];
        for (k = 0; k < i; k++) {
            Tindex[j][k] = -1;
        }
    }

    fftw_complex *BufferVis = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * MaxDim);
    fftw_complex *AUX = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * MaxDim);
    fftw_complex *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * MaxDim);
    fftw_plan p = fftw_plan_dft_2d(NData[0], Nchan, BufferVis, out, FFTW_FORWARD, FFTW_MEASURE);

    cplx64 *Temp;
    cplx64 *BufferC;

    Temp = reinterpret_cast<cplx64*>(out);
    BufferC = reinterpret_cast<cplx64*>(BufferVis);

    int a1, a2, a1sel, a2sel, BNum;
    a1= -1;
    a2= -1;
    DtMin = 10000.0;

    for (j = 0; j < NBas; j++) {

        BLRates[j] = 0.0;
        BLDelays[j] = 0.0;
        BLPhases[j] = 0.0;

        int NcurrVis = 0;

        BasWgt = 0.0;
        tmin = 1.e9;
        tmax = 0.;

        if (NData[j] != 0) {
            // Arrange data for this baseline:
            for (k = 0; k < Ntime; k++) {
                a1sel = Ant1[k];
                a2sel = Ant2[k];

                BNum = BasNum[a1sel][a2sel];

                if (BNum==j && fittable[k]) {
                    a1 = a1sel;
                    a2 = a2sel;
                    inu = NcurrVis*Nchan;
                    knu = k*Nchan;
                    if (DOMODEL == 0) {
                        memcpy(&BufferC[inu],&ObsVis[knu],Nchan*sizeof(cplx64));
                        for (i = 0; i < Nchan; i++) {
                            BasWgt += DataWeights[knu + i];
                            if (DataWeights[knu+i]<=0.0) {
                                BufferC[inu + i] = 0.0;
                            }
                        }
                    } else {
                        for (i = 0; i < Nchan; i++) {
                            BufferC[inu + i] = ObsVis[knu + i]/ModVis[knu + i];
                            BasWgt += std::abs(ModVis[knu + i])*DataWeights[knu + i];
                            if (DataWeights[knu+i]<=0.0) {
                                BufferC[inu + i] = 0.0;
                            }
                        }
                    }
                    Tindex[j][NcurrVis] = k;
                    NcurrVis += 1;
                    if (Times[k]>tmax) {
                        tmax=Times[k];
                    }
                    if (Times[k]<tmin) {
                        tmin=Times[k];
                    }
                }
            }
            BasWgt /= ((double) Nchan);
            Dtime = tmax-tmin;
            if (Dtime>0.0 && DtMin>Dtime) {
                DtMin = Dtime;
            }

            // FFT the fringe and find the peak:
            if (NcurrVis > 1) {

                // Re-define the FFTW plan if dimensions changed:
                if (NcurrVis != prevNvis) {
                    prevNvis = NcurrVis;
                    memcpy(&AUX[0],&BufferVis[0],NcurrVis*Nchan*sizeof(fftw_complex));
                    fftw_destroy_plan(p);
                    p = fftw_plan_dft_2d(NcurrVis, Nchan, BufferVis, out, FFTW_FORWARD, FFTW_MEASURE);
                    memcpy(&BufferVis[0],&AUX[0],NcurrVis*Nchan*sizeof(fftw_complex));
                }
                fftw_execute(p);

                double Peak = 0.0;
                double AbsP;
                double rmsFringe = 0.0;
                double AvgFringe = 0.0;
                // double LastPeak = 0.0;
                double Dnpix = (double) Nchan*NcurrVis;
                double FringeSNR;

                int nu[3],ti[3],row;
                int Chi, Chf, t0, tf;

                if (Nchan>npix) {
                    Chi = npix/2;
                    Chf = Nchan-npix/2;
                }
                else  {
                    Chi = Nchan/2;
                    Chf = Nchan/2;
                }

                if (NcurrVis>npix) {
                    t0 = npix/2;
                    tf = NcurrVis-npix/2;
                }
                else  {
                    t0 = NcurrVis/2;
                    tf = NcurrVis/2;
                }

                nu[1] = 0;
                ti[1] = 0;

                // First Quadrant:
                for (l = 0; l < t0; l++) {
                    row = l*Nchan;
                    for (k = 0; k < Chi; k++) {
                        AbsP = std::abs(Temp[row + k]);
                        if (AbsP>Peak) {
                            Peak = AbsP;
                            nu[1] = k;
                            ti[1] = l;
                        }
                        rmsFringe += Peak*Peak;
                        AvgFringe += Peak;
                    }
                }
                // Second Quadrant:
                for (l = tf; l < NcurrVis; l++) {
                    row = l*Nchan;
                    for (k = 0; k < Chi; k++) {
                        AbsP = std::abs(Temp[row + k]);
                        if (AbsP>Peak) {
                            Peak = AbsP;
                            nu[1] = k;
                            ti[1] = l;
                        }
                    }
                    rmsFringe += Peak*Peak;
                    AvgFringe += Peak;
                }
                // Third Quadrant:
                for (l = 0; l < t0; l++) {
                    row = l*Nchan;
                    for (k = Chf; k < Nchan; k++) {
                        AbsP = std::abs(Temp[row + k]);
                        if (AbsP>Peak) {
                            Peak = AbsP;
                            nu[1] = k;
                            ti[1] = l;
                        }
                        rmsFringe += Peak*Peak;
                        AvgFringe += Peak;
                    }
                }
                // Fourth Quadrant:
                for (l = tf; l < NcurrVis; l++) {
                    row = l*Nchan;
                    for (k = Chf; k < Nchan; k++) {
                        AbsP = std::abs(Temp[row + k]);
                        if (AbsP>Peak) {
                            Peak = AbsP;
                            nu[1] = k;
                            ti[1] = l;
                        }
                        rmsFringe += Peak*Peak;
                        AvgFringe += Peak;
                    }
                }

                // Unwrap:
                if (nu[1]==0) {
                    nu[0]=Nchan-1;
                }
                else {
                    nu[0]=nu[1]-1;
                }
                if (nu[1]==Nchan-1) {
                    nu[2]=0;
                }
                else {
                    nu[2]=nu[1]+1;
                }
                if (ti[1]==0) {
                    ti[0]=NcurrVis-1;
                }
                else {
                    ti[0]=ti[1]-1;
                }
                if (ti[1]==NcurrVis-1) {
                    ti[2]=0;
                }
                else {
                    ti[2]=ti[1]+1;
                }

                // Get the SNR of the fringe (i.e., peak over RMS, but without the peak):
                AbsP = std::abs(Temp[ti[1]*Nchan + nu[1]]);
                rmsFringe -= AbsP*AbsP;
                AvgFringe -= AbsP;

                AbsP = std::abs(Temp[ti[0]*Nchan + nu[1]]);
                rmsFringe -= AbsP*AbsP;
                AvgFringe -= AbsP;
                AbsP = std::abs(Temp[ti[2]*Nchan + nu[1]]);
                rmsFringe -= AbsP*AbsP;
                AvgFringe -= AbsP;

                AbsP = std::abs(Temp[ti[1]*Nchan + nu[0]]);
                rmsFringe -= AbsP*AbsP;
                AvgFringe -= AbsP;
                AbsP = std::abs(Temp[ti[1]*Nchan + nu[2]]);
                rmsFringe -= AbsP*AbsP;
                AvgFringe -= AbsP;

                FringeSNR = Peak/pow(rmsFringe/(Dnpix-5.) - pow(AvgFringe/(Dnpix-5.),2.),0.5);
                //  printf("BAS %i - %i. SNR: %.3e  %i\n",a1,a2,FringeSNR, NcurrVis);

                // SNR cutoff:
                if (FringeSNR<minSNR) {
                    FringeSNR=0.0;
                }

                // Estimate the rate with sub-bin precision:
                aroundPeak[0] = Temp[nu[1] + Nchan*(ti[0])];
                aroundPeak[1] = Temp[nu[1] + Nchan*(ti[1])];
                aroundPeak[2] = Temp[nu[1] + Nchan*(ti[2])];

                BLRates[BasNum[a1][a2]] = ((double) ti[1]);
                BLRates[BasNum[a1][a2]] += Estimate(aroundPeak);

                if (BLRates[BasNum[a1][a2]] > ((double) NcurrVis)/2.) {
                    BLRates[BasNum[a1][a2]] = BLRates[BasNum[a1][a2]] - (double) NcurrVis;
                }

                BLRates[BasNum[a1][a2]] /= Dtime;

                // Estimate the delay with sub-bin precision:
                aroundPeak[0] = Temp[nu[0] + Nchan*(ti[1])];
                aroundPeak[1] = Temp[nu[1] + Nchan*(ti[1])];
                aroundPeak[2] = Temp[nu[2] + Nchan*(ti[1])];

                BLDelays[BasNum[a1][a2]] = ((double) nu[1]);
                BLDelays[BasNum[a1][a2]] += Estimate(aroundPeak);

                if (BLDelays[BasNum[a1][a2]] > ((double) Nchan)/2.) {
                    BLDelays[BasNum[a1][a2]] = BLDelays[BasNum[a1][a2]] - (double) Nchan;
                }

                BLDelays[BasNum[a1][a2]] /= Dnu;

                Weights[BasNum[a1][a2]] = BasWgt*FringeSNR;
            } else {   // Comes from if (NcurrVis > 1)
                BLRates[BasNum[a1][a2]] = 0.0;
                BLRates[BasNum[a1][a2]] = 0.0;
                BLDelays[BasNum[a1][a2]] = 0.0;
                BLDelays[BasNum[a1][a2]] = 0.0;
                Weights[BasNum[a1][a2]] = 0.0;
                printf("WARNING! BASELINE %i-%i HAS NO DATA!\n",i,k);
            }

        }
    }

    // Tell UVMULTIFIT about the bad data (i.e., low-SNR fringes):
    for (i = 0; i < NBas; i++) {
        if (Weights[i] == 0.0) {
            for (j = 0; j < NData[i]; j++) {
                if (Tindex[i][j]>=0) {
                    fittable[Tindex[i][j]]=false;
                }
            }
        }
    }

    // Globalize the rate and delay solutions:
    NcalAnt = 0;
    calAnt = new int[Nant];
    bool GoodRef = false;

    for (i = 0; i < Nant; i++) {
        for (j = 0; j < Nant; j++) {
            if ((NData[BasNum[i][j]]>0 || NData[BasNum[j][i]]>0)
                && i!=j && (Weights[BasNum[i][j]]>0. ||  Weights[BasNum[j][i]]>0.)) {
                if (i==REFANT || j==REFANT) {
                    GoodRef=true;
                }
                if (i!=REFANT) {
                    calAnt[NcalAnt] = i;
                    NcalAnt ++;
                    break;
                }
            }
        }
    }

    if (!GoodRef) { // Refant is bad. Take the last in list as new REFANT:
        REFANT = calAnt[NcalAnt-1];
        NcalAnt -= 1;
    }
    printf("# of free antenna gains: %i\n",NcalAnt);

    int NBasFit = (Nant-1)*(Nant-1);

    double *Hessian = new double[NBasFit];
    double *RateResVec = new double[NcalAnt];
    double *DelResVec = new double[NcalAnt];
    double *CovMat = new double[NBasFit];

    for (i = 0; i < NBasFit; i++) {
        Hessian[i] = 0.0;
        CovMat[i] = 0.0;
    }

    for (i = 0; i < NcalAnt; i++) {
        RateResVec[i] = 0.0;
        DelResVec[i] = 0.0;
    }

    int ca1, ca2;
    ca1 = -1;
    ca2 = -1;

    if (Nant>2) {
        for (i = i0; i < i1; i++) {
            for (j = 0; j < Nant; j++) {
                if (NData[BasNum[i][j]]>0 && i != j) {
                    a1=i;
                    a2=j;
                    BNum = BasNum[i][j];
                    if (a1 != REFANT) {
                        for (ca1 = 0; ca1 < NcalAnt; ca1++) {
                            if (calAnt[ca1]==a1) {
                                break;
                            }
                        }
                        RateResVec[ca1] += Weights[BNum]*BLRates[BNum];
                        DelResVec[ca1] += Weights[BNum]*BLDelays[BNum];
                        Hessian[ca1*NcalAnt + ca1] += Weights[BNum];
                    }
                    if (a2 != REFANT) {
                        for (ca2 = 0; ca2 < NcalAnt; ca2++) {
                            if (calAnt[ca2]==a2) {
                                break;
                            }
                        }
                        RateResVec[ca2] -= Weights[BNum]*BLRates[BNum];
                        DelResVec[ca2] -= Weights[BNum]*BLDelays[BNum];
                        Hessian[ca2*NcalAnt + ca2] += Weights[BNum];
                    }
                    if (a1 != REFANT && a2 != REFANT) {
                        Hessian[ca1*NcalAnt + ca2] += -Weights[BNum];
                        Hessian[ca2*NcalAnt + ca1] += -Weights[BNum];
                    }
                }
            }
        }

    } else {
        Hessian[0] += Weights[0];
        RateResVec[0] += Weights[0]*BLRates[0];
        DelResVec[0] += Weights[0]*BLDelays[0];
    }

    /*
    printf("\n\n GLOBALIZATION HESSIAN:\n");
    for (i = 0;i < NcalAnt;i++){
      printf("\n");
      for (j = 0;j < NcalAnt;j++){
        printf("%.2e ",Hessian[i*NcalAnt + j]);
      }
    }

    printf("\n DONE \n");
    */


    // The Hessian's inverse can be reused for rates and delays!
    gsl_matrix_view m = gsl_matrix_view_array (Hessian, NcalAnt, NcalAnt);
    gsl_matrix_view inv = gsl_matrix_view_array(CovMat,NcalAnt,NcalAnt);

    int s;

    gsl_permutation *perm = gsl_permutation_alloc (NcalAnt);

    gsl_linalg_LU_decomp (&m.matrix, perm, &s);
    gsl_linalg_LU_invert (&m.matrix, perm, &inv.matrix);


    // Derive the rates as CovMat*RateVec and delays as CovMat*DelVec:
    for (i = 0; i < Nant; i++) {
        Phases[i] = 0.0;
        Rates[i] = 0.0;
        Delays[i] = 0.0;
    }

    if (Nant>2) {
        for (i = 0; i < NcalAnt; i++) {
            for (j = 0; j < NcalAnt; j++) {
                Rates[calAnt[i]] += RateResVec[j]*gsl_matrix_get(&inv.matrix,i,j);
                Delays[calAnt[i]] += DelResVec[j]*gsl_matrix_get(&inv.matrix,i,j);
            }
        }
    }

    for (i = 0; i < Nant; i++) {
        printf("Antenna %i -> Rate: %+.3e Hz; Delay: %+.3e s.\n",i,Rates[i],Delays[i]);
    }

    printf("Estimating phases (not globalized)...\n");

    double Phase, DTau, DRate;
    cplx64 CalPhase;

    DTau = Dnu/((double) Nchan);
    DRate = DtMin/((double) (MaxDim/Nchan));

    for (i = 0; i < NcalAnt; i++) {
        CalPhase = cplx64(0.0);

        j = BasNum[REFANT][calAnt[i]];
        if (NData[j] > 0) {
            for (t = 0; t < NData[j]; t++) {
                for (k = 0; k < Nchan; k++) {
                    Phase = + Rates[calAnt[i]]*((double) t)*DRate + Delays[calAnt[i]]*((double) k)*DTau;
                    CalPhase += ObsVis[Tindex[j][t]*Nchan + k] * std::polar(1., 6.283185*Phase);
                }
            }
            Phases[calAnt[i]] = -std::arg(CalPhase);
        } else {
            j = BasNum[calAnt[i]][REFANT];
            for (t = 0; t < NData[j]; t++) {
                for (k = 0; k < Nchan; k++) {
                    Phase = - Rates[calAnt[i]]*((double) t)*DRate - Delays[calAnt[i]]*((double) k)*DTau;
                    CalPhase += ObsVis[Tindex[j][t]*Nchan + k] * std::polar(1.,6.283185*Phase);
                }
            }
            Phases[calAnt[i]] = std::arg(CalPhase);
        }

    }

    delete[] aroundPeak;
    delete[] BLRates;
    delete[] BLDelays;
    delete[] BLPhases;
    delete[] Weights;

    for (j = 1; j < NBas; j++) {
        delete[] Tindex[j];
    }

    delete[] calAnt;

    delete[] Hessian;
    delete[] RateResVec;
    delete[] DelResVec;
    delete[] CovMat;

    fftw_free(BufferVis);
    fftw_free(out);

    return 0;
}
