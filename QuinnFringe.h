#include <sys/types.h>
#include <stdlib.h>
#include <complex>

typedef std::complex<double> cplx64;

class QuinnFringe {

 public:
    QuinnFringe(int Nant, int Ntime, int nnu,
                cplx64 *ObsVis,
                cplx64 *ModVis,
                int *Ant1,
                int *Ant2,
                double *t,
                int8_t *dofit,
                double *freqs,
                double *wgts);
    ~QuinnFringe();

    int GFF(int refant, int doGlobal, int doModel);

    int getRates(double *Rates);
    int getDelays(double *Delays);
    int getPhases(double *Phases);
    int getBins(double *Bins);

    static double Tau(double x);
    static double Estimate(cplx64 *FFTVec);

  private:
    int Nant, Ntime, Nchan;
    cplx64 *ObsVis;
    cplx64 *ModVis;
    int *Ant1;
    int *Ant2;
    double *Freqs;
    int8_t *fittable;
    double *Times;
    double *DataWeights;

    int *NData;
    int **BasNum;
    double *Phases, *Rates, *Delays;
    double Dnu, Dtime,DtMin;
    int NcalAnt, NBas;
    int *calAnt;
};
