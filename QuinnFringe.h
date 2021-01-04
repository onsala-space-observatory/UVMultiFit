#include <sys/types.h>
#include <stdlib.h>
#include <complex>

typedef std::complex<double> cplx64f;

class QuinnFringe {

 public:
    QuinnFringe(int Nant, int Ntime, int nnu, cplx64f *ObsVis, cplx64f *ModVis,
                int *Ant1, int *Ant2, double *t, char *dofit, double *freqs, double *wgts);
    ~QuinnFringe();

    int GFF(int refant, int doGlobal, int doModel);

    int getRates(double *Rates);
    int getDelays(double *Delays);
    int getPhases(double *Phases);
    int getBins(double *Bins);

  private:
    cplx64f *ObsVis, *ModVis;
    int *Ant1, *Ant2, *calAnt;

    int Nant, NcalAnt, Nchan, Ntime, NBas;

    int *NData;
    int **BasNum;
    char *fittable;
    double *Phases, *Rates, *Delays, *Freqs, *Times, *DataWeights;
    double Dnu, Dtime,DtMin;
};
