import sys
import os
import numpy as np
import logging

from casatools import image

from NordicARC import measurementset, modeler

import matplotlib.pyplot as plt

def save_results(outfile, results, mdl, ms):
    Nvis = results["Degrees of Freedom"]
    if isinstance(Nvis, list):
        DOG = int(np.average(Nvis))
    else:
        DOG = int(Nvis)
    ChiSq = results["Reduced Chi squared"]
    fitparams = results["Parameters"]
    fiterrors = results["Uncertainties"]
    if not mdl.OneFitPerChannel:
        prtpars = []
        for pp, ppar in enumerate(fitparams):
            prtpars.append(ppar)
            prtpars.append(fiterrors[pp])
    else:
        nspwtot = ms.spwlist[-1][3] + len(ms.spwlist[-1][2])
        prtpars = [[[] for spi in range(len(fitparams[sp]))] for sp in range(nspwtot)]
        for sp in range(nspwtot):
            for pp, ppar in enumerate(fitparams[sp]):
                for ppi, ppari in enumerate(ppar):
                    prtpars[sp][pp].append(ppari)
                    prtpars[sp][pp].append(fiterrors[sp][pp][ppi])

    with open(outfile, 'w') as outf:
        outf.write("# MODELFITTING RESULTS FOR MS: " + ','.join(ms.vis) + "\n")
        outf.write("# NOTICE THAT THE REDUCED CHI SQUARE SHOWN HERE IS THE VALUE\n")
        outf.write("# *BEFORE* THE RE-SCALING OF THE VISIBILITY WEIGHTS.\n")
        outf.write("# HOWEVER, THE PARAMETER UNCERTAINTIES ARE *ALREADY* GIVEN\n")
        outf.write("# FOR A REDUCED CHI SQUARE OF 1.\n")
        outf.write(f"# AVG. NUMBER OF DEGREES OF FREEDOM: {DOG}\n")
        if ms.pbeam:
            if isinstance(ms.dish_diameter, float):
                outf.write("# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. "
                           "USING A DISH DIAMETER OF: {ms.dish_diameter:.3f} METERS\n")
            else:
                outf.write("# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. USING: \n")
                for i, name in enumerate(ms.antnames):
                    outf.write(f"#    FOR ANTENNA {name} A DIAMETER OF {ms.userDiameters[i]:.2f}m \n")

        else:
            outf.write("# PRIMARY-BEAM CORRECTION HAS NOT BEEN APPLIED.\n")

        outf.write("###########################################\n")
        outf.write("#\n# MODEL CONSISTS OF:\n")
        for m, mod in enumerate(mdl.model):
            outf.write("# '" + mod + "' with variables: " + mdl.var[m] + "\n")
        if len(mdl.fixed) > 0:
            outf.write("#\n#\n# FIXED MODEL CONSISTS OF:\n")
            if mdl.takeModel:
                outf.write("# THE MODEL COLUMN FOUND IN THE DATA\n")
            else:
                for m, mod in enumerate(mdl.fixed):
                    var2print = list(map(float, mdl.fixedvar[m].split(',')))
                    outf.write(
                        "# '" + mod + "' with variables: " + ' '.join(['%.3e'] * len(var2print)) % tuple(var2print))
                outf.write(f"#\n#  - AND SCALING FACTOR: {mdl.scalefix}\n")

        outf.write("#\n# INITIAL PARAMETER VALUES:\n")
        for p0i, p0 in enumerate(mdl.p_ini):
            if mdl.bounds is not None:
                outf.write(f"#  p[{p0i}] = {p0:.5e} with bounds: {str(mdl.bounds[p0i])}\n")
            else:
                outf.write(f"#  p[{p0i}] = {p0:.5e} with no bounds\n")
        outf.write("#\n##########################################\n")

        N = len(mdl.p_ini)
        parshead = [n for n in range(N) for i in range(2)]

        headstr = ("# Frequency (Hz)     " + "  p[%i]       err[%i]   " * N + "   red. ChiSq\n") \
            % tuple(parshead)
        outf.write(headstr)

        formatting = "%.12e   " + "%.4e " * (2 * N) + "   %.4e \n"
        if not mdl.OneFitPerChannel:
            toprint = tuple([np.average(ms.averfreqs)] + prtpars + [ChiSq])
            outf.write(formatting % toprint)

        else:
            for spwvi in ms.spwlist:
                for r, rr in enumerate(spwvi[2]):
                    k = spwvi[3] + r
                    # rang = rr[1]
                    for nu, freq in enumerate(ms.averfreqs[k]):
                        toprint = tuple([freq] + prtpars[k][nu] + [ChiSq[k][nu]])
                        outf.write(formatting % toprint)

def uvmultifit(vis='', spw='0', column='data', field=0, scans=[],
               uniform=False, chanwidth=1, timewidth=1, stokes='I',
               write='', MJDrange=[-1.0, -1.0], model=['delta'],
               var=['p[0],p[1],p[2]'], p_ini=[0.0, 0.0, 1.0], phase_center='',
               fixed=[], fixedvar=[], scalefix='1.0', outfile='modelfit.dat',
               NCPU=4, pbeam=False, ldfac=1.22, dish_diameter=0.0,
               OneFitPerChannel=True, bounds=None, cov_return=False,
               finetune=False, uvtaper=0.0, method='levenberg', wgt_power=1.0,
               LMtune=[1.0e-3, 10.0, 1.0e-5, 200, 1.0e-3], SMPtune=[1.0e-4, 1.0e-1, 200],
               only_flux=False, proper_motion=0.0, HankelOrder=80,
               amp_gains={}, phase_gains={}):
    logger = logging.getLogger("uvmultifit")
    logging.basicConfig(level=logging.INFO,
                        format='%(name)s - %(levelname)s - %(message)s')
    logger.info("started")
    mdl = modeler.Modeler(model, var, p_ini, bounds, OneFitPerChannel, fixed, fixedvar,
                            scalefix, phase_gains, amp_gains, method, LMtune, SMPtune, proper_motion)
    if not isinstance(only_flux, bool):
        logger.error("'only_flux' must be boolean")
        mdl.fluxOnly = False
    else:
        mdl.fluxOnly = only_flux

    corrected = column == "corrected"
    ms = measurementset.MeasurementSet(vis, spw, field, scans, corrected, uniform,
                                       uvtaper, chanwidth, timewidth, stokes,
                                       MJDrange, phase_center, pbeam, wgt_power, dish_diameter)
    if not ms.check_measurementset():
        logger.error("failed to check measurement set")
        sys.exit(1)

    if not ms.read_data():
        logger.error("failed to read data")
        sys.exit(1)

    nspw = ms.Nspw
    if not mdl.init_data(nspw, ms):
        logger.error("failed to init data")
        sys.exit(1)

    if not mdl.init_model(nspw, only_flux, ms.refpos):
        logger.error("failed to init model")
        sys.exit(1)

    write_model = {'': 0, 'model': 1, 'residuals': 2, 'calibrated': 3}
    write_model_index = write_model.get(write, -1)
    if write_model_index == -1:
        logger.error("keyword 'write' should be set to either '', 'model', 'residuals' or 'calibrated'")

    results = mdl.fit(ms, write_model, cov_return)
    if not results:
        logger.error("failed to fit model")
        sys.exit(1)

    if not finetune:
        #       if not self.fit():
        #             logger.error("failed fit!")
        #             return False
        if write_model_index in [1, 2, 3]:
            if timewidth == 1 and stokes not in ['Q', 'U', 'V']:
                logger.info("writing into measurement set(s)")
                # writeModel()
            else:
                msg = "cannot write into measurement set!\n"
                msg += "If you want to fill-in the model (or corrected) column:\n"
                msg += "    1.- 'timewidth' and 'chanwidth' should be set to 1\n"
                msg += "    2.- 'stokes' should be set to either 'I', 'PI', or a corr. product."
                logger.error(msg)

        logger.info("fit done!")

    column = "MODEL_DATA" if write_model_index == 1 else "CORRECTED_DATA"
    ms.write_model(column, mdl)

    # print(mdl)
    # mdl.dump()
    #
    # print(ms)
    # ms.dump()
    save_results(outfile, results, mdl, ms)

    return results

vis = 'Disc/Disc.alma.out10.noisy.ms'
model = ['disc']
Nu = '50.0GHz'
NuF = float(Nu.split('G')[0])*1.e9
modvars = "0,0, p[0]*(nu/%.4e)**p[1], p[2], p[3], p[4]" % NuF
modbounds = [[0.0, None], [-2.0, 2.0], [0.0, None], [0.1, 0.9], [0.0, 180.0]]

si = ['disc', 1.0, '0.2arcsec', '0.1arcsec', '60.0deg',
      Nu, 1.0, "J2000 10h00m00.0s -30d00m00.0s"]
size = float(si[2].split('a')[0])
minor = float(si[3].split('a')[0])
initial = [0.8, 0.0, size*1.2, minor/size*0.8, 45.0]

if False:
    from casatasks import clearcal
    from casatasks import tclean

    Cfile = 'TEST1.CLEAN'
    Rfile = 'TEST1.RESIDUALS'
    clearcal(vis)

    os.system('rm -rf %s.*' % Cfile)
    os.system('rm -rf %s.*' % Rfile)

    Npix = 1000
    cell = '0.01arcsec'
    cellN = float(cell.replace("arcsec", ""))
    tclean(vis, imagename=Cfile, cell=cell, imsize=2*Npix, niter=0)

    ia = image()
    ia.open('%s.image' % Cfile)
    resdat = ia.getchunk()[:, :, 0, 0]
    ia.close()

    fig = plt.figure()
    sub = fig.add_subplot(111)
    IMsize = cellN*Npix/2.
    ims = sub.imshow(np.transpose(resdat), origin='lower',
                     extent=[IMsize, -IMsize, -IMsize, IMsize])
    sub.set_xlabel('RA offset (arcsec)')
    sub.set_ylabel('Dec offset (arcsec)')
    cb = plt.colorbar(ims)
    cb.set_label('Jy/beam')
    # plt.savefig('%s.png' % Cfile)
    plt.show()
    impeak = np.max(resdat)

pbeam = False
r = uvmultifit(vis=vis, spw='0',
               model=model, OneFitPerChannel=False,
               var=modvars, write='residuals',
               p_ini=initial, pbeam=pbeam,
               bounds=modbounds)

S = np.array([2.435e-01, 1.000e+00, 2.000e-01, 5.000e-01, 6.000e+01])

print(f"disc flux at 50GHz (Jy): {r['Parameters'][0]:7.3f} +/- {r['Uncertainties'][0]:.3f}, true: {S[0]:7.3f}")
print(f"disc spectral index:     {r['Parameters'][1]:7.3f} +/- {r['Uncertainties'][1]:.3f}, true: {S[1]:7.3f}")
print(f"disc size (as):          {r['Parameters'][2]:7.3f} +/- {r['Uncertainties'][2]:.3f}, true: {S[2]:7.3f}")
print(f"disc axis ratio:         {r['Parameters'][3]:7.3f} +/- {r['Uncertainties'][3]:.3f}, true: {S[3]:7.3f}")
print(f"disc pos.angle (deg):    {r['Parameters'][4]:7.3f} +/- {r['Uncertainties'][4]:.3f}, true: {S[4]:7.3f}")
