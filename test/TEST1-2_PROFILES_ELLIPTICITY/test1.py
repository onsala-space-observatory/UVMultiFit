import time

tic = time.time()
myfit = uvm.uvmultifit(vis=visname, spw='0',
                       model=modelshape, OneFitPerChannel=False,
                       var=modvars, write='residuals',
                       p_ini=pini,
                       bounds=parbound)

msg = "\nDisc Flux at 50GHz (Jy): %.4f +/- %.4f ; True: %.3f " % \
      (myfit.result['Parameters'][0], myfit.result['Uncertainties'][0], S[0])
msg += "\nDisc Spectral Index:  %.4f +/- %.4f ; True: %.3f " % \
       (myfit.result['Parameters'][1], myfit.result['Uncertainties'][1], S[1])
msg += "\nDisc Size (as):  %.4f +/- %.4f ; True: %.3f " % \
       (myfit.result['Parameters'][2], myfit.result['Uncertainties'][2], S[2])
msg += "\nDisc Axis Ratio:  %.4f +/- %.4f ; True: %.3f " % \
       (myfit.result['Parameters'][3], myfit.result['Uncertainties'][3], S[3])
msg += "\nDisc Pos. Ang (deg.):  %.4f +/- %.4f ; True: %.3f " % \
       (myfit.result['Parameters'][4], myfit.result['Uncertainties'][4], S[4])

tac = time.time()
msg += "\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n" % (tac-tic)
with open('test1.dat', 'w') as resf:
    print("\n\n\nTEST 1: ELLIPTICITY\n", file=resf)
    print(msg, file=resf)
