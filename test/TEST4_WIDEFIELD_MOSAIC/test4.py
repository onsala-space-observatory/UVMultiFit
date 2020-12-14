import time
from NordicARC import uvmultifit as uvm

tic = time.time()

myfit = uvm.uvmultifit(vis=msout, model=model, var=modvars, field=field,
                       p_ini=pini, NCPU=4, column='data', write='residuals',
                       OneFitPerChannel=False, phase_center=phref, pbeam=True,
                       dish_diameter=25.0, ldfac=1.0, LMtune=[1.e-3, 10., 1.e-5, 10], bounds=bounds)

msg = ''
for pi in range(len(pini)/3):
    msg += '\n Flux of delta %i: %.4f +/- %.4f Jy | True value: %.2f Jy' % (
        pi, myfit.result['Parameters'][3*pi+2], myfit.result['Uncertainties'][3*pi+2], S[3*pi+2])

msg += '\n\n OFFSETS (TAKE INTO ACCOUNT THE COARSE GRIDDING IN SIMOBSERVE!)\n:'
for pi in range(len(pini)/3):
    msg += '\n RA offset delta %i: %.4f +/- %.4f as | True value: %.2f as' % (
        pi, myfit.result['Parameters'][3*pi], myfit.result['Uncertainties'][3*pi], S[3*pi])
    msg += '\n Dec offset delta %i: %.4f +/- %.4f as | True value: %.2f as' % (
        pi, myfit.result['Parameters'][3*pi+1], myfit.result['Uncertainties'][3*pi+1], S[3*pi+1])

    # msg += '\nGlobal offset RA: %.4f +/- %.4f as | True value: %.2f as' % (
    # myfit.result['Parameters'][-2],myfit.result['Uncertainties'][-2],S[-2])
    # msg += '\nGlobal offset Dec: %.4f +/- %.4f as | True value: %.2f as' % (
    # myfit.result['Parameters'][-1],myfit.result['Uncertainties'][-1],S[-1])

tac = time.time()
msg += '\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n' % (tac-tic)
resf = open('test4.dat', 'w')
print >> resf, '\n\n\nTEST 4: WIDE FIELD AND MOSAIC.\n'
print >> resf, msg
resf.close()
