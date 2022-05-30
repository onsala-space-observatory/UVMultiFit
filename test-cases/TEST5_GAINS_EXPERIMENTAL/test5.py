import time
from NordicARC import uvmultifit as uvm

tic = time.time()
myfit = uvm.uvmultifit(vis=visname, spw='0', stokes=STK,
                       model=modelshape, OneFitPerChannel=False,
                       var=modvars, write='residuals', LMtune=LMtune,
                       p_ini=pini, amp_gains=amp_gains,
                       bounds=parbound)
msg = ''
for pi in range(len(pini)):
    msg += '\n Parameter %i: %.4f +/- %.4f | True value %.2f' % (
        pi, myfit.result['Parameters'][pi], myfit.result['Uncertainties'][pi], S[pi])

tac = time.time()
msg += '\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n' % (tac-tic)
with open('test5.dat', 'w') as resf:
    print('\n\n\nTEST 5: SIMULTANEOUS GAIN CALIBRATION \n(TIME-VARYING AMP. GAIN, USING A PIECE-WISE TIME FUNCTION)\n', file=resf)
    print(msg, file=resf)
