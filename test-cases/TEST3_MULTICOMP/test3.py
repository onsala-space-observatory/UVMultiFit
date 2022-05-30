import time
from NordicARC import uvmultifit as uvm

tic = time.time()
myfit = uvm.uvmultifit(vis=visname, spw='0',
                       model=modelshape, OneFitPerChannel=False,
                       var=modvars, write='residuals',
                       p_ini=pini,
                       bounds=parbound)
msg = ''
for pi in range(len(pini)):
    msg += '\n Parameter %i: %.4f +/- %.4f | True value %.2f' % (
        pi, myfit.result['Parameters'][pi], myfit.result['Uncertainties'][pi], S[pi])

tac = time.time()
msg += '\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n' % (tac-tic)
with open('test3.dat', 'w') as resf:
    print('\n\n\nTEST 3: MULTI-COMPONENT WITH CORRELATED VARIABLES\n', file=resf)
    print(msg, file=resf)
