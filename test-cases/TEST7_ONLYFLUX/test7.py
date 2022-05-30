import time
from NordicARC import uvmultifit as uvm

tic = time.time()
myfit = uvm.uvmultifit(vis=visname, spw='0',
                       model=modelshape, OneFitPerChannel=False,
                       var=modvars, write='residuals',
                       p_ini=pini, only_flux=True,
                       bounds=parbound)
msg = 'NOTICE THAT, USING ONLY-FLUX, THE ESTIMATED UNCERTAINTIES\nARE UNREALISTICALLY SMALL'
msg += 'NOTICE ALSO THAT FIXING THE SOURCE\nPOSITIONS LEAVES SOME SMALL SYSTEMATIC RESIDUALS (DUE TO A\nSUB-PIXEL SHIFT RELATED TO SIMOBSERVE).\n\n '

for pi in range(len(pini)):
    msg += '\n Parameter %i: %.4f +/- %.3e | True value %.2f' % (
        pi, myfit.result['Parameters'][pi], myfit.result['Uncertainties'][pi], S[pi])

tac = time.time()
msg += '\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n' % (tac-tic)
resf = open('test7.dat', 'w')
print >> resf, '\n\n\nTEST 7: ONLY-FLUX\n'
print >> resf, msg
resf.close()
