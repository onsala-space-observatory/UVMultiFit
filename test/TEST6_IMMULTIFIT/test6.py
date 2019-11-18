import time

tic = time.time()

myfit = uvm.immultifit(residual=immname, psf=psfname,
                       model=modelshape, OneFitPerChannel=False,
                       var=modvars, write='residuals',
                       p_ini=pini, LMtune=[1.e-3, 10., 1.e-8, 200, 1.e-5],
                       bounds=parbound)
msg = ''
for pi in range(len(pini)):
    msg += '\n Parameter %i: %.4f +/- %.4f | True value %.2f' % (
        pi, myfit.result['Parameters'][pi], myfit.result['Uncertainties'][pi], S[pi])

tac = time.time()
msg += '\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n' % (tac-tic)
resf = open('test6.dat', 'w')
print >> resf, '\n\n\nTEST 6: IMMULTIFIT\n'
print >> resf, msg
resf.close()
