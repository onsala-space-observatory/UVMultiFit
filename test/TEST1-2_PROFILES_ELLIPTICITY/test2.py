import numpy as np
import pylab as pl
import time

tic = time.time()

rf = open('test2.%s.dat' % modshape, 'w')

print >> rf, '\n\n\nTEST 2: RADIAL PROFILES. SHAPE: %s\n' % modshape

myfit = uvm.uvmultifit(vis=visname, write='residuals',
                       spw='0', model=modshape, OneFitPerChannel=True,
                       var=shapevar, p_ini=pi, bounds=parbound, HankelOrder=80)

maxdev = 100.*np.max(np.abs(myfit.result['Parameters'][0][:, 0] - 1.0)/1.0)
msg = '\n Model %s. Maximum Flux deviation: %.4f %%' % (modshape, maxdev)
ModD = Diameter + np.linspace(0, DiamDelta,
                              np.shape(myfit.result['Frequency'][0])[0])
maxdev = 100.*np.max(np.abs(myfit.result['Parameters'][0][:, 1] - ModD)/ModD)
msg += '\n Model %s. Maximum Size deviation: %.4f %%' % (modshape, maxdev)

fig = pl.figure()
sub = fig.add_subplot(211)
sub2 = fig.add_subplot(212, sharex=sub)
fig.subplots_adjust(hspace=0.08)

nu0, nu1 = myfit.result['Frequency'][0][0] / \
    1.e9, myfit.result['Frequency'][0][-1]/1.e9
sub.cla()
sub2.cla()
sub.plot(myfit.result['Frequency'][0]/1.e9, myfit.result['Parameters'][0][:, 0], 'or', label='Fit')
sub.plot(np.array([nu0, nu1]), np.array([1.0, 1.0]), '--b', label='True Source')
sub.set_ylim((np.min(myfit.result['Parameters'][0][:, 0]) * 0.9, np.max(myfit.result['Parameters'][0][:, 0])*1.1))
sub.set_ylabel('Flux Density (Jy)')
pl.setp(sub.get_xticklabels(), visible=False)
sub2.plot(myfit.result['Frequency'][0]/1.e9, myfit.result['Parameters'][0][:, 1], 'or', label='Fit')
sub2.plot(np.array([nu0, nu1]), np.array([Diameter, Diameter+DiamDelta]), '--b', label='True Source')
sub2.set_ylabel('Size (arcsec)')
sub2.set_xlabel('Frequency (GHz)')
sub2.set_ylim((0.0, (Diameter+DiamDelta)*1.1))

if modshape == 'GaussianRing':
    sub2.plot(myfit.result['Frequency'][0]/1.e9, myfit.result['Parameters'][0][:, 2], 'og', label='Fit')
    sub2.plot(np.array([nu0, nu1]), np.array([Sigma, Sigma]), '--k', label='Fit')

    maxdev = 100. * np.max(np.abs(myfit.result['Parameters'][0][:, 2] - Sigma)/Sigma)
    msg += '\n Model %s. Maximum Sigma deviation: %.4f %%' % (modshape, maxdev)

tac = time.time()
msg += '\n\n DATA READ AND FIT LASTED %.2f SECONDS.\n' % (tac-tic)
print >> rf, msg
sub.legend(numpoints=1)
pl.savefig('TEST2.%s.modelfit.png' % modshape)
rf.close()
