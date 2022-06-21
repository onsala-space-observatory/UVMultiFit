import sys
import uvmultimodel as uvmod
from NordicARC import uvmultifit as uvm

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

foo = uvm.uvmultifit()
foo.select_data(vis)
foo.select_model(model, var=modvars, p_ini=initial, bounds=modbounds)

foo.mymodel = uvm.modeler()
foo.mymodel.Ccompmodel = uvmod.modelcomp

if True:
    if not foo.start_fit():
        print("failed to fit model")
        sys.exit(1)
    else:
        print("success!")
else:
    if not foo.checkInputs():
        print("failed to check inputs")
        sys.exit(1)

    if not foo.readData(del_data=False):
        print("failed to read data")
        sys.exit(1)

    if not foo.initData():
        print("failed to init data")
        sys.exit(1)

    if not foo.initModel():
        print("failed to init model")
        sys.exit(1)

    if not foo.fit():
        print("failed to fit model")
    else:
        print("success!")

# exact solution
S = [2.435e-01, 1.000e+00, 2.000e-01, 5.000e-01, 6.000e+01]

r = foo.result
print(f"disc flux at 50GHz (Jy): {r['Parameters'][0]:7.3f} +/- {r['Uncertainties'][0]:.3f}, true: {S[0]:7.3f}")
print(f"disc spectral index:     {r['Parameters'][1]:7.3f} +/- {r['Uncertainties'][1]:.3f}, true: {S[1]:7.3f}")
print(f"disc size (as):          {r['Parameters'][2]:7.3f} +/- {r['Uncertainties'][2]:.3f}, true: {S[2]:7.3f}")
print(f"disc axis ratio:         {r['Parameters'][3]:7.3f} +/- {r['Uncertainties'][3]:.3f}, true: {S[3]:7.3f}")
print(f"disc pos.angle (deg):    {r['Parameters'][4]:7.3f} +/- {r['Uncertainties'][4]:.3f}, true: {S[4]:7.3f}")
