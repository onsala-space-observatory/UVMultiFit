import pytest
import numpy as np

import uvmultimodel as uvmod
from NordicARC.measurementset import MeasurementSet
from NordicARC.modeler import Modeler
from NordicARC.utils import is_list_of_floats

# @pytest.fixture
# def uwm():
#     return uvm.uvmultifit(vis='')

# @pytest.fixture
# def mdl():
#     return Modeler()
#
# @pytest.fixture
# def ms():
#     return MeasurementSet()

@pytest.fixture
def vis():
    return 'Disc/Disc.alma.out10.noisy.ms'

@pytest.fixture
def model():
    return ['disc']

@pytest.fixture
def modvars():
    Nu = '50.0GHz'
    NuF = float(Nu.split('G')[0])*1.e9
    return "0,0, p[0]*(nu/%.4e)**p[1], p[2], p[3], p[4]" % NuF

@pytest.fixture
def initial():
    Nu = '50.0GHz'
    si = ['disc', 1.0, '0.2arcsec', '0.1arcsec', '60.0deg',
          Nu, 1.0, "J2000 10h00m00.0s -30d00m00.0s"]
    size = float(si[2].split('a')[0])
    minor = float(si[3].split('a')[0])
    return [0.8, 0.0, size*1.2, minor/size*0.8, 45.0]

@pytest.fixture
def bounds():
    Nu = '50.0GHz'
    si = ['disc', 1.0, '0.2arcsec', '0.1arcsec', '60.0deg',
          Nu, 1.0, "J2000 10h00m00.0s -30d00m00.0s"]
    size = float(si[2].split('a')[0])
    minor = float(si[3].split('a')[0])
    bounds = []
    bounds.append([0.0, 1.0])
    bounds.append([-1.0, 1.0])
    bounds.append([size*1.0, size*1.5])
    bounds.append([minor/size*0.5, minor/size*1.0])
    bounds.append([40.0, 50.0])
    return bounds

def test_clearPointers():
    # test clearPointers with argument 2
    result = uvmod.clearPointers(2)
    # we expect 2 back
    assert result == 2

def test_setNspw():
    # function takes one integer
    result = uvmod.setNspw(10)
    # we expect 10 back
    assert result == 10

def test_setNCPU():
    # function takes one integer
    result = uvmod.setNCPU(4)
    # we expect 4 back
    assert result == 4

def test_measurementset_init(vis):
    ms = MeasurementSet(vis)
    assert ms.check_measurementset()

def test_checkInputs(model, modvars, initial):
    mdl = Modeler(model=model, var=modvars, p_ini=initial)
    indices = mdl.check_model_consistency()
    assert isinstance(indices, list)

def test_checkBounds(model, modvars, initial, bounds):
    mdl = Modeler(model=model, var=modvars, p_ini=initial, bounds=bounds)
    bounds = mdl.check_bounds(bounds, initial)
    print(bounds)
    assert all(is_list_of_floats(b, 2) for b in bounds)

def test_fitRange_cont(model, modvars, initial):
    mdl = Modeler(model=model, var=modvars, p_ini=initial)
    nui, spwrange = mdl.fit_range(-1, -1, [1.0, 2.0])
    assert nui == -1 and spwrange == [0, 1]

def test_fitRange_spec(model, modvars, initial):
    mdl = Modeler(model=model, var=modvars, p_ini=initial)
    nui, spwrange = mdl.fit_range(1, 1, [1.0, 2.0])
    assert nui == 1 and spwrange == [1]

def test_phase_center(vis):
    phase_center = 'J2000 12h34m56.0s 01d02m03.0s'
    ms = MeasurementSet(vis)
    refpos = np.array(ms.check_phase_center(phase_center))
    delta = refpos - np.array([3.29401807, 0.01804961])
    assert (np.abs(delta) < 1.0e-6).all()

# def test_initModel(uwm, vis, model, modvars, initial, modeler):
#     uwm.select_data(vis)
#     uwm.select_model(model, modvars, initial)
#     uwm.mymodel = modeler
#     uwm.checkInputs()
#     uwm.readData(del_data=False)
#     uwm.initData()
#     assert uwm.initModel() is True
#
# def test_fitModel(uwm, vis, model, modvars, initial, modeler):
#     uwm.select_data(vis)
#     uwm.select_model(model, modvars, initial)
#     uwm.mymodel = modeler
#     uwm.checkInputs()
#     uwm.readData(del_data=False)
#     uwm.initData()
#     uwm.initModel()
#     assert uwm.fit() is True

### def test_modeler(uwm, vis, model, modeler):
###     uwm.select_data(vis)
###     uwm.select_model(model)
###     uwm.freqs = [[0.0]]
###
###     for ii, component in enumerate(uwm.model):
###         tempstr = uwm.var[ii].replace(
###             'LorentzLine(', 'self.LorentLine(nu,').replace(
###                 'GaussLine(', 'self.GaussLine(nu, ').replace(
###                     'nu0', '%.12f' % uwm.freqs[0][0])
###
###     modeler._setup(uwm.model, uwm.var, uwm.fixed, uwm.fixedvar, uwm.scalefix, uwm.NCPU, uwm.only_flux,
###                    0, uwm.isNumerical, set(),
###                    [{}, {}, {}, {}, {}, {}], False)
###     modeler._compileAllModels()
###     assert 1 == 1
