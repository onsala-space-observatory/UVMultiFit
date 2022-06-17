import pytest

import uvmultimodel as uvmod
from NordicARC import uvmultifit as uvm

@pytest.fixture
def uwm():
    return uvm.uvmultifit()

@pytest.fixture
def modeler():
    mymodel = uvm.modeler()
    mymodel.Ccompmodel = uvmod.modelcomp
    return mymodel

@pytest.fixture
def vis():
    return 'casatest/Disc/Disc.alma.out10.noisy.ms'

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

def test_uvmultifit_init(uwm):
    assert uwm.NCPU == 4

def test_uvmultifit_select_data(uwm, vis):
    uwm.select_data(vis)
    assert uwm.column == 'data' and uwm.field == 0

def test_uvmultifit_select_model(uwm, model):
    uwm.select_model(model)
    assert uwm.var == ['p[0], p[1], p[2]']
    assert len(uwm.model) == 1 and uwm.model[0] == model[0]

def test_checkInputs(uwm, vis, model, modvars, initial, modeler):
    uwm.select_data(vis)
    uwm.select_model(model, modvars, initial)
    uwm.mymodel = modeler
    assert uwm.checkInputs() is True

def test_readData(uwm, vis, model, modvars, initial, modeler):
    uwm.select_data(vis)
    uwm.select_model(model, modvars, initial)
    uwm.mymodel = modeler
    uwm.checkInputs()
    assert uwm.readData(del_data=False) is True

def test_initData(uwm, vis, model, modvars, initial, modeler):
    uwm.select_data(vis)
    uwm.select_model(model, modvars, initial)
    uwm.mymodel = modeler
    uwm.checkInputs()
    uwm.readData(del_data=False)
    assert uwm.initData() is True

def test_initModel(uwm, vis, model, modvars, initial, modeler):
    uwm.select_data(vis)
    uwm.select_model(model, modvars, initial)
    uwm.mymodel = modeler
    uwm.checkInputs()
    uwm.readData(del_data=False)
    uwm.initData()
    assert uwm.initModel() is True

def test_fitModel(uwm, vis, model, modvars, initial, modeler):
    uwm.select_data(vis)
    uwm.select_model(model, modvars, initial)
    uwm.mymodel = modeler
    uwm.checkInputs()
    uwm.readData(del_data=False)
    uwm.initData()
    uwm.initModel()
    assert uwm.fit() is True

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
