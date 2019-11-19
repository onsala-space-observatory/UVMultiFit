import unittest
import imp
# import numpy as np

uvm = imp.load_dynamic("_uvmultimodel", "_uvmultimodel.so")

class TestUVMultiModel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):           # runs before all tests
        pass

    @classmethod
    def tearDownClass(cls):        # runs after all tests
        pass

    def setUp(self):               # runs before every single test
        pass

    def tearDown(self):            # runs after every single test
        pass

    def test_clearPointers(self):
        # test clearPointers with argument 2
        result = uvm.clearPointers(2)
        # we expect 0 back
        self.assertEqual(result, 0)

    def test_clearData(self):
        # we can clear data by calling clearPointers with argument 0
        result = uvm.clearPointers(0)
        # we expect 0 back
        self.assertEqual(result, 0)

    def test_clearModel(self):
        # we can clear model by calling clearPointers with argument 1
        result = uvm.clearPointers(1)
        # we expect 0 back
        self.assertEqual(result, 0)

    def test_setNspw(self):
        # function takes one integer
        result = uvm.setNspw(10)
        # we expect 0 back
        self.assertEqual(result, 0)

#     def test_modelcomp(self):
#         IF = 10
#         nui = 100
#         mode = 0
#         result = uvm.modelcomp(IF, nui, mode)

#     def test_setData(self):
#         IF = 1
#
#         pu = np.zeros(10, dtype=np.double)
#         pv = np.zeros(10, dtype=np.double)
#         pw = np.zeros(10, dtype=np.double)
#
#         pwgt = np.ones(10, dtype=np.double)
#         preal = np.zeros(10, dtype=np.complex)
#         poreal = np.zeros(10, dtype=np.complex)
#
#         pfreqs  = np.zeros(10, dtype=np.double)
#         pfittable = np.zeros(10, dtype=np.int8)
#         pwgtcorr  = np.zeros(10, dtype=np.double)
#         dtime  = np.zeros(10, dtype=np.double)
#         tArr = np.zeros(10, dtype=np.double)
#         tIdx = np.zeros(10, dtype=np.int32)
#
#         RAoffset  = np.zeros(10, dtype=np.double)
#         Decoffset  = np.zeros(10, dtype=np.double)
#         Stretchoff  = np.zeros(10, dtype=np.double)
#
#         ant1l = np.zeros(10, dtype=np.int32)
#         ant2l = np.zeros(10, dtype=np.int32)
#         iG = np.zeros(10, dtype=np.int8)
#
#         Nants = 10
#
#         result = uvm.setData(IF, pu, pv, pw, pwgt, preal, poreal, pfreqs, pfittable, pwgtcorr,
#                              dtime, tArr, tIdx, RAoffset, Decoffset, Stretchoff,
#                              ant1l, ant2l, iG, Nants)
#
#         # ASSERT
#         # we expect 0 back
#         self.assertEqual(result, 10)

#     def test_setNCPU(self):
#         # function takes one integer
#         i = 4
#         result = uvm.setNCPU(i)
#         # we expect 0 back
#         self.assertEqual(result, 0)

if __name__ == "__main__":
    unittest.main()
