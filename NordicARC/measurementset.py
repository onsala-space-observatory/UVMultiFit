import os
import re
import time
import logging
from contextlib import contextmanager

from typing import List

import numpy as np

from casatools import ms             # type: ignore
from casatools import table          # type: ignore
from casatools import coordsys       # type: ignore

from .utils import get_list_of_strings, is_list_of_int, is_list_of_floats
from .utils import is_casa_position, is_valid_stokes

class FieldData():

    def __init__(self, datascanAv, weightscan, uscan, vscan, wscan, tscan, tArray, tIndex,
                 RAscan, Decscan, Stretchscan, ant1scan, ant2scan, modelscanAv=None):
        self.datascanAv = datascanAv
        self.weightscan = weightscan
        self.uscan = uscan
        self.vscan = vscan
        self.wscan = wscan
        self.tscan = tscan
        self.tArray = tArray
        self.tIndex = tIndex
        self.RAscan = RAscan
        self.Decscan = Decscan
        self.Stretchscan = Stretchscan
        self.ant1scan = ant1scan
        self.ant2scan = ant2scan
        self.modelscanAv = modelscanAv

    def __repr__(self):
        return f"FieldData(tIndex = {self.tIndex[0]}...{self.tIndex[-1]})"

@contextmanager
def open_ms(msname: str):
    """A context manager to reattach to a measurement set table."""

    mstool = ms()
    try:
        logging.debug(f"open measurementset {msname}")
        mstool.open(msname)
        yield mstool
    finally:
        logging.debug(f"close measurementset {msname}")
        mstool.close()

@contextmanager
def open_tb(msname: str):
    """A context manager to open a disk file containing an existing casa table."""

    tbtool = table()
    try:
        logging.debug(f"open measurementset {msname}")
        tbtool.open(msname)
        yield tbtool
    finally:
        logging.debug(f"close measurementset {msname}")
        tbtool.close()


class MeasurementSet():
    """Class to deal with model equations and fitting procedures.

    It convert strings representing models and parameters into
    compiled equations, to be used in a ChiSq visibility fitting. It
    also interacts with the C++ extension of UVMultiFit.

    This class should NOT be instantiated by the user (it is called
    from the UVMultiFit class.

    """

    logger = logging.getLogger("measurement")

    def __init__(self, vis: str, spw: str = '0', field: int = 0, scans: List = [],
                 corrected: bool = False, uniform: bool = False, uvtaper: float = 0.0,
                 chanwidth: int = 1, timewidth: int = 1, stokes: str = 'I',
                 MJDrange: List[float] = [-1.0, -1.0], ldfac: float = 1.22,
                 phase_center: str = '', pbeam: bool = False, wgt_power: float = 1.0,
                 dish_diameter: float = 0.0) -> None:
        self.vis = get_list_of_strings(vis)
        self.spw = get_list_of_strings(spw)
        if len(self.spw) > 1 and len(self.vis) != len(self.spw):
            self.logger.error("the length of 'spw' is not equal to the length of 'vis'!")

        self.field = 0
        if isinstance(field, (int, str)):
            self.field = field
        else:
            self.logger.warning("parameter 'field' needs to be an integer or string, resetting to 0!")
            self.field = 0
        if isinstance(self.field, str) and self.field.isdigit():
            self.field = int(self.field)
        self.field_id = []                       # type: List
        self.pointing = []                       # type: List
        self.sourscans = []                      # type: List
        if is_list_of_int(scans):
            self.scans = scans
        else:
            self.logger.warning("'scans' should be list of (lists of) integers, resetting to empty list")
            self.scans = []
        if isinstance(self.scans, list) and len(self.scans) == 0:
            self.scans = [[] for v in self.vis]
        self.column = 'corrected_data' if corrected else 'data'
        self.uniform = uniform
        self.uvtaper = uvtaper
        self.wgt_power = wgt_power
        self.dish_diameter = dish_diameter
        self.chanwidth = chanwidth
        if timewidth != 1:
            self.logger.warning("timewdith>1 cannot be currently set, due to issues with new MS tool")
        self.timewidth = 1   # timewidth
        if is_valid_stokes(stokes):
            self.stokes = stokes
        else:
            self.logger.warning("unknown polarization product 'stokes', resetting to 'I'")
            self.stokes = 'I'
        if is_list_of_floats(MJDrange):
            self.MJDrange = MJDrange
        else:
            self.logger.warning("invalid MJD range, resetting to full range")
            self.MJDrange = [-1.0, -1.0]
        self.Nants = 0
        self.refpos = self.check_phase_center(phase_center)
        self.ldfac = ldfac
        self.pbeam = pbeam
        self.spwlist = []                        # type: List
        self.pol2aver = []                       # type: List
        self.polmod = []                         # type: List
        self.polii = []                          # type: List

    def dump(self):
        temp = vars(self)
        print(f"{self.__class__.__name__}(")
        for item in sorted(temp.keys()):
            print(f"  {item}: {temp[item]}")
        print(")")

    def __repr__(self):
        txt = "MeasurementSet("
        txt += f"vis = '{self.vis}', "
        txt += f"spw = {self.spw}, "
        txt += f"column = '{self.column}', "
        txt += f"field = {self.field}, "
        txt += f"scans = {self.scans}, "
        txt += f"uniform = {self.uniform}, "
        txt += f"uvtaper = {self.uvtaper}, "
        txt += f"chanwidth = {self.chanwidth}, "
        txt += f"timewidth = {self.timewidth}, "
        txt += f"stokes = '{self.stokes}', "
        txt += f"MJDrange = {self.MJDrange})"
        return txt

    def check_phase_center(self, phase_center):
        self.logger.debug("MeasurementSet::check_phase_center")
        if not isinstance(phase_center, str):
            self.logger.error("'phase_center' must be a string!")
            return None
        if len(phase_center) == 0:
            return None

        if not is_casa_position(phase_center):
            self.logger.warning("'phase_center' is not a CASA-formatted sky coordinate!")
            return None

        csys = coordsys().newcoordsys(direction=True)
        dirstr = phase_center.split()
        if len(dirstr) == 2:
            csys.setdirection(refcode="J2000", refval=phase_center)
        else:
            csys.setdirection(refcode=dirstr[0], refval=" ".join(dirstr[1:]))
            csys.convertdirection("J2000")
        refpos = np.copy(csys.torecord()['direction0']['crval'])
        return refpos

    def get_field_index(self, vi, msname):
        field_index = -1
        phasedir = None
        with open_ms(msname) as ms:
            allfields = list(ms.range('fields')['fields'])
            if isinstance(self.field, int):
                field_index = self.field
            elif is_list_of_int(self.field):
                field_index = int(self.field[vi])
            else:
                aux = str(self.field)
                self.field = [aux for v in self.vis]
                field_found = False
                for f, field in enumerate(allfields):
                    if self.field[vi] in field:
                        field_found = True
                        field_index = f
                        break
                    if not field_found:
                        self.logger.error(f"field '{self.field[vi]}' is not in '{msname}'")
                        return field_index, None
            phasedir = ms.range('phase_dir')['phase_dir']['direction'][:, field_index]
        return field_index, phasedir

    def check_measurementset(self):
        self.logger.debug("MeasurementSet::check_measurementset")
        for v in self.vis:
            if not os.path.exists(v):
                self.logger.error(f"measurement set '{v}' does not exist!")
                return False
        phasedirs = [{} for v in self.vis]

        # Open MS and look for the selected data:
        for vi, v in enumerate(self.vis):
            self.field_id.append([])
            field_index, phasedir = self.get_field_index(vi, v)
            if field_index != -1:
                self.field_id[-1].append(field_index)
                phasedirs[vi][field_index] = phasedir
            else:
                return False
        if self.refpos is None:
            self.refpos = phasedirs[0][min(phasedirs[0].keys())]
        success = (self.find_observing_positions(phasedirs)
                   and self.get_spectral_configuration()
                   and self.get_polarization_configuration()
                   and self.set_weight_equation(self.pbeam))
        self.phasedirs = phasedirs
        return success

    def find_observing_positions(self, phasedirs):
        self.logger.debug("MeasurementSet::find_observing_positions")
        for vi, v in enumerate(self.vis):
            with open_ms(v) as ms:
                info = ms.getscansummary()

            self.sourscans.append([])
            self.pointing.append([])
            for key in info.keys():
                if info[key]['0']['FieldId'] in self.field_id[vi]:
                    self.sourscans[-1].append(int(key))
                for fieldid in info[key].keys():
                    myfield = info[key][fieldid]['FieldId']
                    if myfield in self.field_id[vi]:
                        # fi = self.field_id[vi].index(myfield)
                        self.pointing[-1].append(phasedirs[vi][myfield])

            self.pointing[-1] = np.array(self.pointing[-1])

        for vi, v in enumerate(self.vis):
            if len(self.scans[vi]) > 0:
                goodscid = [x for x in self.scans[vi] if x in self.sourscans[vi]]
                if len(goodscid) != len(self.scans[vi]):
                    badscid = [x for x in self.scans[vi] if x not in self.sourscans[vi]]
                    msg = f"the following scans do NOT correspond to source {str(self.field)} "
                    msg += str(badscid)
                    self.logger.error(msg)
                    return False
                self.sourscans[vi] = list(goodscid)
        return True

    def get_spectral_configuration(self):
        self.logger.debug("MeasurementSet::get_spectral_configuration")
        spwi = 0
        for vi, v in enumerate(self.vis):
            j = vi if len(self.spw) != 1 else 0

            with open_ms(v) as ms:
                # freqdic = ms.getspectralwindowinfo()
                spwchans = ms.range(['num_chan'])['num_chan']

            aux = MeasurementSet.channeler(self.spw[j], width=self.chanwidth, maxchans=spwchans)
            if aux[0]:
                ranges = list(aux[1])
            else:
                self.logger.error(aux[1] + f"something seems to be wrong with the 'spw' number {vi}.")
                return False

            nspws = range(len(spwchans))
            selspws = [[i, ranges[i]] for i in nspws if len(ranges[i]) > 0]
            self.spwlist.append([j, vi, selspws, spwi])
            # spwlist[vi][2] es una lista con [spwid, chanranges]

        return True

    def get_polarization_configuration(self):
        self.logger.debug("MeasurementSet::get_polarization_configuration")
        for v in self.vis:
            with open_ms(v) as ms:
                polprods = [x[0] for x in list(ms.range(['corr_names'])['corr_names'])]

            self.pol2aver.append(np.zeros(len(polprods)))
            # 0: normal,   1: multiply by i,   2: pol. independent, 3: just one product
            self.polmod.append(0)
            if self.stokes not in ['I', 'Q', 'U', 'V'] + polprods:
                self.logger.error(f"stokes parameter '{self.stokes}' not available")
                return False
            if self.stokes in polprods:
                self.pol2aver[-1][polprods.index(self.stokes)] = 1.0
                self.polii.append([polprods.index(self.stokes)])
                self.polmod[-1] = 3
            # User asks for a Stokes parameter:
            # CASE 1: Circular feeds.
            elif 'RR' in polprods:
                if self.stokes == 'I':
                    self.polii.append([polprods.index('RR'), polprods.index('LL')])
                    self.pol2aver[-1][polprods.index('RR')] = 0.5
                    self.pol2aver[-1][polprods.index('LL')] = 0.5
                if self.stokes == 'PI':
                    self.polii.append([polprods.index('RR'), polprods.index('LL')])
                    self.pol2aver[-1][polprods.index('RR')] = 0.5
                    self.pol2aver[-1][polprods.index('LL')] = 0.5
                    self.polmod[-1] = 2
                if self.stokes == 'Q':
                    self.polii.append([polprods.index('RL'), polprods.index('LR')])
                    self.pol2aver[-1][polprods.index('RL')] = 0.5
                    self.pol2aver[-1][polprods.index('LR')] = 0.5
                if self.stokes == 'U':
                    self.polii.append([polprods.index('RL'), polprods.index('LR')])
                    self.pol2aver[-1][polprods.index('RL')] = 0.5
                    self.pol2aver[-1][polprods.index('LR')] = 0.5
                    self.polmod[-1] = 1
                if self.stokes == 'V':
                    self.polii.append([polprods.index('RR'), polprods.index('LL')])
                    self.pol2aver[-1][polprods.index('RR')] = 0.5
                    self.pol2aver[-1][polprods.index('LL')] = -0.5
            #  CASE 2: Linear feeds.
            elif 'XX' in polprods:
                if self.stokes == 'I':
                    self.polii.append([polprods.index('XX'), polprods.index('YY')])
                    self.pol2aver[-1][polprods.index('XX')] = 0.5
                    self.pol2aver[-1][polprods.index('YY')] = 0.5
                if self.stokes == 'PI':
                    self.polii.append([polprods.index('XX'), polprods.index('YY')])
                    self.pol2aver[-1][polprods.index('XX')] = 0.5
                    self.pol2aver[-1][polprods.index('YY')] = 0.5
                    self.polmod[-1] = 2
                if self.stokes == 'Q':
                    self.polii.append([polprods.index('XX'), polprods.index('YY')])
                    self.pol2aver[-1][polprods.index('XX')] = 0.5
                    self.pol2aver[-1][polprods.index('YY')] = -0.5
                if self.stokes == 'U':
                    self.polii.append([polprods.index('XY'), polprods.index('YX')])
                    self.pol2aver[-1][polprods.index('XY')] = 0.5
                    self.pol2aver[-1][polprods.index('YX')] = 0.5
                if self.stokes == 'V':
                    self.polii.append([polprods.index('YX'), polprods.index('XY')])
                    self.pol2aver[-1][polprods.index('YX')] = 0.5
                    self.pol2aver[-1][polprods.index('XY')] = -0.5
                    self.polmod[-1] = 1
            else:
                self.logger.error(f"polarization '{self.stokes}' not understood")
                return False

        return True

    @staticmethod
    def get_cross_correlations(msname):
        with open_tb(msname) as tb:
            SPW = tb.getcol('DATA_DESC_ID')
            crosscorr = tb.getcol('ANTENNA1') != tb.getcol('ANTENNA2')
        return SPW, crosscorr

    @staticmethod
    def get_data_description(msname):
        subdir = os.path.join(msname, 'DATA_DESCRIPTION')
        with open_tb(subdir) as tb:
            DDSC = tb.getcol('SPECTRAL_WINDOW_ID')
        return DDSC

    @staticmethod
    def get_spectral_window(msname):
        subdir = os.path.join(msname, 'SPECTRAL_WINDOW')
        with open_tb(subdir) as tb:
            origfreqs = tb.getcol('CHAN_FREQ')
        return origfreqs

    @staticmethod
    def get_scan_mask(msname, scan, maskspw):
        with open_tb(msname) as tb:
            masksc = maskspw * (tb.getcol('SCAN_NUMBER') == int(scan))
            fieldids = list(np.sort(np.unique(tb.getcol('FIELD_ID')[masksc])))
        return masksc, fieldids

    @staticmethod
    def get_field_id(msname, masksc, fieldid, takeModel, column):
        with open_tb(msname) as tb:
            maskfld = np.where(masksc * (tb.getcol('FIELD_ID') == int(fieldid)))[0]
            uvscan = None
            times = None
            datascan = None
            if len(maskfld) != 0:
                tb2 = tb.selectrows(maskfld)
                uvscan = {'uvw': tb2.getcol('UVW'),
                          'antenna1': tb2.getcol('ANTENNA1'),
                          'antenna2': tb2.getcol('ANTENNA2'),
                          'time': tb2.getcol('TIME')}
                times = np.unique(uvscan['time'])
                if takeModel:
                    datascan = {column: tb2.getcol((column).upper()),
                                'model_data': tb2.getcol('MODEL_DATA'),
                                'weight': tb2.getcol('WEIGHT'),
                                'flag': tb2.getcol('FLAG')}
                else:
                    datascan = {column: tb2.getcol((column).upper()),
                                'weight': tb2.getcol('WEIGHT'),
                                'flag': tb2.getcol('FLAG')}

        return maskfld, uvscan, times, datascan

    def get_field_data(self, fieldid, scan, masksc, msname, vis_index, rang,
                       sp, si, max_tIndex, i0scan, takeModel):
        self.logger.info(f"reading scan #{scan}, field: {fieldid}")
        F = self.get_field_id(msname, masksc, fieldid, takeModel, self.column)
        maskfld, uvscan, times, datascan = F
        if len(maskfld) == 0:
            return None

        # Compute the polarization product:
        # Bad data has zero weight:
        datascan['weight'][np.logical_not(np.isfinite(datascan['weight']))] = 0.0

        nfreq = len(rang)
        # All unflagged weights set to equal (if uniform):
        if self.uniform:
            datascan['weight'][datascan['weight'] > 0.0] = 1.0

        datascan['weight'][datascan['weight'] < 0.0] = 0.0

        copyweight = np.copy(datascan['weight'])
        totalmask = datascan['flag']
        origmasked = np.ma.array(datascan[self.column], mask=totalmask, dtype=np.complex128)

        if takeModel:
            origmodmasked = np.ma.array(datascan['model_data'],
                                        mask=totalmask,
                                        dtype=np.complex128)

        # The weights are weighting the RESIDUALS, and not the ChiSq terms.
        # Hence, we divide wgt_power by 2.:
        origweight = np.power(copyweight, self.wgt_power / 2.)

        # Completely flagged times/baselines:
        origweight[np.sum(np.logical_not(totalmask), axis=1) == 0] = 0.0

        datamask = 0.0
        weightmask = 0.0
        flagmask = 0.0
        modelmask = 0.0

        # Construct the required polarization from the correlation products:
        polavg = [pol != 0.0 for pol in self.pol2aver[vis_index]]

        if self.polmod[vis_index] == 2:
            flagmask = np.ma.logical_and(
                totalmask[self.polii[vis_index][0], :, :],
                totalmask[self.polii[vis_index][1], :, :])
            datamask = np.ma.average(origmasked[polavg, :].real, axis=0) + \
                1.j * np.ma.average(origmasked[polavg, :].imag, axis=0)

            if takeModel:
                modelmask = np.ma.average(origmodmasked[polavg, :].real, axis=0) + \
                    1.j * np.ma.average(origmasked[polavg, :].imag, axis=0)

            weightmask = np.sum(origweight[polavg, :], axis=0)
        else:
            if self.polmod[vis_index] == 3:
                flagmask = totalmask[self.polii[vis_index][0], :, :]
            else:
                flagmask = np.ma.logical_or(totalmask[self.polii[vis_index][0], :, :],
                                            totalmask[self.polii[vis_index][1], :, :])

            for pol in self.polii[vis_index]:
                datamask += origmasked[pol, :] * self.pol2aver[vis_index][pol]
                if takeModel:
                    modelmask += origmodmasked[pol, :] * self.pol2aver[vis_index][pol]
                weightmask += origweight[pol, :]

            if self.polmod[vis_index] == 1:
                datamask *= 1.j
                if takeModel:
                    modelmask *= 1.j

        # Free some memory:
        for key in list(datascan):
            del datascan[key]
        del datascan, origmasked, origweight
        del copyweight, totalmask, flagmask
        if takeModel:
            del origmodmasked

        # Mosaic-related corrections:
        phshift = 3600. * 180. / np.pi * (self.phasedirs[vis_index][fieldid] - self.refpos)
        strcos = np.cos(self.phasedirs[vis_index][fieldid][1])

        if phshift[0] != 0.0 or phshift[1] != 0.0:
            self.logger.info(f"offset: {phshift[0]/15.0:.2e} RA (tsec) "
                             f"{phshift[1]:.2e} Dec (asec)")

        # Average spectral channels:
        _, ndata = np.shape(datamask)

        # ntimes = len(times)
        # ntav = int(max([1, round(float(ntimes) / self.timewidth)]))
        datatemp = np.ma.zeros((nfreq, ndata), dtype=np.complex128)
        if takeModel:
            modeltemp = np.ma.zeros((nfreq, ndata), dtype=np.complex128)
        weighttemp = np.ma.zeros((nfreq, ndata))

        if self.chanwidth == 1:
            concRan = [c[0] for c in rang]
            datatemp[:, :] = datamask[concRan, :]
            if takeModel:
                modeltemp[:, :] = modelmask[concRan, :]
            weighttemp[:, :] = weightmask[np.newaxis, :]
        else:
            for nu in range(nfreq):
                datatemp[nu, :] = np.ma.average(
                    datamask[rang[nu], :].real, axis=0) + \
                    1.j * np.ma.average(datamask[rang[nu], :].imag, axis=0)
                if takeModel:
                    modeltemp[nu, :] = np.ma.average(
                        modelmask[rang[nu], :].real, axis=0) + \
                        1.j * np.ma.average(modelmask[rang[nu], :].imag, axis=0)
                weighttemp[nu, :] = weightmask

        shape = np.shape(uvscan['time'])
        RAoff = np.full(shape, float(phshift[0]))
        Decoff = np.full(shape, float(phshift[1]))
        Stretch = np.full(shape, float(strcos))

        ant1scan = np.copy(uvscan['antenna1'][:])
        ant2scan = np.copy(uvscan['antenna2'][:])
        uscan = np.copy(uvscan['uvw'][0, :])
        vscan = np.copy(uvscan['uvw'][1, :])
        wscan = np.copy(uvscan['uvw'][2, :])
        tscan = np.copy(uvscan['time'][:])
        tArray = np.copy(times)

        tIndex = np.zeros(np.shape(uvscan['time']), dtype=np.int32)
        for tid, tiii in enumerate(times):
            tIndex[uvscan['time'] == tiii] = tid
        tIndex += max_tIndex

        datascanAv = np.transpose(datatemp)
        if self.uniform:
            weightscan = np.transpose(np.ones(np.shape(weighttemp)))
        else:
            weightscan = np.transpose(weighttemp)

        # Useful info for function writeModel() and for pointing correction:
        if scan not in self.iscan[msname][sp].keys():
            self.iscan[msname][sp][scan] = []

        self.iscan[msname][sp][scan].append(
            [int(si), int(i0scan), int(ndata), list(rang), np.copy(maskfld)])
        self.iscancoords[si].append([i0scan, i0scan + ndata, phshift[0], phshift[1]])
        fieldData = FieldData(datascanAv, weightscan, uscan, vscan, wscan, tscan,
                              tArray, tIndex, RAoff, Decoff, Stretch, ant1scan, ant2scan)
        if takeModel:
            fieldData.modelscanAv = np.transpose(modeltemp)
        return fieldData

    def read_data(self, takeModel=False):
        """Reads the data, according to the properties ``vis, column, chanwidth``, etc.

        It then fills in the properties ``averdata, averfreqs, averweights, u, v, w``, etc.

        Each one of these properties is a list with the data (one list item per spectral window/scan).

        A previous successful run of function ``checkInputs()`` is assumed.

        .. note:: Instead of re-reading data from scratch, using the same uvmultifit instance, a
        better approach may be to restart CASA and create a fresh uvmultifit instance with the
        new data, avoiding some memory leakage related to potential hidden references to the
        data in the IPython's *recall* prompt."""

        self.logger.debug("MeasurementSet::read_data")
        tic = time.time()

        # if del_data:  # data_changed:
        #     # self.deleteData()
        #     #    self.clearPointers(0)
        #     #    self.clearPointers(1)
        #     pass

        self.success = False
        # self.logger.debug("inside read_data")

        # Initiate the lists and arrays where the data will be read-in:
        ntotspw = self.spwlist[-1][3] + len(self.spwlist[-1][2])
        nsprang = range(ntotspw)
        self.u = [[] for sp in nsprang]
        self.v = [[] for sp in nsprang]
        self.w = [[] for sp in nsprang]
        self.t = [[] for sp in nsprang]
        self.tArr = [[] for sp in nsprang]
        self.tIdx = [[] for sp in nsprang]
        self.ant1 = [[] for sp in nsprang]
        self.ant2 = [[] for sp in nsprang]
        self.RAshift = [[] for sp in nsprang]
        self.Decshift = [[] for sp in nsprang]
        self.Stretch = [[] for sp in nsprang]
        self.averdata = [[] for sp in nsprang]
        self.avermod = [[] for sp in nsprang]
        self.averweights = [[] for sp in nsprang]
        self.averfreqs = [[] for sp in nsprang]
        self.iscancoords = [[] for sp in nsprang]

        # Read data for each spectral window:
        self.iscan = {}
        for vi in self.vis:
            self.iscan[vi] = {}

        for si in nsprang:
            self.logger.info(f"spectral index #{si} ({si+1} of {len(nsprang)})")
            # These are temporary lists of arrays that will be later concatenated:
            field_data = []
            max_tIndex = 0

            i0scan = 0

            # BEWARE! was sp
            for vis in [x for x in self.spwlist if x[3] <= si]:
                msname = self.vis[vis[1]]
                SPW, crosscorr = self.get_cross_correlations(msname)
                DDSC = self.get_data_description(msname)

                for spidx, spi in enumerate(vis[2]):
                    if vis[3] + spidx != si:
                        continue
                    sp = spi[0]
                    rang = spi[1]
                    self.iscan[msname][sp] = {}

                    DDs = np.where(DDSC == sp)[0]
                    if len(DDs) > 1:
                        self.logger.warning(f"spw {sp} has more than one Data Description ID!")

                    maskspw = np.zeros(np.shape(crosscorr), dtype=bool)
                    for ddi in DDs:
                        maskspw = np.logical_or(maskspw, SPW == ddi)
                    maskspw *= crosscorr

                    # For the first ms in the list, read the frequencies of the spw.
                    # All the other mss will be assumed to have the same frequencies:
                    # if True:
                    origfreqs = self.get_spectral_window(msname)
                    self.averfreqs[si] = np.array([np.average(origfreqs[r]) for r in rang])

                    self.logger.info(f"reading scans for spw {sp}")
                    # Read all scans for this field id:
                    for sc, scan in enumerate(self.sourscans[vis[1]]):
                        masksc, fieldids = self.get_scan_mask(msname, scan, maskspw)
                        for fieldid in fieldids:
                            self.logger.info(f"reading scan #{scan} "
                                             f"({sc+1} of {len(self.sourscans[vis[1]])}), field: {fieldid}")
                            fD = self.get_field_data(fieldid, scan, masksc, msname, vis[1],
                                                     rang, sp, si, max_tIndex, i0scan, takeModel)
                            print(fD)
                            field_data.append(fD)
                            max_tIndex = np.max(fD.tIndex) + 1
                            i0scan += np.shape(fD.datascanAv)[0]

            # Concatenate all the scans in one single array. Notice that we separate real and imag and save them
            # as floats. This is because ctypes doesn't handle complex128.
            self.averdata[si] = np.require(np.concatenate([fd.datascanAv for fd in field_data], axis=0),
                                           requirements=['C', 'A'])

            if takeModel:
                self.avermod[si] = np.require(np.concatenate([fd.modelscanAv for fd in field_data], axis=0),
                                              requirements=['C', 'A'])

            self.averweights[si] = np.require(np.concatenate([fd.weightscan for fd in field_data], axis=0),
                                              dtype=np.float64, requirements=['C', 'A'])
            self.u[si] = np.require(np.concatenate([fd.uscan for fd in field_data], axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.v[si] = np.require(np.concatenate([fd.vscan for fd in field_data], axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.w[si] = np.require(np.concatenate([fd.wscan for fd in field_data], axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.t[si] = np.require(np.concatenate([fd.tscan for fd in field_data], axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.tArr[si] = np.require(np.concatenate([fd.tArray for fd in field_data], axis=0),
                                       dtype=np.float64, requirements=['C', 'A'])
            self.tIdx[si] = np.require(np.concatenate([fd.tIndex for fd in field_data], axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])
            self.RAshift[si] = np.require(np.concatenate([fd.RAscan for fd in field_data], axis=0),
                                          dtype=np.float64, requirements=['C', 'A'])
            self.Decshift[si] = np.require(np.concatenate([fd.Decscan for fd in field_data], axis=0),
                                           dtype=np.float64, requirements=['C', 'A'])
            self.Stretch[si] = np.require(np.concatenate([fd.Stretchscan for fd in field_data], axis=0),
                                          dtype=np.float64, requirements=['C', 'A'])
            self.ant1[si] = np.require(np.concatenate([fd.ant1scan for fd in field_data], axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])
            self.ant2[si] = np.require(np.concatenate([fd.ant2scan for fd in field_data], axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])

        # Initial and final times of observations (time reference for proper motions):
        self.t0 = np.min([np.min(ti) for ti in self.t])
        self.t1 = np.max([np.max(ti) for ti in self.t])

        tac = time.time()
        self.logger.info(f"reading took {(tac - tic):.2f} seconds")

        self.success = True

        NIFs = len(self.averdata)
        self.Nspw = NIFs

        for i in range(self.Nspw):
            self.logger.info(f"there are {len(self.tArr[i])} integrations ({len(self.t[i])} visibs.) in spw {i}")

        return True

    def set_weight_equation(self, primary_beam_correction):
        self.logger.debug("MeasurementSet::set_weight_equation")

        if not isinstance(self.dish_diameter, dict):
            try:
                self.dish_diameter = float(self.dish_diameter)
            except ValueError:
                self.logger.error("The dish diameter must be a number! (in meters)")
                return False

        subdir = os.path.join(self.vis[0], 'ANTENNA')
        with open_tb(subdir) as tb:
            self.antnames = tb.getcol('NAME')
            self.Nants = len(self.antnames)
            self.logger.info(f"number of antennas = {self.Nants}")
            try:
                diameters = np.copy(tb.getcol('DISH_DIAMETER'))
            except Exception:
                self.logger.info("dish diameter column not found in antenna tables!")
                diameters = np.zeros(len(self.antnames))

        if primary_beam_correction:
            self.logger.info("You selected to apply primary-beam correction.\n"
                             "PLEASE, remember that the beam is being approximated\n"
                             "with a Gaussian, so it may not be very accuracte far\n"
                             "from the pointing direction.")
            if isinstance(self.dish_diameter, float):
                if self.dish_diameter != 0.0:
                    self.logger.info(f"an antenna diameter of {self.dish_diameter:.3f}m will be applied")
                    diameters = np.array([self.dish_diameter for a in self.antnames])

            elif isinstance(self.dish_diameter, dict):
                diameters = np.array([0.0 for a in self.antnames])
                for anam in self.dish_diameter.keys():
                    for i, name in enumerate(self.antnames):
                        antids = re.search(anam, name)
                        if 'start' in dir(antids):
                            diameters[i] = self.dish_diameter[anam]
                self.logger.info("manual antenna-size setting")
                for i, name in enumerate(self.antnames):
                    self.logger.info(f"antenna {name} has a diameter of {diameters[i]:.2f}m")

            else:
                self.logger.error("BAD dish_diameter! Should be a float or a dict!")

            if np.max(diameters) == 0.0:
                self.logger.error("The antenna diameters are not set in the ms. "
                                  "Please, set it manually or turn off primary-beam correction.")
                return False
            # Negative means not to apply PB corr for that antenna
            diameters[diameters == 0.0] = -1.0
            FWHM = self.ldfac / diameters * (2.99e8)
            sigma = FWHM / 2.35482 * (180. / np.pi) * 3600.
            self.KfacWgt = 1. / (2. * sigma**2.) * (diameters > 0.0)  # (0.5*(diameters/1.17741)**2.)
            self.userDiameters = diameters
        else:
            self.KfacWgt = np.zeros(len(self.antnames))

        # May refine this function in future releases:
        #  self.wgtEquation = lambda D, Kf: -D*Kf
        return True

    def write_model(self, column, model):
        """Writes the requested information into the measurement sets.

        The information can be either the predictions of the compiled model(s) (i.e., they are
        written into the *model column* of the measurement set(s)), or the post-fit residuals
        (they are written into the *corrected column*) or the calibrated data (they are written
        into the *corrected column* as well). The actual information to write is set by the value
        of the ``write`` keyword of the *UVMultiFit* instance when the ``fit()`` method was called.

        This function is executed only if the ``stokes`` keyword is set to either ``PI``, ``I``
        or an individual correlation product (like, e.g., ``XX`` or ``XY``) *and* if no averaging
        has been performed neither in time nor frequency (i.e., if both ``timewidth`` and
        ``chanwidth`` are set to 1). This function should be called AFTER having run ``fit()``.
        """

        self.logger.debug("MeasurementSet::write_model")
        self.logger.warning("writing to mosaics is experimental and may not work!")

        for v in self.vis:
            # Get the columns of parallel-hand correlations:
            spws = list(map(int, self.iscan[v].keys()))
            with open_ms(v) as ms:
                ms.selectinit(datadescid=spws[0])
                polprods = [x[0] for x in list(ms.range(['corr_names'])['corr_names'])]

            if self.stokes in polprods:
                polii = [polprods.index(self.stokes)]
            elif self.stokes in ['PI', 'I']:
                if 'XX' in polprods:
                    polii = [polprods.index('XX'), polprods.index('YY')]
                elif 'RR' in polprods:
                    polii = [polprods.index('RR'), polprods.index('LL')]
                else:
                    self.logger.error(f"Stokes not understood for '{v}', will not update the model column")
                    return False
            else:
                self.logger.error(f"Stokes not understood for '{v}', will not update the model column")
                return False

            # NEW CODE TO WRITE MODEL, BASED ON TB TOOL:
            for sp in spws:
                for scan in self.iscan[v][sp].keys():
                    for select in self.iscan[v][sp][scan]:
                        self.logger.info(f"doing '{v}': spw {sp}, scan_id {scan}")
                        with open_tb(v) as tb:
                            tb2 = tb.selectrows(select[-1])
                            moddata = tb2.getcol(column)
                            re = np.transpose(model.output[select[0]][select[1]:select[1] + select[2], :])
                            for nui, r in enumerate(select[3]):
                                for poli in polii:
                                    moddata[poli, r, :] = re[nui, :]
                            tb2.putcol(column, moddata)

            self.logger.info(f"{column} written successfully")

        return True

    @staticmethod
    def channeler(spw, width=1, maxchans=[3840, 3840, 3840, 3840]):
        """ Function to convert a string with spw selection into lists of channels to select/average."""
        logging.debug(f"channeler({spw})")

        if spw == '':
            spw = ','.join(list(map(str, range(len(maxchans)))))

        entries = spw.split(',')
        output = [[] for i in maxchans]

        for entry in entries:
            check = entry.split(':')
            if check[0] == '*':
                check[0] = f"0~{len(maxchans)-1}"

            spws = list(map(int, check[0].split('~')))
            if len(spws) == 1:
                selspw = [spws[0]]
            else:
                selspw = list(range(spws[0], spws[1] + 1))

            for sp in selspw:
                if sp + 1 > len(maxchans):
                    errstr = f"There are only {len(maxchans)} spw in the data, please revise the 'spw' parameter."
                    logging.error(errstr)
                    return [False, errstr]

                if len(check) == 1:
                    channel_ranges = [f"0~{maxchans[sp] - 1}"]
                else:
                    chans = check[1]
                    channel_ranges = chans.split(';')
                logging.debug(f"channel_range {sp}: {channel_ranges}")
                ranges = []
                for channel_range in channel_ranges:
                    ch1, ch2 = list(map(int, channel_range.split('~')))
                    if ch1 > ch2:
                        errstr = f"Channel {ch1} is larger than {ch2}, revise channels for spw {sp}."
                        logging.error(errstr)
                        return [False, errstr]
                    ch2 = min([ch2, maxchans[sp] - 1])
                    for i in range(ch1, ch2 + 1, width):
                        ranges.append(list(range(i, min([(i + width), ch2 + 1]))))

                output[sp] = ranges
        return [True, output]

if __name__ == "__main__":
    # data = MeasurementSet('foo', spw=['0'])
    # data.check_measurementset()
    # print(data)

    data = MeasurementSet('../test-cases/Disc/Disc.alma.out10.noisy.ms')  # , stokes='YY')
    data.check_measurementset()
    data.read_data()
    # print(data.channeler('')[0])
    # print(data.channeler('0')[0])
    # print(data.channeler('1')[0])
    # print(data.channeler('0~3')[0])
    # print(data.channeler('0,2,3')[0])
    # print(data.channeler('1,3,5')[0])
    # print(data.channeler('*:10~50')[0])
    # print(data.channeler('1~3:30~40')[0])
    # data.dump()
    # print(data)

    # data = MeasurementSet('../test-cases/Disc/Disc.alma.out10.noisy.ms', field='Disc.alma.out10_0',
    #                       phase_center='J2000 12h34m56.0s 01d02m03.0s')
    # data.check_measurementset()
    # print(data)
