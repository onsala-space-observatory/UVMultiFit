import os
import re
import time
import logging
import numpy as np

from casatools import ms
from casatools import table
from casatools import coordsys

from .utils import get_list_of_strings, is_list_of_int, is_list_of_floats
from .utils import is_casa_position, is_valid_stokes

class MeasurementSet():
    """ Class to deal with model equations and fitting procedures.

    If convert strings representing models and parameters into compiled equations, to be used in a ChiSq
    visibility fitting. It also interacts with the C++ extension of UVMultiFit.

    This class should NOT be instantiated by the user (it is called from the UVMultiFit class."""

    logger = logging.getLogger("measurement")
    tb = table()
    ms = ms()
    cs = coordsys()

    def __init__(self, vis, spw='0', field=0, scans=[], corrected=False,
                 uniform=False, uvtaper=0.0, chanwidth=1, timewidth=1, stokes='I',
                 MJDrange=[-1.0, -1.0], ldfac=1.22, phase_center='', pbeam=False,
                 wgt_power=1.0, dish_diameter=0.0):
        self.vis = get_list_of_strings(vis)
        self.spw = get_list_of_strings(spw)
        if len(self.spw) > 1 and len(self.vis) != len(self.spw):
            self.logger.error("the length of 'spw' is not equal to the length of 'vis'!")

        if not isinstance(field, (int, str)):
            self.logger.warning("parameter 'field' needs to be an integer or string, resetting to 0!")
            self.field = 0
        else:
            self.field = field
        if isinstance(self.field, str) and self.field.isdigit():
            self.field = int(self.field)
        self.field_id = []
        self.pointing = []
        self.sourscans = []
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
        self.spwlist = []
        self.pol2aver = []
        self.polmod = []
        self.polii = []

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

        csys = self.cs.newcoordsys(direction=True)
        dirstr = phase_center.split()
        if len(dirstr) == 2:
            csys.setdirection(refcode="J2000", refval=phase_center)
        else:
            csys.setdirection(refcode=dirstr[0], refval=" ".join(dirstr[1:]))
            csys.convertdirection("J2000")
        refpos = np.copy(csys.torecord()['direction0']['crval'])
        return refpos

    def check_measurementset(self):
        self.logger.debug("MeasurementSet::check_measurementset")
        for v in self.vis:
            if not os.path.exists(v):
                self.logger.error(f"measurement set '{v}' does not exist!")
                return False
        phasedirs = [{} for v in self.vis]

        # Get number of antennas:
        self.tb.open(os.path.join(self.vis[0], 'ANTENNA'))
        self.Nants = len(self.tb.getcol('NAME'))
        self.logger.info(f"number of antennas = {self.Nants}")
        self.tb.close()

        # Open MS and look for the selected data:
        for vi, v in enumerate(self.vis):
            self.field_id.append([])
            success = self.ms.open(v)
            if not success:
                self.logger.error(f"failed to open measurement set '{v}'!")
                return False

            allfields = list(self.ms.range('fields')['fields'])
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
                    self.logger.error(f"field '{self.field[vi]}' is not in '{v}'")
                    self.ms.close()
                    return False

            self.field_id[-1].append(field_index)
            phasedirs[vi][field_index] = self.ms.range('phase_dir')['phase_dir']['direction'][:, field_index]
            self.ms.close()
        if self.refpos is None:
            self.refpos = phasedirs[0][min(phasedirs[0].keys())]
        success = self.find_observing_positions(phasedirs)
        success = success and self.get_spectral_configuration()
        success = success and self.get_polarization_configuration()
        success = success and self.set_weight_equation(self.pbeam)
        self.phasedirs = phasedirs
        return success

    def find_observing_positions(self, phasedirs):
        self.logger.debug("MeasurementSet::find_observing_positions")
        for vi, v in enumerate(self.vis):
            success = self.ms.open(v)
            if not success:
                self.logger.error(f"failed to open measurement set '{v}'!")
                return False

            info = self.ms.getscansummary()
            self.ms.close()

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

            if not self.ms.open(v):
                self.logger.error(f"failed to open measurement set '{v}'!")
                return False

            # freqdic = ms.getspectralwindowinfo()
            spwchans = self.ms.range(['num_chan'])['num_chan']

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

            self.ms.close()
        return True

    def get_polarization_configuration(self):
        self.logger.debug("MeasurementSet::get_polarization_configuration")
        for v in self.vis:
            if not self.ms.open(v):
                self.logger.error(f"failed to open measurement set '{v}'!")
                return False

            polprods = [x[0] for x in list(self.ms.range(['corr_names'])['corr_names'])]

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

            self.ms.close()
        return True

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

        # maxDist = 0.0

        # Read data for each spectral window:
        self.iscan = {}
        for vi in self.vis:
            self.iscan[vi] = {}

        for si in nsprang:
            self.logger.info(f"spectral index #{si} ({si+1} of {len(nsprang)})")
            # These are temporary lists of arrays that will be later concatenated:
            datascanAv = []
            modelscanAv = []
            #   datascanim = []
            weightscan = []
            uscan = []
            vscan = []
            wscan = []
            tscan = []
            tArray = []
            tIndex = []
            ant1scan = []
            ant2scan = []
            RAscan = []
            Decscan = []
            Stretchscan = []

            i0scan = 0

            # BEWARE! was sp
            for vis in [x for x in self.spwlist if x[3] <= si]:
                msname = self.vis[vis[1]]

                self.logger.info(f"opening measurement set '{msname}'")
                self.tb.open(msname)
                SPW = self.tb.getcol('DATA_DESC_ID')
                crosscorr = self.tb.getcol('ANTENNA1') != self.tb.getcol('ANTENNA2')
                self.tb.close()

                self.tb.open(os.path.join(msname, 'DATA_DESCRIPTION'))
                DDSC = self.tb.getcol('SPECTRAL_WINDOW_ID')
                self.tb.close()

                for spidx, spi in enumerate(vis[2]):
                    if vis[3] + spidx == si:

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
                        self.tb.open(os.path.join(msname, 'SPECTRAL_WINDOW'))
                        origfreqs = self.tb.getcol('CHAN_FREQ')
                        self.averfreqs[si] = np.array([np.average(origfreqs[r]) for r in rang])
                        nfreq = len(rang)
                        self.tb.close()
                        self.logger.info(f"reading scans for spw {sp}")

                        # Read all scans for this field id:
                        for sc, scan in enumerate(self.sourscans[vis[1]]):
                            self.tb.open(msname)

                            masksc = maskspw * (self.tb.getcol('SCAN_NUMBER') == int(scan))
                            fieldids = list(np.sort(np.unique(self.tb.getcol('FIELD_ID')[masksc])))

                            self.tb.close()

                            for fieldid in fieldids:
                                self.logger.info(f"reading scan #{scan} "
                                                 f"({sc+1} of {len(self.sourscans[vis[1]])}), field: {fieldid}")
                                self.tb.open(msname)
                                maskfld = np.where(masksc * (self.tb.getcol('FIELD_ID') == int(fieldid)))[0]

                                if len(maskfld) == 0:
                                    self.tb.close()
                                else:
                                    tb2 = self.tb.selectrows(maskfld)
                                    uvscan = {'uvw': tb2.getcol('UVW'), 'antenna1': tb2.getcol('ANTENNA1'),
                                              'antenna2': tb2.getcol('ANTENNA2'), 'time': tb2.getcol('TIME')}
                                    times = np.unique(uvscan['time'])
                                    if takeModel:
                                        # datascan = ms.getdata([self.column, 'model_data', 'weight', 'flag'],
                                        #                       ifraxis=True)
                                        datascan = {self.column: tb2.getcol((self.column).upper()),
                                                    'model_data': tb2.getcol('MODEL_DATA'),
                                                    'weight': tb2.getcol('WEIGHT'),
                                                    'flag': tb2.getcol('FLAG')}
                                    else:
                                        # datascan = ms.getdata([self.column, 'weight', 'flag'], ifraxis=True)
                                        datascan = {self.column: tb2.getcol((self.column).upper()),
                                                    'weight': tb2.getcol('WEIGHT'), 'flag': tb2.getcol('FLAG')}

                                    self.tb.close()

                                    # NOTE: There is a bug in np.ma.array that casts complex
                                    # to float under certain operations (e.g., np.ma.average).
                                    # That's why we average real and imag separately.

                                    # Compute the polarization product:

                                    # Bad data has zero weight:
                                    datascan['weight'][np.logical_not(np.isfinite(datascan['weight']))] = 0.0

                                    # All unflagged weights set to equal (if uniform):
                                    if self.uniform:
                                        datascan['weight'][datascan['weight'] > 0.0] = 1.0

                                    datascan['weight'][datascan['weight'] < 0.0] = 0.0
                                    copyweight = np.copy(datascan['weight'])

                                    totalmask = datascan['flag']
                                    origmasked = np.ma.array(datascan[self.column], mask=totalmask, dtype=np.complex128)

                                    if takeModel:
                                        origmodmasked = np.ma.array(
                                            datascan['model_data'],
                                            mask=totalmask, dtype=np.complex128)

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
                                    polavg = [pol != 0.0 for pol in self.pol2aver[vis[1]]]

                                    if self.polmod[vis[1]] == 2:
                                        flagmask = np.ma.logical_and(
                                            totalmask[self.polii[vis[1]][0], :, :],
                                            totalmask[self.polii[vis[1]][1], :, :])
                                        datamask = np.ma.average(origmasked[polavg, :].real, axis=0) + \
                                            1.j * np.ma.average(origmasked[polavg, :].imag, axis=0)

                                        if takeModel:
                                            modelmask = np.ma.average(origmodmasked[polavg, :].real, axis=0) + \
                                                1.j * np.ma.average(origmasked[polavg, :].imag, axis=0)

                                        weightmask = np.sum(origweight[polavg, :], axis=0)
                                    else:
                                        if self.polmod[vis[1]] == 3:
                                            flagmask = totalmask[self.polii[vis[1]][0], :, :]
                                        else:
                                            flagmask = np.ma.logical_or(totalmask[self.polii[vis[1]][0], :, :],
                                                                        totalmask[self.polii[vis[1]][1], :, :])

                                        for pol in self.polii[vis[1]]:
                                            datamask += origmasked[pol, :] * self.pol2aver[vis[1]][pol]
                                            if takeModel:
                                                modelmask += origmodmasked[pol, :] * self.pol2aver[vis[1]][pol]
                                            weightmask += origweight[pol, :]

                                        if self.polmod[vis[1]] == 1:
                                            datamask *= 1.j
                                            if takeModel:
                                                modelmask *= 1.j
                                                #   weightmask[flagmask] = 0.0

                                    # Free some memory:
                                    for key in list(datascan):
                                        del datascan[key]
                                    del datascan, origmasked, origweight
                                    del copyweight, totalmask, flagmask
                                    if takeModel:
                                        del origmodmasked

                                    # Mosaic-related corrections:
                                    phshift = 3600. * 180. / np.pi * (self.phasedirs[vis[1]][fieldid] - self.refpos)
                                    strcos = np.cos(self.phasedirs[vis[1]][fieldid][1])

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

                                    # Average in time and apply uvtaper:
                                    # GaussWidth = 2. * (self.uvtaper / 1.17741)**2.
                                    RAoffi = np.zeros(np.shape(uvscan['time']))
                                    Decoffi = np.copy(RAoffi)
                                    Stretchi = np.copy(RAoffi)

                                    #  if self.timewidth ==1:

                                    RAoffi[:] = float(phshift[0])
                                    Decoffi[:] = float(phshift[1])
                                    Stretchi[:] = float(strcos)

                #########################
                # CODE FOR TIMEWIDTH>1 HAS TO BE BASED ON TB TOOL. WORK IN PROGRESS
                #     else:

                #       ant1s = uvscan['antenna1'][crosscorr]
                #       ant2s = uvscan['antenna2'][crosscorr]
                #       maskan1 = np.zeros(np.shape(ant1s), dtype=np.bool)
                #       mask = np.copy(maskan1)
                #       mask2 = np.copy(mask)
                #       # Antennas participating in this scan:
                #       allants1 = np.unique(ant1s)
                #       allants2 = np.unique(ant2s)
                #       # Fill in time-averaged visibs:
                #       for nt in range(ntav):
                #         t0 = nt*self.timewidth ; t1 = min([ntimes,(nt+1)*self.timewidth])
                #         mask[:] = (uvscan['time']>=times[t0])*(uvscan['time']<times[t1])*crosscorr
                #         for an1 in allants1:
                #           maskan1[:] = (ant1s==an1)*mask
                #           for an2 in allants2:
                #             if an2>an1:
                #                 mask2[:] = maskan1*(ant2s==an2)
                #                 uu = uvscan['uvw'][0, mask2]
                #                 vv = uvscan['uvw'][1, mask2]
                #                 ww = uvscan['uvw'][2, mask2]
                #                 ui[:, nt] = np.average(uu, axis=1)  # Baseline dimension
                #                 vi[:, nt] = np.average(vv, axis=1)  # Baseline dimension
                #                 wi[:, nt] = np.average(ww, axis=1)  # Baseline dimension
                #                 ant1i[:, nt] = ant1s  # Baseline dimension
                #                 ant2i[:, nt] = ant2s  # Baseline dimension
                #                 timei[:, nt] = np.average(uvscan['time'][t0:t1])/86400.
                #                 tArrayi[nt] = time[0, nt]  # np.average(uvscan['time'][t0:t1])/86400.
                #                 tIndexi[:, nt] = nt
                #                 RAoffi[:, nt] = float(phshift[0])
                #                 Decoffi[:, nt] = float(phshift[1])
                #                 Stretchi[:, nt] = float(strcos)
                #
                #         if self.uvtaper > 0.0:
                #           GaussFact = np.exp(-(uu*uu + vv*vv)/GaussWidth)
                #         else:
                #           GaussFact = np.ones(np.shape(uu))
                #
                #       broadwgt = weighttemp[:, :, t0:t1]
                #       avercompl[:, :, nt] = np.ma.average(datatemp[:, :, t0:t1].real, axis=2, weights=broadwgt)+
                #                             1.j*np.ma.average(datatemp[:, :, t0:t1].imag, axis=2, weights=broadwgt)
                #       if takeModel:
                #         avermodl[:, :, nt] = np.ma.average(modeltemp[:, :, t0:t1].real, axis=2, weights=broadwgt)+
                #                              1.j*np.ma.average(modeltemp[:, :, t0:t1].imag, axis=2, weights=broadwgt)
                #
                #       if self.uniform:
                #         averwgt[:, :, nt] = np.ma.sum(np.ones(np.shape(broadwgt))*GaussFact[np.newaxis, :, :], axis=2)
                #       else:
                #         averwgt[:, :, nt] = np.ma.sum(broadwgt*GaussFact[np.newaxis, :, :], axis=2)
                #########################

                                    ant1scan.append(np.copy(uvscan['antenna1'][:]))
                                    ant2scan.append(np.copy(uvscan['antenna2'][:]))
                                    uscan.append(np.copy(uvscan['uvw'][0, :]))
                                    vscan.append(np.copy(uvscan['uvw'][1, :]))
                                    wscan.append(np.copy(uvscan['uvw'][2, :]))
                                    tscan.append(np.copy(uvscan['time'][:]))
                                    tArray.append(np.copy(times))

                                    tIndexi = np.zeros(np.shape(uvscan['time']), dtype=np.int32)
                                    for tid, tiii in enumerate(times):
                                        tIndexi[uvscan['time'] == tiii] = tid

                                    if len(tIndex) > 1:
                                        tIndexi += np.max(tIndex[-1]) + 1
                                    tIndex.append(tIndexi)

                                    RAscan.append(RAoffi)
                                    Decscan.append(Decoffi)
                                    Stretchscan.append(Stretchi)

                                    datascanAv.append(np.transpose(datatemp))
                                    if takeModel:
                                        modelscanAv.append(np.transpose(modeltemp))
                                    if self.uniform:
                                        weightscan.append(np.transpose(np.ones(np.shape(weighttemp))))
                                    else:
                                        weightscan.append(np.transpose(weighttemp))

                                    # Useful info for function writeModel() and for pointing correction:
                                    if scan not in self.iscan[msname][sp].keys():
                                        self.iscan[msname][sp][scan] = []

                                    self.iscan[msname][sp][scan].append(
                                        [int(si), int(i0scan), int(ndata), list(rang), np.copy(maskfld)])
                                    self.iscancoords[si].append([i0scan, i0scan + ndata, phshift[0], phshift[1]])

                                    i0scan += ndata

            # Concatenate all the scans in one single array. Notice that we separate real and imag and save them
            # as floats. This is because ctypes doesn' t handle complex128.
            self.averdata[si] = np.require(np.concatenate(datascanAv, axis=0),
                                           requirements=['C', 'A'])  # , np.concatenate(datascanim, axis=0)]

            if takeModel:
                self.avermod[si] = np.require(np.concatenate(modelscanAv, axis=0),
                                              requirements=['C', 'A'])

            self.averweights[si] = np.require(np.concatenate(weightscan, axis=0),
                                              dtype=np.float64, requirements=['C', 'A'])
            self.u[si] = np.require(np.concatenate(uscan, axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.v[si] = np.require(np.concatenate(vscan, axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.w[si] = np.require(np.concatenate(wscan, axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.t[si] = np.require(np.concatenate(tscan, axis=0),
                                    dtype=np.float64, requirements=['C', 'A'])
            self.tArr[si] = np.require(np.concatenate(tArray, axis=0),
                                       dtype=np.float64, requirements=['C', 'A'])
            self.tIdx[si] = np.require(np.concatenate(tIndex, axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])
            self.RAshift[si] = np.require(np.concatenate(RAscan, axis=0),
                                          dtype=np.float64, requirements=['C', 'A'])
            self.Decshift[si] = np.require(np.concatenate(Decscan, axis=0),
                                           dtype=np.float64, requirements=['C', 'A'])
            self.Stretch[si] = np.require(np.concatenate(Stretchscan, axis=0),
                                          dtype=np.float64, requirements=['C', 'A'])
            self.ant1[si] = np.require(np.concatenate(ant1scan, axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])
            self.ant2[si] = np.require(np.concatenate(ant2scan, axis=0),
                                       dtype=np.int32, requirements=['C', 'A'])

            # Free some memory:
            #   del uu, vv, ww, avercomplflat, weightflat
            del datatemp, weighttemp, uvscan  # , avercompl, averwgt
            for dda in datascanAv:
                del dda
            #   for dda in datascanim:
            #     del dda
            if takeModel:
                for dda in modelscanAv:
                    del dda
            for dda in weightscan:
                del dda
            for dda in tscan:
                del dda
            for dda in uscan:
                del dda
            for dda in vscan:
                del dda
            for dda in wscan:
                del dda
            for dda in RAscan:
                del dda
            for dda in Decscan:
                del dda
            for dda in Stretchscan:
                del dda
            for dda in ant1scan:
                del dda
            for dda in ant2scan:
                del dda
            for dda in tArray:
                del dda
            for dda in tIndex:
                del dda

            del datascanAv  # , datascanim
            del weightscan, tscan, uscan, vscan, wscan, tArray, tIndex
            del RAscan, Decscan, Stretchscan, ant1scan, ant2scan
            if takeModel:
                del modelscanAv

            # try:
            #     del GaussFact
            # except Exception:
            #     pass
            # gc.collect()

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
        tempfloat = 0.0

        if not isinstance(self.dish_diameter, dict):
            try:
                self.dish_diameter = float(self.dish_diameter)
            except ValueError:
                self.logger.error("The dish diameter must be a number! (in meters)")
                return False

        self.tb.open(os.path.join(self.vis[0], 'ANTENNA'))
        self.antnames = self.tb.getcol('NAME')

        if primary_beam_correction:
            self.logger.info("You selected to apply primary-beam correction.\n"
                             "PLEASE, remember that the beam is being approximated\n"
                             "with a Gaussian, so it may not be very accuracte far\n"
                             "from the pointing direction.")
            if isinstance(self.dish_diameter, float):
                if self.dish_diameter == 0.0:
                    try:
                        tempfloat = np.copy(self.tb.getcol('DISH_DIAMETER'))
                    except Exception:
                        self.logger.info("dish diameter column not found in antenna tables!")
                    tempfloat = np.zeros(len(self.antnames))
                else:
                    self.logger.info(f"an antenna diameter of {self.dish_diameter:.3f}m will be applied")
                    tempfloat = np.array([self.dish_diameter for a in self.antnames])

            elif isinstance(self.dish_diameter, dict):
                tempfloat = np.array([0.0 for a in self.antnames])
                for anam in self.dish_diameter.keys():
                    for i, name in enumerate(self.antnames):
                        antids = re.search(anam, name)
                        if 'start' in dir(antids):
                            tempfloat[i] = self.dish_diameter[anam]
                self.logger.info("manual antenna-size setting")
                for i, name in enumerate(self.antnames):
                    self.logger.info(f"antenna {name} has a diameter of {tempfloat[i]:.2f}m")

            else:
                self.logger.error("BAD dish_diameter! Should be a float or a dict!")

            if np.max(tempfloat) == 0.0:
                self.logger.error("The antenna diameters are not set in the ms. "
                                  "Please, set it manually or turn off primary-beam correction.")
                return False
            # Negative means not to apply PB corr for that antenna
            tempfloat[tempfloat == 0.0] = -1.0
            FWHM = self.ldfac / tempfloat * (2.99e8)
            sigma = FWHM / 2.35482 * (180. / np.pi) * 3600.
            self.KfacWgt = 1. / (2. * sigma**2.) * (tempfloat > 0.0)  # (0.5*(tempfloat/1.17741)**2.)
            self.userDiameters = tempfloat
        else:
            self.KfacWgt = np.zeros(len(self.antnames))

        self.tb.close()
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
            success = self.ms.open(v)
            if not success:
                self.logger.error(f"'{v}' cannot be openned in write mode")
                return False

            spws = list(map(int, self.iscan[v].keys()))

            self.ms.selectinit(datadescid=spws[0])
            polprods = [x[0] for x in list(self.ms.range(['corr_names'])['corr_names'])]
            self.ms.close()

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

            if self.write_model == 1:
                column = 'MODEL_DATA'
            elif self.write_model in [2, 3]:
                column = 'CORRECTED_DATA'

            # NEW CODE TO WRITE MODEL, BASED ON TB TOOL:
            for sp in spws:
                for scan in self.iscan[v][sp].keys():
                    for select in self.iscan[v][sp][scan]:
                        self.logger.info(f"doing '{v}': spw {sp}, scan_id {scan}")
                        self.tb.open(v, nomodify=False)
                        tb2 = self.tb.selectrows(select[-1])
                        moddata = tb2.getcol(column)
                        re = np.transpose(model.output[select[0]][select[1]:select[1] + select[2], :])
                        for nui, r in enumerate(select[3]):
                            for poli in polii:
                                moddata[poli, r, :] = re[nui, :]
                        tb2.putcol(column, moddata)
                        self.tb.close()

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
                    errstr = f"there are only {len(maxchans)} spw in the data, please revise the 'spw' parameter"
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
                        errstr = f"{ch1} is larger than {ch2}, revise channels for spw {sp}"
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
