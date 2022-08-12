import logging

from utils import get_list_of_strings, is_list_of_int

class MeasurementSet():
    """ Class to deal with model equations and fitting procedures.

    If convert strings representing models and parameters into compiled equations, to be used in a ChiSq
    visibility fitting. It also interacts with the C++ extension of UVMultiFit.

    This class should NOT be instantiated by the user (it is called from the UVMultiFit class."""

    logger = logging.getLogger("measurment")

    def __init__(self, vis, spw='0', field=0, scans=[], corrected=False):
        logging.basicConfig(level=logging.INFO,
                            format='%(name)s - %(levelname)s - %(message)s')
        self.vis = get_list_of_strings(vis)
        self.spw = get_list_of_strings(spw)

        self.field = field
        if is_list_of_int(scans):
            self.scans = scans
        else:
            self.logger.warning("'scans' should be list of (lists of) integers, resetting to empty list")
            self.scans = []
        self.column = 'corrected_data' if corrected else 'data'

    def __repr__(self):
        txt = "MeasurementSet("
        txt += f"vis = '{self.vis}', "
        txt += f"spw = {self.spw}, "
        txt += f"column = '{self.column}', "
        txt += f"field = {self.field}, "
        txt += f"scans = {self.scans})"
        return txt
