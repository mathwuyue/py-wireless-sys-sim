import numpy as np
from core.positions import cal_dist_3d


class Satellite(object):
    """ The clasatellite for creating a satellite system. """
    def __init__(self, satellites=None, antennas=None, stations=None):
        """
        Kwargs:
        satellites (numpy array): positions of satellites
        speeds (numpy array): np.array([[1,2,3], [2,3,4]])
        antennas (list): list of antennas
        stations (list): list of earth stations
        .. note::
        satellites should have same shape as speeds
        """
        self.satellites = satellites
        self.antennas = antennas
        self.stations = stations

    def add_satellite(self, satellites, antennas):
        if self.satellites is None:
            self.satellites = satellites
            self.antennas = antennas
        else:
            self.satellites = np.append(self.satellites, satellites)
            self.antennas.append(antennas)

    def add_stations(self, stations):
        self.stations.append(stations)

    def update_pos(self, t):
        pass

    def cal_satellites_topo(self):
        pass


class IridiumSatellite(Satellite):
    def __init__(self):
        super(IridiumSatellite, self).__init__()

    def cal_satellites_topo(self):
