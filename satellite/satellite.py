import numpy as np
from core.positions import cal_dist_3d


class KeplerElement(object):
    """
    See: https://en.wikipedia.org/wiki/Orbital_elements
    """
    def __init__(self, e, a, i, omega, av, epoch):
        """
        Args:
        e (float): eccentricity. Circle is 0.
        a (float): semimajor axis. Radius of circle
        i (float): inclination.
        omega (float): longitude of the ascending node
        av (float): Circle: angular velocity
        epoch (float): Mean anomaly at epoch.
        """
        self.e = e
        self.a = a
        self.i = i
        self.omega = omega
        self.av = av
        self.epoch = epoch


class Satellite(object):
    """
    The class for creating a satellite system.  It is assumed that all satellites have a circular orbit.
    """
    def __init__(self, kes=None, antennas=None, intra_bw=0, earth_bw=0, init_pos=None, stations=None):
        """
        Kwargs:
        kes (list): array of Kepler elements of kes
        antennas (list): list of antennas
        stations (list): list of earth stations
        init_pos (list): initial positions of satellites (in sphercial coordination system)
        .. note::
        kes should have same shape as init_pos
        """
        self.kes = kes
        self.antennas = antennas
        self.stations = stations
        self.intra_bw = intra_bw
        self.earth_bw = earth_bw
        self.init_pos = init_pos

    def add_satellite(self, kes, antennas):
        if self.kes is None:
            self.kes = kes
            self.antennas = antennas
        else:
            self.kes = np.append(self.kes, kes)
            self.antennas.append(antennas)

    def add_stations(self, stations):
        self.stations.append(stations)

    def update_pos(self, t):
        pass

    def cal_satellite_topo(self):
        pass


class IridiumSatellite(Satellite):
    def __init__(self, intra_bw=1, earth_bw=1, stations=None):
        """
        Kwargs:
        intra_bw (float): in MHz
        earth_bw (float): in MHz
        """
        super(IridiumSatellite, self).__init__(intra_bw=intra_bw, earth_bw=earth_bw, stations=stations)

    def cal_satellite_topo(self):
