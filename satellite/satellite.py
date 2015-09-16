import numpy as np
from core.positions import cal_dist_3d
from core.antenna import SatelliteAntenna, ParabolicAntenna


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
        init_pos (numpy array): initial positions of satellites (in sphercial coordination system)
        .. note::
        kes, init_pos and antennas should have same number of rows.
        """
        self.kes = kes
        self.antennas = antennas
        self.stations = stations
        self.intra_bw = intra_bw
        self.earth_bw = earth_bw
        self.pos = init_pos

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
        if self.pos is None:
            return
        if type(kes) is not list:
            self.pos = np.array([self.pos[0], self.pos[1]+self.kes.av*t, self.pos[2]])
        self.pos = np.array([[p[0], p[1]+ke.av*t, p[2]]
                             for p, ke in zip(self.pos, self.kes)])

    def to_cartesian(self, idx=None):
        if self.pos is None:
            return None
        if type(idx) is int:
            p = self.pos
            return np.array([self.p[0]*np.sin(self.p[1])*np.sin(self.p[2]),
                             self.p[0]*np.sin(self.p[1])*np.cos(self.p[2]),
                             self.p[0]*np.cos(self.p[1])])
        if idx is None:
            rpos = self.pos
        else:
            rpos = self.pos[idx, :]
        return np.array([[p[0]*np.sin(p[1])*np.sin(p[2]),
                          p[0]*np.sin(p[1])*np.cos(p[2]),
                          p[0]*np.cos(p[1])] for p in rpos])

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
        self.kes = []
        self.pos = []
        self.antennas = []
        for i in range(6):
            for j in range(11):
                self.kes.append(KeplerElement(0, 7136e3, np.pi/2, i*np.pi/6, 0, 0))
                self.antennas.append(SatelliteAntenna(27.5, 6, 30, 0, 4.2))
                if i % 2 == 0:
                    self.pos.append([7136e3, 2*np.pi/11*j, i*np.pi/6])
                else:
                    self.pos.append([7136e3, np.pi/11+2*np.pi/11*j, i*np.pi/6])
        self.pos = np.array(self.pos)

    def cal_satellite_topo(self):
        topo = []
        cpos = self.to_cartesian()
        for i in range(66):
            self.row = []
            for j in range(66):
                if i == j:
                    self.row.append(np.inf)
                if j == i-1 or j == i+1 or j == i+11 or j == i-11:
                    self.row.append(cal_dist_3d(cpos[i, :], cpos[j, :]))
            topo.append(self.row)
        return np.array(topo)


class GeoSatellite(Satellite):
    def __init__(self, intra_bw=1, earth_bw=1, n=3, stations=None):
        super(GeoSatellite, self).__init__(intra_bw=intra_bw, earth_bw=earth_bw, stations=stations)
        self.kes = [KeplerElement(0, 42164e3, 0, 0, 0, 0) for i in range(n)]
        self.pos = np.array([[42164e3, np.pi/2, i*2*np.pi/n] for i in range(n)])
        self.antennas = [SatelliteAntenna(27.5, 6, 30, 0, 4.2) for i in range(n)]
