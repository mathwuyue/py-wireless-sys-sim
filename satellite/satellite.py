import numpy as np
from core.position import cal_dist_3d, to_cartesian
from core.antenna import SatelliteAntenna, ParabolicAntenna
import networkx as nx


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
    def __init__(self, n=0, kes=None, antennas=None, intra_bw=0, earth_bw=0, init_pos=None, stations=None):
        """
        Kwargs:
        kes (list): array of Kepler elements of kes
        antennas (list): list of antennas
        stations (list): list of earth stations
        init_pos (numpy array): initial positions of satellites (in sphercial coordination system)
        .. note::
        kes, init_pos and antennas should have same number of rows.
        """
        self.n = n
        self.kes = kes
        self.antennas = antennas
        self.stations = stations
        self.intra_bw = intra_bw
        self.earth_bw = earth_bw
        self.pos = init_pos
        self.topo = None

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
        if type(self.kes) is not list:
            self.pos = np.array([self.pos[0], self.pos[1]+self.kes.av*t, self.pos[2]])
        else:
            self.pos = np.array([[p[0], p[1]+ke.av*t, p[2]]
                                 for p, ke in zip(self.pos, self.kes)])

    def satellite_pos(self, idx=None):
        if self.pos is None:
            return None
        if idx is None:
            rpos = self.pos
        else:
            rpos = self.pos[idx, :]
        return to_cartesian(rpos)

    def get_antenna_param(self, idx=None, param=None):
        if self.pos is None:
            return None
        if type(idx) is np.int64 or type(idx) is int:
            return self.antennas[idx][param]
        if idx is None:
            return np.array([a[param] for a in self.antennas])
        else:
            return np.array([self.antennas[i][param] for i in idx])

    def get_satellite_topo(self):
        pass


class IridiumSatellite(Satellite):
    def __init__(self, intra_bw=1, earth_bw=1, stations=None):
        """
        Kwargs:
        intra_bw (float): in MHz
        earth_bw (float): in MHz
        """
        super(IridiumSatellite, self).__init__(n=66, intra_bw=intra_bw, earth_bw=earth_bw, stations=stations)
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

    def get_satellite_topo(self):
        if self.topo is not None:
            return self.topo
        self.topo = nx.Graph()
        self.topo.add_nodes_from(range(66))
        for i in range(66):
            if i % 11 != 0:
                self.topo.add_edge(i-1, i)
            if i+1 != 10:
                self.topo.add_edge(i+1, i)
            if i+11 < 65:
                self.topo.add_edge(i+11, i)
            if i-11 >= 0:
                self.topo.add_edge(i-11, i)


class GeoSatellite(Satellite):
    def __init__(self, n=3, intra_bw=1, earth_bw=1, stations=None):
        super(GeoSatellite, self).__init__(n=n, intra_bw=intra_bw, earth_bw=earth_bw, stations=stations)
        self.kes = [KeplerElement(0, 42164e3, 0, 0, 0, 0) for i in range(n)]
        self.pos = np.array([[42164e3, np.pi/2, i*2*np.pi/n] for i in range(n)])
        self.antennas = [SatelliteAntenna(27.5, 6, 30, 0, 4.2) for i in range(n)]


class EarthStation(object):
    def __init__(self, bw=1, pos=None, antennas=None):
        """
        Kwargs:
        bw (float): MHz
        pos (numpy array): (r, i, theta)
        antennas (list)
        .. note::
        pos have same rows as antennas, except that pos is not None but antennas is None, in which case antennas will be set to the array of default parabolicantennas. 
        """
        self.bw = bw
        self.pos = pos
        if pos is not None and antennas is None:
            self.antennas = [ParabolicAntenna() for p in self.pos]
        else:
            self.antennas = antennas

    def update_pos(self, t):
        pass

    def station_pos(self, idx=None):
        if self.pos is None:
            return None
        if idx is None:
            rpos = self.pos
        else:
            rpos = self.pos[idx, :]
        return to_cartesian(rpos)

    def get_antenna_param(self, idx=None, param=None):
        if self.pos is None:
            return None
        if type(idx) is int or type(idx) is np.int64:
            return self.antennas[idx][param]
        if idx is None:
            return np.array([a[param] for a in self.antennas])
        else:
            return np.array([self.antennas[i][param] for i in idx])
