import numpy as np
from position import to_cartesian
from antenna import ParabolicAntenna


class Station(object):
    def __init__(self, pos, bw, tp, n_channel=1, height=0, antenna=None):
        self.bw = bw
        self.pos = pos,
        self.height = height
        self.n_channel = n_channel
        if n_channel > 1 and (type(tp) is int or tp.shape == ()):
            self.tp = tp * np.ones(n_channel)
        else:
            self.tp = tp
        self.antenna = antenna

    def cal_antenna_gain(self, ue):
        pass

    def get_antenna_param(self, idx=None, param=None):
        if self.pos is None:
            return None
        if type(idx) is int or type(idx) is np.int64:
            return self.antennas[idx][param]
        if idx is None:
            return np.array([a[param] for a in self.antennas])
        else:
            return np.array([self.antennas[i][param] for i in idx])


class LTEeNB(Station):
    def __init__(self, pos, bw=10, n_channel=1, tp=16, height=25, antenna=None):
        super(LTEeNB, self).__init__(pos, bw, tp, n_channel, height, antenna)

    def cal_antenna_gain(self, device):
        phi = np.angle(device.pos)
        theta = np.arctan((self.height-device.height) / np.abs(device.pos))
        return self.antenna.cal_gain(phi, theta)


class EarthStation(Station):
    def __init__(self, pos=None, bw=1, tp=30, n_channel=1, height=0, antennas=None):
        """
        Kwargs:
        bw (float): MHz
        pos (numpy array): (r, i, theta)
        antennas (list)
        .. note::
        pos have same rows as antennas, except that pos is not None but antennas is None, in which case antennas will be set to the array of default parabolicantennas. 
        """
        super(EarthStation, self).__init__(pos, bw, tp, n_channel, height, antennas)
        if pos is not None and antennas is None:
            self.antennas = [ParabolicAntenna() for p in self.pos]
        else:
            self.antennas = antennas

    def station_pos(self, idx=None):
        if self.pos is None:
            return None
        if idx is None:
            rpos = self.pos
        else:
            rpos = self.pos[idx, :]
        return to_cartesian(rpos)
