import numpy as np


class Antenna(object):
    """
    A class to decribe the antenna system
    """
    def __init__(self, f=0, gain=0, theta_3db=0, max_tp=0):
        self.f = f
        self.gain = gain
        self.theta_3db = theta_3db
        self.max_tp = 0

    def cal_gain(self):
        return self.gain


class ParabolicAntenna(Antenna):
    def __init__(self, f=6.0, d=3.0, k=70, ea=0.6, max_tp=20):
        """
        Args:
        f (float): in GHz,
        d (float): in m,
        ea (float): efficiency
        max_tp (float): dBW
        """
        super(ParabolicAntenna, self).__init__(f, max_tp=max_tp)
        self.c = 3.0e8
        self.d = d
        self.gain = (np.pi*d*f*1.0e9/self.c)**2 * ea
        self.theta_3db = (k*self.c/(f*1e9*d))/d

    def cal_gain(self, alpha_t):
        return self.gain - 0.00245*(alpha_t*self.d*self.f*1e9/self.c)**2


class SatelliteAntenna(Antenna):
    def __init__(self, intra_f, earth_f, gain, theta_3db, max_tp):
        super(SatelliteAntenna, self).__init__(gain=gain, theta_3db=theta_3db, max_tp=max_tp)
        self.intra_f = intra_f
        self.earth_f = earth_f
