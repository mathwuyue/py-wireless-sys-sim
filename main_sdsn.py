import numpy as np
from core import gen_rician, gen_logNshadowing, cal_dist_2d, cal_thermal_noise,\
    cal_recv_power, SatelliteAntenna
from satellite import Satellite

G = 6.67e-11
M_EARTH = 5.97e24
DELTA_T = 1e-3
HEIGHT = (np.linspace(160, 2000, 8) + 6371) * 1e3
V = np.sqrt(G*M_EARTH/HEIGHT)
LOOP = 1
PB = 0.127
PB_P = 0.08


class LEOSystem(Satellite):
    def __init__(self, height, n=2, intra_bw=1, earth_bw=1):
        super(LEOSystem, self).__init__(n)
        delta_s = 2 * np.pi * height / 12.0
        self.pos = np.array([-i*delta_s + height*1j for i in range(n)])
        self.antennas = [SatelliteAntenna(27.5, 6, 30, 0, 4.2) for i in range(n)]

    def update_pos(self, v, t):
        return self.pos + v * t


# check handover condition
def check_handover(p1, p2, puser):
    return cal_dist_2d(p2, puser) < cal_dist_2d(p1, puser)


# check whether handover successful
def is_handover_finish(porbit, puser):
    noise = cal_thermal_noise(1e6, 2.73)
    recv_power = cal_recv_power()
