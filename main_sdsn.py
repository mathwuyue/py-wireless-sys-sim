import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from satellite import Satellite, SatelliteComm, GeoSatellite
from core import EarthStation, SatelliteAntenna, UE, cal_dist_3d, to_cartesian

G = 6.67e-11
M_EARTH = 5.97e24
DELTA_T = 1e-3
HEIGHT = (np.linspace(160, 2000, 8) + 6371) * 1e3
AV = np.sqrt(G*M_EARTH/(HEIGHT**3))
LOOP = 1
PB = 0.127
PB_P = 0.08
R_H = 7800*8
R_RE = 15000*8
R_H_SDN = 1500*8
COMM_THRESHOLD = 6e-18


class LEOSystem(Satellite):
    def __init__(self, height, av):
        super(LEOSystem, self).__init__(n=2)
        self.height = height
        self.av = av
        self.pos = np.array([[height, np.pi/2, np.pi/2+i*np.pi/6] for i in range(self.n)])
        self.antennas = [SatelliteAntenna(27.5, 6, 30, 0, 4.2) for i in range(self.n)]

    def update_pos(self, t):
        return self.pos + np.array([[0, 0, self.av*t] for i in range(self.n)])


class LEOUE(UE):
    def __init__(self):
        super(LEOUE, self).__init__(np.array([6371, np.pi/2, np.pi/2]))

    def user_pos(self):
        return to_cartesian(self.pos)


class SDSN(object):
    def __init__(self, i=0, t=0):
        self.earth_station = EarthStation(np.array([[6371, np.pi/2, i*2*np.pi/3]
                                                    for j in range(3)]))
        self.satellites = {'l': LEOSystem(HEIGHT[i], AV[i]),
                           'g': GeoSatellite(stations=self.earth_station)}
        self.s_comm = SatelliteComm(self.satellites)
        self.ue = LEOUE()
        self.is_handover = False
        self.is_success = False
        self.sdsn_latency = 0
        self.no_latency = 0
        self.handover_t = t

    def begin_handover(self):
        th = (np.pi/12) / self.satellites['l'].av
        self.satellites['l'].update_pos(th+self.handover_t)

    def handover_process(self):
        # hard handover
        throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, 1)
        if rp > COMM_THRESHOLD:
            self.is_success = True
            return self.is_success

    # check handover status
    def check_handover(self):
        if self.is_handover and self.is_success:
            self.is_handover = False
        return self.is_handover


def main():
    failed_times = np.zeros(6)
    for i in range(6):
        sdsn_sim = SDSN(i)
        sdsn_sim.begin_handover()
        for j in range(10000):
            if not sdsn_sim.handover_process():
                failed_times[i] = failed_times[i] + 1
    print failed_times

if __name__ == '__main__':
    main()
