import numpy as np
import csv
import matplotlib.pyplot as plt
from satellite import Satellite, SatelliteComm, GeoSatellite
from core import EarthStation, SatelliteAntenna, UE, to_cartesian

G = 6.67e-11
M_EARTH = 5.97e24
DELTA_T = 1e-3
HEIGHT = (np.linspace(160, 2000, 8) + 6371) * 1e3
T_REGION = [600, 100, 0]
UE_LOC = [0, 45, 100]
AV = np.sqrt(G*M_EARTH/(HEIGHT**3))
LOOP = 10000
R_RE = 1500*8.0
SDN_PACKET = 1500*8.0
H_IP = 590*8.0
COMM_THRESHOLD = 6e-18
T_RETRY = 0.25


class LEOSystem(Satellite):
    def __init__(self, height, av, max_tp):
        super(LEOSystem, self).__init__(n=2)
        self.height = height
        self.av = av
        self.pos = np.array([[height, np.pi/2, np.pi/2-i*np.pi/6] for i in range(self.n)])
        self.antennas = [SatelliteAntenna(27.5, 6, 30, 0, max_tp) for i in range(self.n)]

    def update_pos(self, t):
        return self.pos + np.array([[0, 0, self.av*t] for i in range(self.n)])


class LEOUE(UE):
    def __init__(self, user_loc):
        super(LEOUE, self).__init__(np.array([6371, np.pi/2, np.pi/2]))
        self.user_loc = user_loc

    def user_pos(self):
        return to_cartesian(self.pos)

    def __getitem__(self, key):
        return self.__dict__[key]


class SDSN(object):
    def __init__(self, i=0, t=0, user_loc=0):
        self.earth_station = EarthStation(np.array([[6371e3, np.pi/2, j*2*np.pi/3]
                                                    for j in range(3)]))
        self.satellites = {'l': LEOSystem(HEIGHT[i], AV[i], 20),
                           'g': GeoSatellite(stations=self.earth_station)}
        self.s_comm = SatelliteComm(self.satellites)
        self.ue = LEOUE(user_loc)
        self.is_hard_success = False
        self.is_sdsn_success = False
        self.is_hard_first = False
        self.is_sdsn_first = False
        self.is_out_of_region = False
        self.sdsn_latency = 0
        self.hard_latency = 0
        self.handover_t = 0
        self.region_t = t
        self.angle = np.arctan(user_loc/self.satellites['l'].height)

    def begin_handover(self):
        th = (np.pi/12) / self.satellites['l'].av
        self.satellites['l'].update_pos(th+self.handover_t)

    def handover_process(self):
        # hard handover
        throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, self.angle, 1)
        # first try
        if rp > COMM_THRESHOLD:
            self.is_hard_success = True
            self.is_hard_first = True
            self._cal_hard_latency()
        # second try
        else:
            if self.hard_latency + T_RETRY <= self.region_t:
                self.hard_latency = self.hard_latency + T_RETRY
                self.satellites['l'].update_pos(T_RETRY)
                throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, self.angle, 1)
                if rp > COMM_THRESHOLD:
                    self.is_hard_success = True
                    self._cal_hard_latency()
        # sdsn handover
        if not self.is_hard_first:
            self.satellites['l'].update_pos(-T_RETRY)
        # first try
        throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, self.angle, 1)
        if rp > COMM_THRESHOLD:
            self.is_sdsn_first = True
            self.is_sdsn_success = True
            self._cal_sdsn_latency()
        else:
            while True:
                throughput, rp = self.s_comm.comm_ue({'l': 0}, self.ue, self.angle, 1)
                if rp < COMM_THRESHOLD:
                    break
                if self.sdsn_latency+T_RETRY > self.region_t:
                    self.is_out_of_region = True
                    break
                self.satellites['l'].update_pos(T_RETRY)
                self.sdsn_latency = self.sdsn_latency + T_RETRY
                throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, self.angle, 1)
                if rp > COMM_THRESHOLD:
                    self.is_sdsn_first = True
                    self.is_sdsn_success = True
                    self._cal_sdsn_latency()
                    break
            # second try
            if not self.is_sdsn_success and not self.is_out_of_region:
                self.sdsn_latency = self.sdsn_latency + T_RETRY
                self.satellites['l'].update_pos(T_RETRY)
                throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, self.angle, 1)
                if rp > COMM_THRESHOLD:
                    self.is_sdsn_success = True
                    self._cal_hard_latency()

    # check handover status
    def get_handover_status(self):
        return self.is_hard_success,\
            self.is_hard_first,\
            self.is_sdsn_success,\
            self.is_sdsn_first

    def get_latency(self):
        return self.hard_latency, self.sdsn_latency

    def _cal_hard_latency(self):
        for i in range(2):
            throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, self.angle, 1)
            self.hard_latency = self.hard_latency + H_IP / throughput
            throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, self.angle, 2)
            self.hard_latency = self.hard_latency + H_IP / throughput

    def _cal_sdsn_latency(self):
        # ue -> leo
        throughput, rp = self.s_comm.comm_ue({'l': 1}, self.ue, self.angle, 1)
        self.sdsn_latency = self.sdsn_latency + SDN_PACKET/throughput
        # leo -> geo
        throughput = self.s_comm.intra_comm(start={'l': [1]}, dest={'g': [1]})
        self.sdsn_latency = self.sdsn_latency + SDN_PACKET/throughput
        # geo -> station
        throughput = self.s_comm.intra_comm(start={'g': [1]}, dest={'g': [1]}, comm_t=2)
        self.sdsn_latency = self.sdsn_latency + SDN_PACKET/throughput
        # station -> geo
        throughput = self.s_comm.intra_comm(start={'g': [1]}, dest={'g': [1]}, comm_t=1)
        self.sdsn_latency = self.sdsn_latency + SDN_PACKET/throughput
        # geo -> leo
        throughput = self.s_comm.intra_comm(start={'g': [1]}, dest={'l': [1]})
        self.sdsn_latency = self.sdsn_latency + SDN_PACKET/throughput


def main():
    hard_failed_times = np.zeros([3, 8])
    sdsn_failed_times = np.zeros([3, 8])
    hard_first_times = np.zeros([3, 8])
    sdsn_first_times = np.zeros([3, 8])
    hard_latency = np.zeros([3, 8])
    sdsn_latency = np.zeros([3, 8])
    i_hard = np.zeros([3, 8])
    i_sdsn = np.zeros([3, 8])
    for k in range(3):
        for i in range(8):
            for j in range(LOOP):
                sdsn_sim = SDSN(i, T_REGION[k], UE_LOC[k])
                sdsn_sim.begin_handover()
                sdsn_sim.handover_process()
                h, hf, s, sf = sdsn_sim.get_handover_status()
                hl, sl = sdsn_sim.get_latency()
                if not h:
                    hard_failed_times[k, i] = hard_failed_times[k, i] + 1
                else:
                    hard_latency[k, i] = hard_latency[k, i] + hl
                    i_hard[k, i] = i_hard[k, i] + 1
                if not hf:
                    hard_first_times[k, i] = hard_first_times[k, i] + 1
                if not s:
                    sdsn_failed_times[k, i] = sdsn_failed_times[k, i] + 1
                else:
                    sdsn_latency[k, i] = sdsn_latency[k, i] + sl
                    i_sdsn[k, i] = i_sdsn[k, i] + 1
                if not sf:
                    sdsn_first_times[k, i] = sdsn_first_times[k, i] + 1

    print hard_failed_times
    print sdsn_failed_times
    print hard_first_times
    print sdsn_first_times
    print i_hard
    print i_sdsn

    hard_mean_latency = hard_latency / i_hard * 1000.0
    sdsn_mean_latency = sdsn_latency / i_sdsn * 1000.0
    hard_score = 1-np.exp(-1.0/hard_mean_latency)-0.1*hard_first_times/10000-0.1*hard_failed_times/10000
    sdsn_score = 1-np.exp(-1.0/sdsn_mean_latency)-0.1*sdsn_first_times/10000-0.1*sdsn_failed_times/10000

    with open('score.csv', 'wb') as f:
        dw = csv.writer(f, delimiter=' ')
        for i in range(3):
            dw.writerow(hard_score[i, :])
        for i in range(3):
            dw.writerow(sdsn_score[i, :])
    with open('latency.csv', 'wb') as f:
        dw = csv.writer(f, delimiter=' ')
        for i in range(3):
            dw.writerow(hard_mean_latency[i, :])
        for i in range(3):
            dw.writerow(sdsn_mean_latency[i, :])
    with open('complete_failed.csv', 'wb') as f:
        dw = csv.writer(f, delimiter=' ')
        for i in range(3):
            dw.writerow(hard_failed_times[i, :])
        for i in range(3):
            dw.writerow(sdsn_failed_times[i, :])
    with open('first_failed.csv', 'wb') as f:
        dw = csv.writer(f, delimiter=' ')
        for i in range(3):
            dw.writerow(hard_first_times[i, :])
        for i in range(3):
            dw.writerow(sdsn_first_times[i, :])

    x = np.linspace(160, 2000, 8)
    plt.plot(x, hard_mean_latency[0, :], '-o', linewidth=2, label='Hard, t=0')
    plt.plot(x, sdsn_mean_latency[0, :], '--o', linewidth=2, label='SDSN, t=0')
    plt.plot(x, hard_mean_latency[1, :], '-x', linewidth=2, label='Hard, t=250')
    plt.plot(x, sdsn_mean_latency[1, :], '--x', linewidth=2, label='SDSN, t=250')
    plt.plot(x, hard_mean_latency[2, :], '-^', linewidth=2, label='Hard, t=500')
    plt.plot(x, sdsn_mean_latency[2, :], '--^', linewidth=2, label='SDSN, t=500')
    plt.legend()
    plt.savefig('sdsn_latency.eps', format='eps')
    plt.show()

if __name__ == '__main__':
    main()
