import numpy as np
from satellite.communication import SatelliteComm, INTRA_COMM, UPLINK, DOWNLINK
from satellite.satellite import IridiumSatellite, GeoSatellite, EarthStation
from core.ue import UE
import matplotlib.pyplot as plt
import random

tcp_size = 1500
leos = IridiumSatellite()
earth_station = EarthStation(pos=np.array([6376e3, np.pi/2, np.pi/2]))
geos = GeoSatellite(3, stations=earth_station)
s_comm = SatelliteComm({'leo': leos, 'geo': geos})
ue1 = UE(pos=np.array([6376e3, np.pi/2, 0]), tp=14)
ue2 = [UE(pos=np.array([6376e3, np.pi/2, i*np.pi/6]), tp=14)
       for i in range(1, 7)]

leo_topo = leos.get_satellite_topo()


def geo_sim():
    t = []
    ss = s_comm.choose_satellite('geo', np.array([ue1.pos] +
                                                 [ue.pos for ue in ue2]))
    s1 = ss[0]
    for i in range(1, 7):
        s2 = ss[i]
        t1 = s_comm.comm_ue({'geo': s1}, ue1, UPLINK)
        if s1 != s2:
            t2 = s_comm.intra_comm({'geo': [s1]}, {'geo': [s2]}, INTRA_COMM)
        else:
            t2 = np.inf
        t3 = s_comm.comm_ue({'geo': s2}, ue2[i-1], DOWNLINK)
        t.append(min(t1, t2, t3))
    return np.array(t)


def leo_sim():
    t = []
    ss = s_comm.choose_satellite('leo', np.array([ue1.pos] +
                                                 [ue.pos for ue in ue2]))
    s1 = ss[0]
    for i in range(1, 7):
        s2 = ss[i]
        t1 = s_comm.comm_ue({'leo': s1}, ue1, UPLINK)
        # random walk
        if s2 != s1:
            path = []
            sc = s1
            s1_orbit_idx = s1 / 11
            s2_orbit_idx = s2 / 11
            diff_orbit_idx = abs(s1_orbit_idx - s2_orbit_idx)
            if s1_orbit_idx > s2_orbit_idx:
                for j in range(1, diff_orbit_idx+1):
                    path.append(s1-11*j)
                sc = s1 - diff_orbit_idx*11
            elif s1_orbit_idx < s2_orbit_idx:
                for j in range(1, diff_orbit_idx+1):
                    path.append(s1+11*j)
                sc = s1 + diff_orbit_idx*11
            if sc < s2:
                for j in range(1, s2-sc):
                    path.append(sc+j)
            elif sc > s2:
                for j in range(1, sc-s2):
                    path.append(sc-j)
            start = [s1] + path
            end = path + [s2]
            t2 = min(s_comm.intra_comm({'leo': start}, {'leo': end}, INTRA_COMM))
        t3 = s_comm.comm_ue({'leo': s2}, ue2[i-1], DOWNLINK)
        t.append(min(t1, t2, t3))
    return np.array(t)


def get_satellite_comm_topo(topo):
    pass


def leo_sdn_sim(n_packets):
    leo_comm_topo = get_satellite_comm_topo(leo_topo)
    pass


def main():
    geo_throughput = np.zeros(6)
    leo_throughput = np.zeros(6)
    leo_sdn_throughput = np.zeros(6)
    loop = 10000

    for i in xrange(loop):
        geo_throughput = geo_throughput + geo_sim()
        leo_throughput = leo_throughput + leo_sim()
        if i % 1000 == 0:
            print i
#        leo_sdn_throughput = leo_sdn_throughput + leo_sdn_sim()
    x = [i * np.pi/6 for i in range(6)]
    plt.plot(x, geo_throughput/loop)
    plt.plot(x, leo_throughput/loop)
#    plt.plot(x, leo_sdn_throughput/loop)


if __name__ == '__main__':
    main()
