import numpy as np
from satellite.communication import SatelliteComm, INTRA_COMM, UPLINK, DOWNLINK
from satellite.satellite import IridiumSatellite, GeoSatellite, EarthStation
from core.ue import UE

tcp_size = 1500
leos = IridiumSatellite()
earth_station = EarthStation(pos=np.array([6376e3, np.pi/2, np.pi/2]))
geos = GeoSatellite(3, stations=earth_station)
s_comm = SatelliteComm({'leo': leos, 'geo': geos})
ue1 = UE(pos=np.array([6376e3, np.pi/3, np.pi/6]), tp=14)
ue2 = UE(pos=np.array([6376e3, np.pi/2, 7*np.pi/6]), tp=14)


def geo_sim(n_packets):
    t = 0
    s1, s2 = s_comm.choose_satellite('geo', np.array([ue1.pos, ue2.pos]))
    for i in range(n_packets):
        t1 = tcp_size*8 / s_comm.comm_ue({'geo': s1}, ue1, UPLINK)
        if s1 != s2:
            t2 = tcp_size*8 / s_comm.intra_comm({'geo': [s1]}, {'geo': [s2]}, INTRA_COMM)
        else:
            t2 = 0
        t3 = tcp_size*8 / s_comm.comm_ue({'geo': s2}, ue2, DOWNLINK)
        t = t1 + t2[0] + t3
    return t


def main():
    print geo_sim(10)


if __name__ == '__main__':
    main()
