import numpy as np
from core.statchannel import gen_rician, gen_logNshadowing
from core.position import cal_dist_3d
from core.pathloss import cal_fiirs
from core.communicaiton import cal_recv_power, cal_thermal_noise, cal_shannon_cap


INTRA_COMM = 0
UPLINK = 1
DOWNLINK = 2


class SatelliteComm(object):
    """
    This is the class for satellite communicaitons
    """
    def __init__(self, satellite_system):
        """
        Args:
        satellite_system (dict): dict contains all the satellite systems in the simulation
        """
        self.satellites = satellite_system
        self.earth_stations = {}
        for key, ss in self.satellites.iteritems():
            self.earth_stations[key] = ss.stations

    def update_pos(self, t):
        pass

    def comm(self, start, dest, n, comm_t=INTRA_COMM):
        """
        Args:
        start (dict): {system_key: idx_list (list, None for all)},
        end (dict): {system_key: idx_list (list, None for all)}
        .. note::
        The elements in start and end should be consistent, e.g. all elements in start/end are satellites or earth stations.

        Returns:
        Shannon capacity (numpy array): array of the shannon capacity between all trs and rvs.
        """
        tr_pos = np.array([[0, 0, 0]])
        rv_pos = np.array([[0, 0, 0]])
        f = np.array([])
        gt = np.array([])
        gr = np.array([])
        bw = np.array([])
        n = 0
        for key, idx in start.iteritems():
            if comm_t == INTRA_COMM:
                s = self.satellites[key]
                np.append(f, s.get_antenna_param(idx, 'intra_f'))
                np.append(bw, s.intra_bw*np.ones(len(idx)))
            elif comm_t == DOWNLINK:
                s = self.satellites[key]
                np.append(f, s.get_antenna_param(idx, 'earth_f'))
                np.append(bw, s.earth_bw*np.ones(len(idx)))
            else:
                s = self.satellites[key].stations
                np.append(f, s.get_antenna_param(idx, 'f'))
                np.append(bw, s.stations.bw*np.ones(len(idx)))
            np.append(tr_pos, s.to_cartesian(idx))
            np.append(gt, s.get_antenna_param(idx, 'gain'))
            n = n + len(idx)
        for key, idx in dest.iteritems():
            if comm_t == INTRA_COMM or UPLINK:
                e = self.satellites[key]
            else:
                e = self.satellites[key].stations
            np.append(rv_pos, e.to_cartesian(idx))
            np.append(gr, e.get_antenna_param(idx, 'gain'))
        tr_pos = tr_pos[1:, :]
        rv_pos = rv_pos[1:, :]
        rp = cal_recv_power(tr_pos, rv_pos, n, 1,
                            dist_func=cal_dist_3d,
                            pl_func=cal_fiirs, pl_args=[f, gt, gr],
                            fading_func=gen_rician, fading_args=[10, 1],
                            shadowing_func=gen_logNshadowing, shadowing_args=[4])
        if comm_t == INTRA_COMM or comm_t == UPLINK:
            noise = cal_thermal_noise(bw*1e6, 2.73)
        else:
            noise = cal_thermal_noise(bw*1e6, 290)
        throughput = cal_shannon_cap(bw*1e6, rp, noise=noise)
        return rp, throughput
