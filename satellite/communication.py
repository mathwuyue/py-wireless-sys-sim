import numpy as np
from core.statchannel import gen_rician, gen_logNshadowing
from core.position import cal_dist_3d, to_cartesian
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
        tp = np.array([])
        n = 0
        for key, idx in start.iteritems():
            if comm_t == INTRA_COMM:
                s = self.satellites[key]
                np.append(f, s.get_antenna_param(idx, 'intra_f'))
                np.append(bw, s.intra_bw*np.ones(len(idx)))
                np.append(tp, s.get_antenna_param(idx, 'max_tp'))
            elif comm_t == DOWNLINK:
                s = self.satellites[key]
                np.append(f, s.get_antenna_param(idx, 'earth_f'))
                np.append(bw, s.earth_bw*np.ones(len(idx)))
                np.append(tp, s.get_antenna_param(idx, 'max_tp'))
            else:
                s = self.satellites[key].stations
                np.append(f, s.stations.get_antenna_param(idx, 'f'))
                np.append(bw, s.stations.bw*np.ones(len(idx)))
                np.append(tp, s.stations.get_antenna_param(idx, 'max_tp'))
            np.append(tr_pos, s.satellite_pos(idx))
            np.append(gt, s.get_antenna_param(idx, 'gain'))
            n = n + len(idx)
        for key, idx in dest.iteritems():
            if comm_t == INTRA_COMM or UPLINK:
                e = self.satellites[key]
            else:
                e = self.satellites[key].stations
            np.append(rv_pos, e.satellite_pos(idx))
            np.append(gr, e.get_antenna_param(idx, 'gain'))
        tr_pos = tr_pos[1:, :]
        rv_pos = rv_pos[1:, :]
        rp = cal_recv_power(tr_pos, rv_pos, 10**(tp/10), n, 1,
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

    def choose_satellite(self, ss_idx, ue_pos):
        """
        Args:
        ss_idx (string): idx of the satellite system
        ue_pos (numpy array): (r, i, theta)...

        Return:
        Satellites (list of tuple): (1, 2, 3...)
        """
        s = self.satellites[ss_idx]
        s_pos = s.satellite_pos()
        ue_pos_cart = to_cartesian(ue_pos)
        if len(ue_pos.shape) != 1:
            n_ues = ue_pos.shape[0]
            s_pos = np.kron(np.ones((n_ues, 1)), s_pos)
        ue_pos_cart = np.kron(ue_pos_cart, np.ones((s.n, 1)))
        return np.argmin(s_pos - ue_pos_cart)

    def comm_ue(self, satellite, ue, comm_t):
        """
        Args:
        satellite (list of tuple): (s1, 1)
        ue_pos (numpy array): ([r,i,theta], tp)
        comm_t: UPLINK or DOWNLINK

        Return:
        throughput (numpy array): numpy.array([1,2,3,4...])
        """
        # ss_idx: satellite system idx. s_idx: satellite idx in system
        ss_idx, s_idx = satellite
        s_pos = self.satellites[ss_idx].satellite_pos(s_idx)
        bw = self.satellites[ss_idx].earth_bw
        f = self.satellites[ss_idx].get_antenna_param(s_idx, 'earth_f')
        if comm_t == UPLINK:
            noise = cal_thermal_noise(bw, 2.73)
            tp = ue.tp
            tr_pos = ue.pos
            rv_pos = s_pos
            gt = 0
            gr = self.satellites[ss_idx].get_antenna_param(s_idx, 'gain')
        else:
            noise = cal_thermal_noise(bw, 290)
            tp = self.satellites[ss_idx].get_antenna_param(s_idx, 'max_tp')
            tr_pos = s_pos
            rv_pos = ue.pos
            gt = self.satellites[ss_idx].get_antenna_param(s_idx, 'gain')
            gr = 0
        rp = cal_recv_power(tr_pos, rv_pos, 10**(tp/10), 1, 1,
                            dist_func=cal_dist_3d,
                            pl_func=cal_fiirs, pl_args=[f, gt, gr],
                            fading_func=gen_rician, fading_args=[10, 1],
                            shadowing_func=gen_logNshadowing, shadowing_args=[4])
        throughput = cal_shannon_cap(bw, rp, noise=noise)
        return rp, throughput
