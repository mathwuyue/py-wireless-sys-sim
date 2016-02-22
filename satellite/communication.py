import numpy as np
from core import gen_rician, gen_logNshadowing, cal_dist_3d, to_cartesian, \
    cal_fiirs, cal_recv_power, cal_thermal_noise, cal_shannon_cap


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

    def update_pos(self, *args):
        self.satellites.update_pos(*args)

    def intra_comm(self, start, dest, comm_t=INTRA_COMM):
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
                tr_pos = np.append(tr_pos, s.satellite_pos(idx), axis=0)
                f = np.append(f, s.get_antenna_param(idx, 'intra_f'))
                bw = np.append(bw, s.intra_bw*np.ones(len(idx)))
                tp = np.append(tp, s.get_antenna_param(idx, 'max_tp'))
                gt = np.append(gt, s.get_antenna_param(idx, 'gain'))
            elif comm_t == DOWNLINK:
                s = self.satellites[key]
                tr_pos = np.append(tr_pos, s.satellite_pos(idx), axis=0)
                f = np.append(f, s.get_antenna_param(idx, 'earth_f'))
                bw = np.append(bw, s.earth_bw*np.ones(len(idx)))
                tp = np.append(tp, s.get_antenna_param(idx, 'max_tp'))
                gt = np.append(gt, s.get_antenna_param(idx, 'gain'))
            else:
                s = self.satellites[key].stations
                tr_pos = np.append(tr_pos, s.station_pos(idx), axis=0)
                f = np.append(f, s.get_antenna_param(idx, 'f'))
                bw = np.append(bw, s.bw*np.ones(len(idx)))
                tp = np.append(tp, s.get_antenna_param(idx, 'max_tp'))
                gt = np.append(gt, s.get_antenna_param(idx, 'gain'))
            n = n + len(idx)
        for key, idx in dest.iteritems():
            if comm_t == INTRA_COMM or comm_t == UPLINK:
                e = self.satellites[key]
                rv_pos = np.append(rv_pos, e.satellite_pos(idx), axis=0)
                gr = np.append(gr, e.get_antenna_param(idx, 'gain'))
            else:
                e = self.satellites[key].stations
                rv_pos = np.append(rv_pos, e.station_pos(idx), axis=0)
                gr = np.append(gr, e.get_antenna_param(idx, 'gain'))
        tr_pos = tr_pos[1:, :]
        rv_pos = rv_pos[1:, :]
        #print tr_pos, rv_pos
        rp = cal_recv_power(tr_pos, rv_pos, 10**(tp/10), n, 1,
                            dist_func=cal_dist_3d, dist_args=[],
                            pl_func=cal_fiirs, pl_args=[f, gt, gr],
                            fading_func=gen_rician, fading_args=[10, 1],
                            shadowing_func=gen_logNshadowing, shadowing_args=[4])
        if comm_t == INTRA_COMM or comm_t == UPLINK:
            noise = cal_thermal_noise(bw*1e6, 2.73)
        else:
            noise = cal_thermal_noise(bw*1e6, 290)
        throughput = cal_shannon_cap(bw*1e6, rp, 0, noise=noise)
        return throughput

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
        ue_pos_cart = np.kron(ue_pos_cart, np.ones((s.n, 1)))
        if len(ue_pos.shape) != 1:
            n_ues = ue_pos.shape[0]
            s_pos = np.kron(np.ones((n_ues, 1)), s_pos)
            dist = cal_dist_3d(s_pos, ue_pos_cart)
            return np.argmin(np.reshape(dist, (n_ues, s.n)), axis=1)
        else:
            return np.argmin(s_pos-ue_pos_cart)

    def comm_ue(self, satellite, ue, comm_t):
        """
        Args:
        satellite (list of tuple): {'s1': [1,2,..])
        ue (numpy array): ([r,i,theta], tp)
        comm_t: UPLINK or DOWNLINK

        Return:
        throughput (numpy array): numpy.array([1,2,3,4...])

        .. note::
        TODO: multiple UEs and satellite
        """
        # ss_idx: satellite system idx. s_idx: satellite idx in system
        for ss_idx, s_idx in satellite.iteritems():
            s_pos = self.satellites[ss_idx].satellite_pos(s_idx)
        bw = self.satellites[ss_idx].earth_bw
        f = self.satellites[ss_idx].get_antenna_param(s_idx, 'earth_f')
        if comm_t == UPLINK:
            noise = cal_thermal_noise(bw*1e6, 2.73)
            tp = ue.tp
            tr_pos = ue.pos
            rv_pos = s_pos
            gt = 0
            gr = self.satellites[ss_idx].get_antenna_param(s_idx, 'gain')
        else:
            noise = cal_thermal_noise(bw*1e6, 290.0)
            tp = self.satellites[ss_idx].get_antenna_param(s_idx, 'max_tp')
            tr_pos = s_pos
            rv_pos = ue.pos
            gt = self.satellites[ss_idx].get_antenna_param(s_idx, 'gain')
            gr = 0
        rp = cal_recv_power(tr_pos, rv_pos, 10**(tp/10.0), 1, 1,
                            dist_func=cal_dist_3d, dist_args=[],
                            pl_func=cal_fiirs, pl_args=[f, gt, gr],
                            fading_func=gen_rician, fading_args=[10.0, 1.0],
                            shadowing_func=gen_logNshadowing, shadowing_args=[4.0])
        throughput = cal_shannon_cap(bw*1e6, rp, ip=0, noise=noise)
        return throughput, rp
