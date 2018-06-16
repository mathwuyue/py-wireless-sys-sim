import numpy as np
import numpy.matlib as npmt
from d2d.model import D2DSystemModel
from core.communication import cal_recv_power, cal_shannon_cap, cal_thermal_noise
from core.pathloss import cal_umi_nlos
from core.statchannel import gen_rayleigh, gen_logNshadowing
from core.position import gen_uni_circ_pos, cal_dist_2d


class Rui2016(D2DSystemModel):
    """Benchmark from Rui2016 paper"""
    def __init__(self, n_cc, n_pairs, n_rb=50, ue_tp=-10, total_bw = 10e6, cell_r=100, d2d_r=20):
        super(Rui2016, self).__init__(n_cc, n_pairs, ue_tp, total_bw, cell_r, d2d_r)
        self.d2d_tps = ue_tp
        self.cc_tps = ue_tp
        self.n_rb = n_rb

    def gen_cc_ues(self, n_cc=None):
        if not n_cc:
            n_cc = self.n_cc
        super(Rui2016, self).gen_cc_ues(n_cc)

    def gen_d2d_pairs(self, n_pairs=None):
        if not n_pairs:
            n_pairs = self.n_pairs
        super(Rui2016, self).gen_d2d_pairs(n_pairs)

    def cal_interference(self, cc_tps, d2d_tps, alpha=3.5):
        # cal cc interference
        cc_inter_tr = np.array([d.pos for d in self.d2d_tr])
        cc_inter_rv = 0
        cc_inter_tp = d2d_tps
        c_ip = cal_recv_power(cc_inter_tr, cc_inter_rv, d2d_tps, self.n_pairs, self.n_rb,
                              cal_dist_2d, [],
                              cal_umi_nlos, [alpha, self.d2d_tr[0].freq],
                              fading_func=gen_rayleigh, fading_args=[1.0],
                              shadowing_func=gen_logNshadowing, shadowing_args=[4])
        cc_interference = np.sum(c_ip, axis=0)
        # cal d2d interference
        # cc -> d2d
        cc_tr = np.array([c.pos for c in self.cc_ue])
        d2d_rv = np.array([d.pos for d in self.d2d_rc])
        d2d_cc_tr = np.kron(cc_tr, ones(self.n_rb / self.n_cc, 1))
        d2d_cc_tr = npmt.repmat(d2d_cc_tr, self.n_pairs, 1)
        d2d_cc_rv = npmt.repmat(np.reshape(d2d_rv, self.n_pairs, 1), 1, self.n_rb)
        cc_d2d_ip = cal_recv_power(d2d_cc_tr, d2d_cc_rv, self.n_pairs, self.n_rb,
                                   cal_dist_2d, [],
                                   cal_umi_nlos, [alpha, self.cc_ue[0].freq],
                                   fading_func=gen_rayleigh, fading_args=[1.0],
                                   shadowing_func=gen_logNshadowing, shadowing_args=[4])
        # other d2d -> d2d
        d2d_d2d_ip = []
        for idx, d in enumerate(d2d_rv):
            d2d_inter_tr = [d.pos for d in self.d2d_tr]
            del d2d_inter_tr[idx]
            d2d_inter_tr = np.array(d2d_inter_tr)
            d_d2d_ip = cal_recv_power(d2d_inter_tr, d, self.n_pairs-1, self.n_rb,
                                      cal_dist_2d, [],
                                      cal_umi_nlos, [alpha, self.cc_ue[0].freq],
                                      fading_func=gen_rayleigh, fading_args=[1.0],
                                      shadowing_func=gen_logNshadowing, shadowing_args=[4])
            d2d_d2d_ip.append(d_d2d_ip)
        d2d_interference = cc_d2d_ip + np.array(d2d_d2d_ip)
        return cc_interference, d2d_interference

    def cal_sinr(self):
