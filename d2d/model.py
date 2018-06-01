from core.position import gen_uni_circ_pos, cal_dist_2d
from core.device import UE, BS
from core.communication import cal_recv_power, cal_shannon_cap, cal_thermal_noise
from core.pathloss import cal_umi_nlos
from core.statchannel import gen_rayleigh, gen_logNshadowing
import numpy as np

UL = 0
DL = 1


class D2DSystemModel(object):
    def __init__(self, n_cc, n_pairs, ue_tp=-10, total_bw=10e6, cell_r=250, d2d_r=20):
        self.n_cc = n_cc
        self.n_pairs = n_pairs
        self.gen_cc_ues(n_cc)
        self.gen_d2d_pairs(n_pairs)
        self.bs = BS(0)
        self.total_bw = total_bw
        self.cell_r = cell_r
        self.d2d_r = d2d_r
        self.ue_tp = ue_tp
        self.cc_ue = None
        self.d2d_tr = None
        self.d2d_rc = None

    def __setitem__(self, key, value):
        self.__dict__[key] = value
        if key == 'd2d_r':
            self.gen_d2d_pairs(self.n_pairs)
        elif key == 'cell_r':
            self.gen_cc_ues(self.n_cc)
            self.gen_d2d_pairs(self.n_pairs)

    def gen_cc_ues(self, n_cc):
        self.cc_ue = [UE(p, self.ue_tp)
                      for p in gen_uni_circ_pos(0, self.cell_r, n_cc)]

    def gen_d2d_pairs(self, n_pairs):
        d2d_tr_pos = gen_uni_circ_pos(0, self.cell_r, n_pairs)
        self.d2d_tr = [UE(p, self.ue_tp) for p in d2d_tr_pos]
        self.d2d_rc = [UE(gen_uni_circ_pos(p, self.d2d_r, 1)[0], self.ue_tp)
                       for p in d2d_tr_pos]

    def set_cc_tps(self, idx_cc, tps):
        for i in idx_cc:
            self.cc_ue[i]['tp'] = tps[i]

    def set_d2d_tps(self, idx_d2d, tps):
        for i in idx_d2d:
            self.d2d_tr[i]['tp'] = tps[i]

    # def set_bs_tps(self, idx_channel, tps):
    #     for i in idx_channel:
    #         self.bs.tp[i] = tps[i]

    def cal_reuse(self, idx_d2d):
        """
        Return:
        idx_cc
        bw
        """
        pass

    def cal_interference(self, idx_d2d, idx_cc, method=UL):
        s = self.d2d_tr[idx_d2d]
        d = self.d2d_rc[idx_d2d]
        cc_ues = self.cc_ue[idx_cc]
        if method == 0:
            cc_inter_tr = np.array([ue.pos for ue in cc_ues])
            cc_inter_tp = np.array([ue.tp for ue in cc_ues])
            n_d2d_inter = self.n_cc
            n_channel = 1
        else:
            cc_inter_tr = self.bs.pos
            cc_inter_tp = self.bs.tp
            n_d2d_inter = 1
            n_channel = self.n_cc
        # !FIXME:  only [cc->d2d]
        return cal_recv_power(cc_inter_tr, d.pos, 10**(cc_inter_tp/10.0), n_d2d_inter, n_channel,
                              cal_dist_2d, [],
                              cal_umi_nlos, [3.5, s.freq],
                              fading_func=gen_rayleigh, fading_args=[],
                              shadowing_func=gen_logNshadowing, shadowing_args=[4])

    def cal_throughput(self, idx_d2d, inter=0, method=UL):
        tr = self.d2d_tr[idx_d2d]
        rv = self.d2d_rc[idx_d2d]
        d2d_intra = self.cal_interference(idx_d2d, list(range(self.n_cc)), method)
        bw = self.total_bw / self.n_cc
        noise = cal_thermal_noise(bw, 290.0)
        sp = cal_recv_power(tr.pos, rv.pos, 10**(tr.tp/10.0), 1, 1,
                            cal_dist_2d, [],
                            cal_umi_nlos, [3.5, tr.freq],
                            fading_func=gen_rayleigh, fading_args=[],
                            shadowing_func=gen_logNshadowing, shadowing_args=[4])
        d2d_cap = sum(cal_shannon_cap(bw, sp, d2d_intra+inter, noise))
        return d2d_cap
