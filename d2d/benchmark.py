import numpy as np
from d2d.model import D2DSystemModel
from core.communication import cal_recv_power, cal_SINR, cal_thermal_noise
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
        self.n_cc_rb = self.n_rb / self.n_cc
        self.n_d2d_rb = self.n_rb

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
        cc_inter_tr = np.reshape(cc_inter_tr, (self.n_pairs, 1))
        cc_inter_rv = np.zeros(1, self.n_cc)
        cc_inter_tp = np.reshape(d2d_tps, (self.n_pairs, 1))
        for i in range(self.n_cc):
            c_ip = cal_recv_power(cc_inter_tr, cc_inter_rv, cc_inter_tp, self.n_cc_rb,
                                  cal_dist_2d, [],
                                  cal_umi_nlos, [alpha, self.d2d_tr[0].freq],
                                  fading_func=gen_rayleigh, fading_args=[1.0],
                                  shadowing_func=gen_logNshadowing, shadowing_args=[4])
        cc_interference = np.sum(c_ip, axis=0)
        # cal d2d interference
        # cc -> d2d
        cc_d2d_ip = []
        d2d_rvs = np.array([d.pos for d in self.d2d_rc])
        for d2d_rv in d2d_rvs:
            cc_tr = np.array([c.pos for c in self.cc_ue])
            cc_tr = np.reshape(cc_tr, (self.n_cc, 1))
            cc_inter_tps = np.reshape(np.array(cc_tps), (self.n_cc, 1))
            tmp = cal_recv_power(d2d_cc_tr, d2d_rv, cc_inter_tps, self.n_cc_rb,
                                 cal_dist_2d, [],
                                 cal_umi_nlos, [alpha, self.cc_ue[0].freq],
                                 fading_func=gen_rayleigh, fading_args=[1.0],
                                 shadowing_func=gen_logNshadowing, shadowing_args=[4])
            cc_d2d_ip.append(np.reshape(tmp, (1, self.n_rb)))
        # other d2d -> d2d
        d2d_d2d_ip = []
        for idx, d in enumerate(d2d_rvs):
            d2d_inter_tr = d2d_rvs.copy()
            del d2d_inter_tr[idx]
            d2d_inter_tr = np.reshape(d2d_inter_tr, (self.n_pairs-1, 1))
            d2d_inter_tps = d2d_tps.copy()
            del d2d_inter_tps[idx]
            d2d_inter_tps = np.reshape(d2d_inter_tps, (self.n_pairs-1, 1))
            d_d2d_ip = cal_recv_power(d2d_inter_tr, d, d2d_inter_tps, self.n_rb,
                                      cal_dist_2d, [],
                                      cal_umi_nlos, [alpha, self.cc_ue[0].freq],
                                      fading_func=gen_rayleigh, fading_args=[1.0],
                                      shadowing_func=gen_logNshadowing, shadowing_args=[4])
            d_d2d_ip = np.sum(d_d2d_ip, axis=0)
            d2d_d2d_ip.append(d_d2d_ip)
        d2d_interference = np.array(cc_d2d_ip) + np.array(d2d_d2d_ip)
        return cc_interference, d2d_interference

    def cal_signal_power(self, cc_tps, d2d_tps):
        # cc signal recv power
        cc_trs = np.reshape(np.array([c.pos for c in self.cc_ue]), (self.n_cc, 1))
        cc_tps = np.reshape(cc_tps, (self.n_cc, 1))
        cc_rvs = np.zeros(1, self.n_cc)
        cc_signal_power = cal_recv_power(cc_trs, cc_rvs, cc_tps, self.n_cc_rb,
                                         cal_dist_2d, [],
                                         cal_umi_nlos, [alpha, self.cc_ue[0].freq],
                                         fading_func=gen_rayleigh, fading_args=[1.0],
                                         shadowing_func=gen_logNshadowing, shadowing_args=[4])
        # d2d signal recv power
        d2d_trs = [d.pos for d in self.d2d_tr]
        d2d_rvs = [d.pos for d in self.d2d_rc]
        d2d_signal_power = []
        for d2d_tr, d2d_rv, d2d_tp in zip(d2d_trs, d2d_rvs, d2d_tps):
            d_signal = cal_recv_power(d2d_tr, d2d_rv, d2d_tp, self.n_cc_rb,
                                      cal_dist_2d, [],
                                      cal_umi_nlos, [alpha, self.cc_ue[0].freq],
                                      fading_func=gen_rayleigh, fading_args=[1.0],
                                      shadowing_func=gen_logNshadowing, shadowing_args=[4])
            d2d_signal_power.append(d_signal)
        return cc_signal_power, d2d_signal_power

    def cal_ak_bk(self, d2d_s, d2d_interference):
        rb_bw = self.total_bw / self.n_rb
        noise = cal_thermal_noise(rb_bw, 290)
        akn = []
        bkn = []
        # d2d sinr
        for i in range(self.n_pairs):
            d_interference = d2d_interference[i, :]
            d_s = d2d_s[i, :]
            d2d_sinr = cal_SINR(d_s, d_interference, noise)
            ak = d2d_sinr / (1+d2d_sinr)
            bk = np.log(1+d2d_sinr) - (d2d_sinr/(1+d2d_sinr))
            akn.append(ak)
            bkn.append(bk)
        return akn, bkn

    def run(self, **kwargs):
        # update simulation parameters
        for key in kwargs:
            self.__dict__[key] = kwargs[key]
        # start simulation
        self.gen_cc_ues(self.n_cc)
        self.gen_d2d_pairs(self.n_pairs)
