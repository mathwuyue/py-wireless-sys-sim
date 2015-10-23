from core.position import cal_dist_2d
from core.communication import cal_recv_power, cal_shannon_cap, cal_thermal_noise
from core.pathloss import cal_umi_nlos
from core.statchannel import gen_rayleigh, gen_logNshadowing
import numpy as np


class LTEUSystem(object):
    def __init__(self, device=None, total_bw=20e6):
        self.device = device
        self.total_bw = total_bw
        self.channel_bw = 10e6
        self.n = self.total_bw / self.channel_bw

    def lbt(self):
        pass


class LTEUSystemSimple(LTEUSystem):
    def __init__(self, total_bw=20e6, prob=0.5):
        super(LTEUSystemSimple, self).__init__(None, total_bw)
        self.prob = prob

    def lbt(self):
        lbt_seed = np.random.rand(self.n)
        return self.channel_bw * sum(lbt_seed <= self.prob)

    def set_prob(self, prob):
        self.prob = prob

    def cal_lteu_throughput(self, tr, rv, interference=0):
        bw = self.lbt()
        if bw == 0:
            return 0
        noise = cal_thermal_noise(bw, 290.0)
        sp = cal_recv_power(tr.pos, rv.pos, 10**(tr.tp/10.0), 1, 1,
                            cal_dist_2d, [],
                            cal_umi_nlos, [3.5, tr.freq],
                            fading_func=gen_rayleigh, fading_args=[],
                            shadowing_func=gen_logNshadowing, shadowing_args=[4])
        return cal_shannon_cap(bw, sp, interference, noise)
