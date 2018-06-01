from d2d.model import D2DSystemModel


class Rui2016(D2DSystemModel):
    """Benchmark from Rui2016 paper"""
    def __init__(self, n_cc, n_pairs, ue_tp=-10, total_bw = 10e6, cell_r=100, d2d_r=20):
        super(Rui2016, self).__init__(n_cc, n_pairs, ue_tp, total_bw, cell_r, d2d_r)

    def gen_cc_ues(self, n_cc=None):
        if not n_cc:
            n_cc = self.n_cc
        super(Rui2016, self).gen_cc_ues(n_cc)

    def gen_d2d_pairs(self, n_pairs=None):
        if not n_pairs:
            n_pairs = self.n_pairs
        super(Rui2016, self).gen_d2d_pairs(n_pairs)

    def cal_sinr(self):
