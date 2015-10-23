import numpy as np
from d2d.model import D2DSystemModel
from core.position import gen_uni_circ_pos, gen_uni_ring_pos
from core.device import UE
from lte.lteU import LTEUSystemSimple


class D2DDLModel(D2DSystemModel):
    def __init__(self, n_cc, n_pairs):
        super(D2DDLModel, self).__init__(n_cc, n_pairs)

    def gen_d2d_pairs(self, n_pairs):
        d2d_tr_pos = gen_uni_circ_pos(0, 100, n_pairs)
        self.d2d_tr = np.array([UE(p) for p in d2d_tr_pos])
        self.d2d_rc = [UE(gen_uni_circ_pos(p, 20, 1)[0]) for p in d2d_tr_pos]


class D2DULModel(D2DSystemModel):
    def __init__(self, n_cc, n_pairs):
        super(D2DULModel, self).__init__(n_cc, n_pairs)

    def gen_d2d_pairs(self, n_pairs):
        d2d_tr_pos = gen_uni_ring_pos(0, 150, 250, n_pairs)
        self.d2d_tr = np.array([UE(p) for p in d2d_tr_pos])
        self.d2d_rc = [UE(gen_uni_circ_pos(p, 20, 1)[0]) for p in d2d_tr_pos]


def main():
    lteu_system = LTEUSystemSimple()
    dl_system = D2DDLModel(5, 1)
    ul_system = D2DULModel(5, 1)
    loop = 10000

    # lte-u
    lteu_ue1 = UE(0)
    lteu_ue2 = [UE(p) for p in gen_uni_circ_pos(0, 20, loop)]
    # free
    throughput = 0
    lteu_system.set_prob(1.8)
    for i in range(loop):
        throughput = throughput + lteu_system.cal_lteu_throughput(lteu_ue1, lteu_ue2[i], 1e-8)
    print "LTE-U Free: ", throughput / loop

    # busy
    throughput = 0
    lteu_system.set_prob(0.2)
    for i in range(loop):
        throughput = throughput + lteu_system.cal_lteu_throughput(lteu_ue1, lteu_ue2[i], 1e-6)
    print "LTE-U Busy: ", throughput / loop


    # conventional
    # busy, UL
    throughput = 0
    for i in range(loop):
        ul_system.gen_d2d_pairs(1)
        ul_system.gen_cc_ues(5)
        throughput = throughput + ul_system.cal_throughput(0, 1e-6)
    print 'UL, busy: ', throughput/loop
    
    # busy, DL
    throughput = 0
    for i in range(loop):
        dl_system.gen_d2d_pairs(1)
        dl_system.gen_cc_ues(5)
        throughput = throughput + dl_system.cal_throughput(0, 1e-6, method=1)
    print 'DL, busy: ', throughput/loop
    
    # free, UL
    throughput = 0
    ul_system.set_cc_tps([0, 1, 2, 3], -80)
    for i in range(loop):
        ul_system.gen_cc_ues(5)
        ul_system.gen_d2d_pairs(1)
        throughput = throughput + ul_system.cal_throughput(0, 1e-8)
    print 'UL, free: ', throughput/loop

    # free, DL
    throughput = 0
    dl_system.set_bs_tps([0, 1, 2, 3], -80)
    for i in range(loop):
        dl_system.gen_cc_ues(5)
        dl_system.gen_d2d_pairs(1)
        throughput = throughput + dl_system.cal_throughput(0, 1e-8, method=1)
    print 'DL, free: ', throughput/loop


if __name__ == '__main__':
    main()
