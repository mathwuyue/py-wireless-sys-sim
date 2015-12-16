import pickle

import numpy as np
from d2d import cal_D2D_ergodic_tp, D2DSystemModel, cal_ergodic_subopt_tp

#import ipdb; ipdb.set_trace()

d2d_radius = np.linspace(10, 120, 12)
t_period = np.linspace(10, 1e5, 12)

L = 82.5
ALPHAL = 2.09
ALPHANL = 3.75
KL = 10**10.38
KNL = 10**14.54

RC = 1.8
RC_RATIO = 0.8

def main():
    # # period
#    d2dsystem = D2DSystemModel(1, 1)     #
#    d2d_tr = d2dsystem.d2d_tr[0].pos
#    d2d_rc = d2dsystem.d2d_rc[0].pos
#    cc_ue = d2dsystem.cc_ue[0].pos
    with open('test.pkl') as f:
        d2d_tr = pickle.load(f)
        cc_ue = pickle.load(f)
        d2d_rc = pickle.load(f)
    beta, delta = cal_D2D_ergodic_tp(d2d_tr, d2d_rc, cc_ue, RC, 1, 1,
                                     KL, KNL, ALPHAL, ALPHANL, L)
    print beta, delta

    # sub-opt
    delta = cal_ergodic_subopt_tp(d2d_tr, d2d_rc, cc_ue, RC_RATIO, 1, 1, 0.8,
                          KL, KNL, ALPHAL, ALPHANL, L)
    print delta

    with open('test.pkl', 'wb') as f:
        pickle.dump(d2d_tr, f)
        pickle.dump(cc_ue, f)
        pickle.dump(d2d_rc, f)


if __name__ == '__main__':
    main()
