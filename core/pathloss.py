import numpy as np
from core.probLOS import prob_exp_los


def cal_fiirs(d, f, gt, gr):
    """ d (m)
    f (GHz)
    """
    return 20*np.log10(d) + 20*np.log10(f) + 32.45 - gt - gr


def cal_umi_nlos(d, alpha, f):
    """
    d (m)
    f (GHz)
    """
    return 10*alpha*np.log10(d) + 22.7 + 26*np.log10(f)


def cal_umi_exp_los(d, l, alpha_los, alpha_nlos, k_los, k_nlos, f):
    los = prob_exp_los(d, l)
    if los:
        return 10*alpha_los*np.log10(d) + k_los + 26*np.log10(f)
    else:
        return 10*alpha_nlos*np.log10(d) + k_nlos + 26*np.log10(f)
