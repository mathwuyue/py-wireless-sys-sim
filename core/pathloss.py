import numpy as np


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
