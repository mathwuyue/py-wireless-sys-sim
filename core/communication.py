import numpy as np


def cal_channel_gain(tr, rv, n, n_channel=1,
                     dist_func=None, dist_args=[],
                     pl_func=None, pl_args=[],
                     fading_func=None, fading_args=[],
                     shadowing_func=None, shadowing_args=[]):
    if tr.shape != rv.shape:
        return None
    d = dist_func(tr, rv, *dist_args)
    pl = pl_func(d, *pl_args)
    fading_args = fading_args + [n*n_channel] if n*n_channel > 1 else fading_args
    h = fading_func(*fading_args)
    shadowing_args = shadowing_args + [n] if n > 1 else shadowing_args
    s = shadowing_func(*shadowing_args)
    return np.kron(10**((-pl + s)/10), np.ones((1, n_channel))) * (abs(h)**2)


def cal_recv_power(tr, rv, tp, n_channel,
                   dist_func, dist_args,
                   pl_func, pl_args,
                   fading_func, fading_args,
                   shadowing_func, shadowing_args):
    return np.kron(tp, np.ones((1, n_channel))) *\
        cal_channel_gain(tr, rv, n_channel,
                         dist_func, dist_args,
                         pl_func, pl_args,
                         fading_func, fading_args,
                         shadowing_func, shadowing_args)
