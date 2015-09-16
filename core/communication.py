import numpy as np


def cal_channel_gain(tr, rv, n_channel=1,
                     dist_func=None, dist_args=None,
                     pl_func=None, pl_args=None,
                     fading_func=None, fading_args=None,
                     shadowing_func=None, shadowing_args=None):
    if tr.shape != rv.shape:
        return None
    if len(tr.shape) == 1:
        n = 1
    else:
        n = tr.shape[0]
    d = dist_func(tr, rv)
    pl = pl_func(d, pl_args)
    h = fading_func(fading_args+[1, n_channel*n])
    s = shadowing_func(shadowing_args+[1, n])
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
