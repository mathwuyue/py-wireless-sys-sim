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
    tmp = 10**((-pl + s)/10) if n_channel == 1 else np.kron(10**((-pl + s)/10), np.ones(n_channel))
    return tmp * (abs(h)**2)


def cal_recv_power(tr, rv, tp, n, n_channel,
                   dist_func, dist_args,
                   pl_func, pl_args,
                   fading_func, fading_args,
                   shadowing_func, shadowing_args):
    tp = tp if n_channel == 1 else np.kron(tp, np.ones(n_channel))
    return tp * cal_channel_gain(tr, rv, n, n_channel,
                                 dist_func, dist_args,
                                 pl_func, pl_args,
                                 fading_func, fading_args,
                                 shadowing_func, shadowing_args)


def cal_thermal_noise(bw, t):
    K = 1.3806488e-23
    return bw*t*K


def cal_SINR(sp, ip, noise):
    """
    Args:
    sp (float or numpy array): signal power
    ip (float or numpy array): interference power
    """
    return sp / (ip+noise)


def cal_shannon_cap(bw, sp, ip=0, noise=0):
    sinr = cal_SINR(sp, ip, noise)
    return bw * np.log2(1+sinr)
