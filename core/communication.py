import numpy as np
import numpy.matlib


def cal_channel_gain(tr, rv, n_row, n_col, n_channel,
                     dist_func=None, dist_args=[],
                     pl_func=None, pl_args=[],
                     fading_func=None, fading_args=[],
                     shadowing_func=None, shadowing_args=[]):
    d = dist_func(tr, rv, *dist_args)
    pl = pl_func(d, *pl_args)
    fading_args = fading_args + [n_row, n_col]
    h = fading_func(*fading_args)
    shadowing_args = shadowing_args + [n_row, n_col/n_channel]
    s = shadowing_func(*shadowing_args)
    return np.kron(10**((-pl+s)/10.0), np.ones(n_channel)) * (abs(h)**2)


def cal_recv_power(tr, rv, tp, n_channel,
                   dist_func, dist_args,
                   pl_func, pl_args,
                   fading_func, fading_args,
                   shadowing_func, shadowing_args):
    n_row = len(tr) if hasattr(tr, '__len__') else 1
    n_col = len(rv) * n_channel if hasattr(tr, '__len__') else n_channel
    if n_channel != 1:
        if type(tp) is np.ndarray:
            tp = np.matlib.repmat(tp, 1, n_col)
        else:
            tp = tp * np.ones(n_row, n_col)
    tr = np.matlib.repmat(tr, 1, len(rv))
    rv = np.matlib.repmat(rv, n_row, 1)
    return tp * cal_channel_gain(tr, rv, n_row, n_col, n_channel,
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


def cal_shannon_cap(bw, sp, ip, noise):
    sinr = cal_SINR(sp, ip, noise)
    return bw * np.log2(1+sinr)


def cal_transmission_time(throughput, packet=1e3, period=1e-3):
    """
    Args:
    throughput (float): bps
    packet (float): bits
    period (float): time slot (s)

    Return:
    Total transmission time (float): in s
    """
    pass
