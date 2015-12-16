import operator
import itertools

import numpy as np
import scipy.optimize
from core import cal_thermal_noise, cal_umi_nlos, cal_umi_exp_los


def _sum(func, *args):
    return reduce(operator.add, itertools.imap(func, *args), 0)


def cal_D2D_basic_tp(d2d_ues, g_d2d_bs, kappa, bw, alpha, freq):
    """
    This function calculates the transmit power for D2D UEs (Spectrum Sharing Scheme Between Cellular Users and Ad-hoc Device-to-Device Users)
    Args:
        d2d_ues (numpy array): d2d_ues positions
        g_d2d_cc (): channel gain between d2d and cc ues
        kappa (float): scale param for cc
        bw (float): bandwidth for d2d_ues
        alpha (float): pathloss parameter
        freq (float): frequency

    Returns:
        numpy array. The transmit power of D2D UEs.
    """
    noise = cal_thermal_noise(bw, 273)
    pathloss = cal_umi_nlos(np.abs(d2d_ues), alpha, freq)
    return (kappa - 1) * pathloss * noise / g_d2d_bs


def cal_D2D_opt_tp(d2d_ues, cc_ues,
                   pmax_d, pmax_c,
                   g_d2d_bs, g_cc, g_d2d, g_cc_d2d,
                   sinr_d2d, sinr_cc,
                   bw, alpha, freq):
    """
    This function calculates the RRM for D2D UEs (Device-to-Device Communications
Underlaying Cellular Networks)
    Args:
        d2d_ues (numpy array): d2d_ues positions
        g_d2d_cc (): channel gain between d2d and cc ues
        kappa (float): scale param for cc
        bw (float): bandwidth for d2d_ues
        alpha (float): pathloss parameter
        freq (float): frequency

    Returns:
        list of numpy array. The transmit power of D2D UEs and CC UEs.

    ::TODO: only consider one D2D
    """
    noise = cal_thermal_noise(bw, 273)
    # set up reuse array
    idx_avail = []
    p_c = (g_d2d*sinr_cc+g_d2d_bs*sinr_cc*sinr_d2d)*noise / \
          (g_d2d*g_cc-sinr_d2d*sinr_cc*g_cc_d2d*g_d2d_bs)
    p_d2d = (g_cc_d2d*sinr_cc*sinr_d2d+g_cc*sinr_d2d)*noise / \
            (g_d2d*g_cc-sinr_cc*sinr_d2d*g_cc_d2d*g_d2d_bs)
    for i in range(cc_ues.size):
        if (p_d2d > 0 and p_d2d <= pmax_c) and (p_c > 0 and p_c <= pmax_c):
            idx_avail.append(i)

    # calculate optimal transmit power
    # FIXME: one D2D
    def _argmax(tp_pairs):
        f = 0
        idx = 0
        for i, (pc, pd) in enumerate(tp_pairs):
            fc = np.log2(1+pc*g_cc/(pd*g_d2d_bs+noise))+np.log2(1+pd*g_d2d/(pc*g_cc_d2d+noise))
            if fc > f:
                f = fc
                idx = i
        return tp_pairs[idx]

    p1 = (pmax_c*g_cc_d2d[idx_avail]+noise)*sinr_d2d/g_d2d
    p2 = (pmax_c*g_cc[idx_avail]-sinr_cc*noise)/(sinr_cc*g_d2d_bs)
    p3 = (pmax_d*g_d2d-sinr_d2d*noise)/(sinr_d2d*g_cc_d2d[idx_avail])
    p4 = (pmax_d*g_d2d_bs+noise)*sinr_cc/g_cc[idx_avail]
    opt_tp_pairs = []
    for i, j in enumerate(idx_avail):
        if (pmax_c*g_cc[i])/(noise+pmax_d*g_d2d_bs) <= sinr_cc:
            opt_tp_pairs.append(_argmax([(pmax_c, p1[j]), (pmax_c, p2[j])]))
        elif pmax_d*g_d2d/(noise+pmax_c*g_cc_d2d[i]) < sinr_d2d:
            opt_tp_pairs.append(_argmax([(p3[j], pmax_d), (p4[j], pmax_d)]))
        else:
            opt_tp_pairs.append(_argmax([(pmax_c, p1[j]), (pmax_c, pmax_d), (p4[j], pmax_d)]))

    # calculate channel allocation.
    return _argmax(opt_tp_pairs)


def cal_D2D_ergodic_tp(d2d_tr, d2d_rc, cc_ue, rc, a_gain_c, a_gain_d,
                       k_los, k_nlos, alpha_los, alpha_nlos, l):
    def _f(x):
        return x*np.log2(x)/(np.log(2)*(x-1))
    a_c = a_gain_c/a_gain_d       # antenna gain from CC to D2D
    a_d = 1                       # antenna gain from D2D to BS
    d1_d = np.abs(d2d_tr - d2d_rc)
    d2_c = np.abs(d2d_tr)
    d1_c = np.abs(cc_ue)
    d2_d = np.abs(cc_ue - d2d_rc)

    # M, N
    def _m1(a, d1, d2):
        return a*(d1/d2)**(-alpha_los)

    def _m2(a, d1, d2):
        return a*k_los*d1**(-alpha_los) / (k_nlos*d2**(-alpha_nlos))

    def _m3(a, d1, d2):
        return a*k_nlos*d1**(-alpha_nlos) / (k_los*d2**(-alpha_los))

    def _m4(a, d1, d2):
        return a*(d1/d2)**(-alpha_nlos)

    def _n1(d1, d2):
        return np.exp(-(d1**2+d2**2)/l**2)

    def _n2(d1, d2):
        return np.exp(-d1**2/l**2) * (1-np.exp(-d2**2/l**2))

    def _n3(d1, d2):
        return np.exp(-d2**2/l**2) * (1-np.exp(-d1**2/l**2))

    def _n4(d1, d2):
        return (1-np.exp(-d1**2/l**2)) * (1-np.exp(-d2**2/l**2))

    # equation
    def _f_beta_delta(beta):
        delta = (rc - _sum(lambda x, y: y*_f(x/beta),
                           [_m1(a_c, d1_c, d2_c), _m2(a_c, d1_c, d2_c),
                            _m3(a_c, d1_c, d2_c), _m4(a_c, d1_c, d2_c)],
                           [_n1(d1_c, d2_c), _n2(d1_c, d2_c), _n3(d1_c, d2_c), _n4(d1_c, d2_c)])) / \
                _sum(lambda x, y: y*_f(beta*x)-y*_f(x/beta),
                     [_m1(a_c, d1_c, d2_c), _m2(a_c, d1_c, d2_c),
                      _m3(a_c, d1_c, d2_c), _m4(a_c, d1_c, d2_c)],
                     [_n1(d1_c, d2_c), _n2(d1_c, d2_c), _n3(d1_c, d2_c), _n4(d1_c, d2_c)])
        # lambda1 = _sum(lambda x, y: y*_f(beta*x)-y*_f(x/beta),
        #                [_m1(a_d, d1_d, d2_d), _m2(a_d, d1_d, d2_d),
        #                 _m3(a_d, d1_d, d2_d), _m4(a_d, d1_d, d2_d)],
        #                [_n1(d1_d, d2_d), _n2(d1_d, d2_d), _n3(d1_d, d2_d), _n4(d1_d, d2_d)]) / \
        #         _sum(lambda x, y: y*_f(x/beta)+y*_f(beta*x),
        #              [_m1(a_c, d1_c, d2_c), _m2(a_c, d1_c, d2_c),
        #               _m3(a_c, d1_c, d2_c), _m4(a_c, d1_c, d2_c)],
        #              [_n1(d1_c, d2_c), _n2(d1_c, d2_c), _n3(d1_c, d2_c), _n4(d1_c, d2_c)])
        lambda1 = 0

        # h1, h2
        def _h1(x, y):
            return x*y*(1-x/beta)/np.log(2)+np.log2(x/beta)/(np.log(2)*(x-beta)**2)

        def _h2(x, y):
            return y*_f(beta*x)/beta + beta*x*y*((1-1.0/(beta*x))/np.log(2)-x*np.log2(beta*x)) / \
                (np.log(2)*(beta*x-1)**2)

        return _sum(lambda xc, yc, xd, yd: delta*_h1(xd, yd)+(1-delta)*_h2(xd, yd)-lambda1*((1-delta)*_h1(xc, yc)+delta*_h2(xc, yc)),
                     [_m1(a_c, d1_c, d2_c), _m2(a_c, d1_c, d2_c),
                      _m3(a_c, d1_c, d2_c), _m4(a_c, d1_c, d2_c)],
                    [_n1(d1_c, d2_c), _n2(d1_c, d2_c), _n3(d1_c, d2_c), _n4(d1_c, d2_c)],
                    [_m1(a_d, d1_d, d2_d), _m2(a_d, d1_d, d2_d),
                     _m3(a_d, d1_d, d2_d), _m4(a_d, d1_d, d2_d)],
                    [_n1(d1_d, d2_d), _n2(d1_d, d2_d), _n3(d1_d, d2_d), _n4(d1_d, d2_d)])

    opt_beta = scipy.optimize.brentq(_f_beta_delta, 1e-5, 1-1e-5)
    opt_delta = (rc - _sum(lambda x, y: y*_f(x/opt_beta),
                       [_m1(a_c, d1_c, d2_c), _m2(a_c, d1_c, d2_c),
                        _m3(a_c, d1_c, d2_c), _m4(a_c, d1_c, d2_c)],
                       [_n1(d1_c, d2_c), _n2(d1_c, d2_c), _n3(d1_c, d2_c), _n4(d1_c, d2_c)])) / \
            _sum(lambda x, y: y*_f(opt_beta*x)-y*_f(x/opt_beta),
                 [_m1(a_c, d1_c, d2_c), _m2(a_c, d1_c, d2_c),
                  _m3(a_c, d1_c, d2_c), _m4(a_c, d1_c, d2_c)],
                 [_n1(d1_c, d2_c), _n2(d1_c, d2_c), _n3(d1_c, d2_c), _n4(d1_c, d2_c)])
    return opt_beta, opt_delta


def cal_ergodic_subopt_tp(d2d_tr, d2d_rc, cc_ue, rc_ratio, a_gain_c, a_gain_d, beta,
                          k_los, k_nlos, alpha_los, alpha_nlos, l):
    def _f(x):
        return x*np.log2(x)/(np.log(2)*(x-1))
    a_c = a_gain_c/a_gain_d       # antenna gain from CC to D2D
    d2_c = np.abs(d2d_tr)
    d1_c = np.abs(cc_ue)

    # M, N
    def _m1(a, d1, d2):
        return a*(d1/d2)**(-alpha_los)

    def _m2(a, d1, d2):
        return a*k_los*d1**(-alpha_los) / (k_nlos*d2**(-alpha_nlos))

    def _m3(a, d1, d2):
        return a*k_nlos*d1**(-alpha_nlos) / (k_los*d2**(-alpha_los))

    def _m4(a, d1, d2):
        return a*(d1/d2)**(-alpha_nlos)

    def _n1(d1, d2):
        return np.exp(-(d1**2+d2**2)/l**2)

    def _n2(d1, d2):
        return np.exp(-d1**2/l**2) * (1-np.exp(-d2**2/l**2))

    def _n3(d1, d2):
        return np.exp(-d2**2/l**2) * (1-np.exp(-d1**2/l**2))

    def _n4(d1, d2):
        return (1-np.exp(-d1**2/l**2)) * (1-np.exp(-d2**2/l**2))

    # rc
    max_rc = _sum(lambda x, y: y*_f(x/beta),
                  [_m1(a_c, d1_c, d2_c), _m2(a_c, d1_c, d2_c),
                   _m3(a_c, d1_c, d2_c), _m4(a_c, d1_c, d2_c)],
                  [_n1(d1_c, d2_c), _n2(d1_c, d2_c), _n3(d1_c, d2_c), _n4(d1_c, d2_c)])
    min_rc = _sum(lambda x, y: y*_f(x*beta),
                  [_m1(a_c, d1_c, d2_c), _m2(a_c, d1_c, d2_c),
                   _m3(a_c, d1_c, d2_c), _m4(a_c, d1_c, d2_c)],
                  [_n1(d1_c, d2_c), _n2(d1_c, d2_c), _n3(d1_c, d2_c), _n4(d1_c, d2_c)])
    rc = rc_ratio * (max_rc - min_rc) + min_rc
    print rc
    delta = (rc - max_rc) / (min_rc - max_rc)
    return delta
