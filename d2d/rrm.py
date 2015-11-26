import numpy as np
from core import cal_thermal_noise, cal_umi_nlos, cal_umi_exp_los


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


def cal_D2D_ergodic_tp()
