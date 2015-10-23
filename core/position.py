import numpy as np
from scipy.spatial.distance import euclidean


def gen_uni_circ_pos(center, r, n):
    theta = np.random.rand(n) * 2 * np.pi
    r1 = np.random.rand(n) * r
    r2 = np.random.rand(n) * r
    rpos = np.maximum(r1, r2)
    return rpos * np.exp(theta*1j) + center


def gen_uni_ring_pos(center, r1, r2, n):
    theta = np.random.rand(n) * 2 * np.pi
    r = np.random.rand(n) * (r2-r1) + r1
    return r * np.exp(theta*1j) + center


def cal_dist_2d(p1, p2):
    """
    Args:
    p1 (numpy array or scalar): positions represented by complex number.
    p2 (numpy array or scalar): positions represented by complex number.
    .. note::
    If both p1 and p2 are numpy arrays, they should be same size.

    Returns:
    numpy array of distances.
    """
    return np.abs(p1 - p2)


def cal_dist_3d(p1, p2):
    """
    Args:
    p1 (numpy array or scalar): x, y, z.
    p2 (numpy array or scalar): x, y, z.
    .. note::
    If both p1 and p2 are numpy arrays, they should be same size.

    Returns:
    numpy array of distances.
    """
    # check scalar
    if len(p1.shape) == 1 and len(p2.shape) == 1:
        return euclidean(p1, p2)
    elif len(p1.shape) == 1:
        return np.array(map(lambda x: euclidean(x, p1), p2))
    elif len(p2.shape) == 1:
        return np.array(map(lambda x: euclidean(x, p2), p1))
    else:
        return np.array(map(euclidean, p1, p2))


def to_cartesian(ps):
    if len(ps.shape) == 1:
        return np.array([ps[0]*np.sin(ps[1])*np.sin(ps[2]),
                         ps[0]*np.sin(ps[1])*np.cos(ps[2]),
                         ps[0]*np.cos(ps[1])])
    else:
        return np.array([[p[0]*np.sin(p[1])*np.sin(p[2]),
                          p[0]*np.sin(p[1])*np.cos(p[2]),
                          p[0]*np.cos(p[1])] for p in ps])
