from numpy import np


def gen_uni_circ_pos(center, r, n):
    theta = np.random.rand(n) * 2 * np.pi
    r1 = np.random.rand(n) * r
    r2 = np.random.rand(n) * r
    rpos = np.maximum(r1, r2)
    return rpos * np.exp(theta*1j)
