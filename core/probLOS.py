import numpy as np


def prob_exp_los(d, l):
    """ Calculate the probability of LOS and NLOS.
    Returns:
        int. 0 for NLOS, 1 for LOS
    """
    rnd = np.random.rand(1)[0]
    if rnd > np.exp(-(d/l)**2):
        return 0
    else:
        return 1
