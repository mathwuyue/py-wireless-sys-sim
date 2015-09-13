import numpy as np


def gen_rayleigh(sigma=1, *args):
    return np.sqrt(sigma/2) * np.random.randn(*args) +\
        np.sqrt(sigma/2) * np.random.randn(*args) * 1j


def gen_rician(k=10, sigma=1, *args):
    return np.sqrt(k/k+1)*sigma*np.exp(1j) +\
        np.sqrt((1/k+1)*sigma/2) * np.random.randn(*args) +\
        np.sqrt((1/k+1)*sigma/2) * np.random.randn(*args) * 1j
