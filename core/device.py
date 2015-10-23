import numpy as np


class UE(object):
    def __init__(self, pos, tp=-10, height=1.5, freq=2.4):
        self.pos = pos
        self.tp = tp
        self.height = height
        self.freq = freq


class BS(object):
    def __init__(self, pos, n_channel=1, tp=16, height=25, freq=2.4, antenna=0):
        self.pos = pos
        self.n_channel = n_channel
        if n_channel > 1 and (type(tp) is int or tp.shape == ()):
            self.tp = tp * np.ones(n_channel)
        else:
            self.tp = tp
        self.height = height
        self.freq = freq
        self.antenna = antenna
