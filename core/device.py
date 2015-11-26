import numpy as np


class Device(object):
    def __init__(self, pos, tp, height, antenna):
        self.pos = pos
        self.tp = tp
        self.height = height
        self.antenna = antenna


class UE(Device):
    def __init__(self, pos, tp=-10, height=1.5, freq=2.4):
        super(UE, self).__init__(pos, tp, height, None)
        self.freq = freq


class BS(Device):
    def __init__(self, pos, n_channel=1, tp=16, height=25, antenna=0):
        super(BS, self).__init__(pos, tp, height, antenna)
        self.n_channel = n_channel
        if n_channel > 1 and (type(tp) is int or tp.shape == ()):
            self.tp = tp * np.ones(n_channel)
        self.antenna = antenna

    def cal_antenna_gain(self, device):
        phi = np.angle(device.pos)
        theta = np.arctan((self.height-device.height) / np.abs(device.pos))
        return self.antenna.cal_gain(phi, theta)
