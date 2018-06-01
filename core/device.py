
class Device(object):
    def __init__(self, pos, tp, height, antenna, freq):
        self.pos = pos
        self.tp = tp
        self.height = height
        self.antenna = antenna
        self.freq = freq

    def __getitem__(self, key):
        return self.__dict___.get(key)

    def __setitem__(self, key, value):
        self.__dict__[key] = value


class UE(Device):
    def __init__(self, pos, tp=-10, freq=2.4):
        super(UE, self).__init__(pos, tp, 1.5, None, freq)


class BS(Device):
    def __init__(self, pos, tp=23, height=25, antenna=None, freq=2.4):
        super(BS, self).__init__(pos, tp, height, antenna, freq)
