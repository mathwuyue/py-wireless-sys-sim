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
    def __init__(self, pos, tp=23, height=25, freq=2.4):
        super(BS, self).__init__(pos, tp, height, None)
        self.freq = freq
