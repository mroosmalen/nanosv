#!/usr/bin/python


class Segment:
    def __init__(self, qname, flag, rname, pos, mapq, length, clip, clip_2, pid, rlength):
        self.id = [rname, pos, qname + ";" + str(clip)]
        self.qname = qname
        self.rlength = rlength
        self.flag = int(flag)
        self.rname = rname
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.length = int(length)
        self.end = (int(pos) + int(length))
        self.clip = clip
        self.clip_2 = clip_2
        self.pid = pid

    def setPlength(self, rlength):
        """
        Calculates plength
        :param rlength is used to calculate the plength:
        """
        self.plength = format(self.length / rlength, '.3f')
