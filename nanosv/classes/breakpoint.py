#!/usr/bin/python

class Breakpoint:
    def __init__(self, id, segment_1, segment_2):
        self.id = id
        if segment_1.flag & 16:
            flag_1 = 16
        else:
            flag_1 = 0
        if segment_2.flag & 16:
            flag_2 = 16
        else:
            flag_2 = 0
        self.segment_1 = {"id": segment_1.id, "rname": segment_1.rname, "flag": flag_1, "mapq":segment_1.mapq}
        self.segment_2 = {"id": segment_2.id, "rname": segment_2.rname, "flag": flag_2, "mapq":segment_2.mapq}

    def setGap(self, gap):
        """
        Saves gap in class
        :param gap:
        """
        self.gap = gap

    def setBreakpoint(self, segment_1, segment_2):
        """
        Saves start positions of segments depending on the flag
        :param segment_1 used to access and save data:
        :param segment_2 used to access and save data:
        """
        self.segment_1["pos"] = segment_1.pos if segment_1.flag & 16 else segment_1.end
        self.segment_2["pos"] = segment_2.end if segment_2.flag & 16 else segment_2.pos

    def switchSegments(self):
        """
        Switches segments in case segment 2 comes before segment 1 in the genome
        """
        self.segment_1, self.segment_2 = self.segment_2, self.segment_1
        self.segment_1["flag"] = 0 if self.segment_1["flag"] & 16 else 16
        self.segment_2["flag"] = 0 if self.segment_2["flag"] & 16 else 16

    def setSVtype(self):
        """
        Sets SV-type to INS in case the difference between the two segments is smaller than the gap setting
        """
        self.svtype = "BND"
        if abs(self.segment_2["pos"] - self.segment_1["pos"]) < self.gap:
            self.svtype = "INS"
            self.segment_2["pos"] = self.segment_1["pos"]+1
            self.segment_1["flag"] = 0
            self.segment_2["flag"] = 0
