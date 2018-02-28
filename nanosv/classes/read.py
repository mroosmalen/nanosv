#!/usr/bin/python

class Read:
    def __init__(self, qname, length):
        self.qname = qname
        self.length = length
        self.segments = dict()

    def addSegment(self, segment):
        """
        Adds segment to read
        :param segment:
        """
        self.segments[segment.clip] = segment.id
