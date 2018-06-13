
class Variant:
    def __init__(self, chr, pos):
        self.chr = chr
        self.pos = pos
        self.segments = {}

    def add_segment(self, segment_id, variant_info):
        """
        Adds segment to variant object that contains a similar variant
        :param segment_id: 
        :param variant_info: 
        :return: 
        """
        self.segments[segment_id[2]] = variant_info
