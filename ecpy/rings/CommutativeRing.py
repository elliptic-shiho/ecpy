from .Ring import Ring, RingElement


class CommutativeRing(Ring):
    def __init__(s, element_class):
        Ring.__init__(s, element_class)


class CommutativeRingElement(RingElement):
    def __init__(s, ring, x):
        RingElement.__init__(s, ring, x)

    def __radd__(s, lhs):
        return s + lhs

    def __rsub__(s, lhs):
        return -s + lhs

    def __rmul__(s, lhs):
        return s * lhs
