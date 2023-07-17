from ecpy.fields.FractionField import FractionField, FractionFieldElement
from ecpy.rings.Integer import ZZ
from ecpy.utils import gcd


class RationalField(FractionField):
    def __init__(s):
        FractionField.__init__(s, ZZ)
        s.element_class = RationalFieldElement


class RationalFieldElement(FractionFieldElement):
    def __init__(s, *args):
        FractionFieldElement.__init__(s, *args)
        p, q = s.x, s.y
        s.p = int(p // gcd(p, q))
        s.q = int(q // gcd(p, q))

    def __int__(s):
        if s.p % s.q != 0:
            raise ValueError("Can't divide a value: %s" % s)
        return s.p // s.q


QQ = RationalField()
