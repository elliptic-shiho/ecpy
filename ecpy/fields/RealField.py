from .Field import Field, FieldElement


class RealField(Field):
    def __init__(s):
        Field.__init__(s, RealFieldElement)

    def __repr__(s):
        return "RealField()"

    def __str__(s):
        return "RealField()"

    def _add(s, a, b):
        return s.element_class(s, a[0] + b[0])

    def _mul(s, a, b):
        return s.element_class(s, a[0] * b[0])

    def _neg(s, a):
        return s.element_class(s, -a[0])

    def _inv(s, a):
        return s.element_class(s, 1.0 / a[0])

    def _equ(s, a, b):
        d = a[0] - b[0]
        if isinstance(d, RealFieldElement):
            d = d.x
        return abs(d) < 0.00001


class RealFieldElement(FieldElement):
    pass


RR = RealField()
