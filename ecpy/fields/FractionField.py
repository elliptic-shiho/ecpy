from ecpy.utils import gcd, lcm
from .Field import Field, FieldElement


class FractionField(Field):
    def __init__(s, ring):
        from ecpy.rings.Ring import Ring

        assert issubclass(ring.__class__, Ring)
        super(FractionField, s).__init__(FractionFieldElement)
        s.ring = ring

    def _add(s, a, b):
        if a[1] == b[1]:
            return s.element_class(s, a[0] + b[0], a[1])
        else:
            l = a[1] * b[1]
            p1 = a[0] * l / a[1]
            p2 = b[0] * l / b[1]
            return s.element_class(s, p1 + p2, l)

    def _mul(s, a, b):
        return s.element_class(s, a[0] * b[0], a[1] * b[1])

    def _neg(s, a):
        return s.element_class(s, -a[0], a[1])

    def _inv(s, a):
        return s.element_class(s, a[1], a[0])

    def _equ(s, a, b):
        return a[0] == b[0] and a[1] == b[1]

    def _mod(s, a, b):
        return s.element_class(s, a[0] % b[0], a[1])

    def __str__(s):
        return "Fraction Field of %s" % s.ring

    def __repr__(s):
        return "%s(%r)" % (s.__class__.__name__, s.ring)


class FractionFieldElement(FieldElement):
    def __init__(s, *args):
        p, q = None, 1
        if len(args) >= 1:
            if isinstance(args[0], s.__class__):
                v = args[0]
                p = v.p
                q = v.q
                field = v.field
            elif len(args) >= 2:
                field = args[0]
                p = args[1]
                q = 1
                if len(args) == 3:
                    q = args[2]
        if p == None or q == None:
            raise TypeError("Invalid Argument: '%s'" % (tuple(map(str, args)),))
        super(FractionFieldElement, s).__init__(field, p)
        s.y = q
        if hasattr(s.ring.ring, "_div"):
            s.p = p / gcd(p, q)
            s.q = q / gcd(p, q)
        else:
            s.p = p
            s.q = q

    def __repr__(s):
        return "%r(%r, %r)" % (s.field, s.p, s.q)

    def __str__(s):
        if s.q == 1:
            return "%s" % s.p
        return "%s/%s" % (s.p, s.q)

    def __iter__(s):
        return (s.p, s.q).__iter__()

    def _to_tuple(s, d):
        if isinstance(d, s.__class__):
            return tuple(d)
        elif isinstance(d, tuple):
            return d
        else:
            return (d, 1)
