from .Ring import Ring, RingElement


class QuotientRing(Ring):
    def __init__(s, base_ring, modulo):
        assert isinstance(base_ring, Ring)
        Ring.__init__(s, QuotientRingElement)
        s.base_ring = base_ring
        s.modulo = modulo

    def __str__(s):
        return "Quotient ring of %s modulo %s" % (s.base_ring, s.modulo)

    def __repr__(s):
        return "%s(%r, %r)" % (s.__class__.__name__, s.base_ring, s.modulo)


class QuotientRingElement(RingElement):
    def __init__(s, ring, x):
        RingElement.__init__(s, ring, x)
        s.x = s.x % ring.modulo

    def lift(s):
        return s._lift(s)

    def _lift(s, t):
        if isinstance(t, QuotientRingElement):
            assert t.ring == s.ring
            x = t.x
        else:
            x = t
        return x

    def __add__(s, rhs):
        x = s.x
        y = s._lift(rhs)
        return s.__class__(s.ring, x + y)

    def __sub__(s, rhs):
        x = s.x
        y = s._lift(rhs)
        return s.__class__(s.ring, x - y)

    def __mul__(s, rhs):
        x = s.x
        y = s._lift(rhs)
        return s.__class__(s.ring, x * y)

    def __rmul__(s, lhs):
        x = s._lift(lhs)
        y = s.x
        return s.__class__(s.ring, x * y)

    def __rtruediv__(s, lhs):
        return s.__rdiv__(lhs)

    def __rfloordiv__(s, lhs):
        return s.__rdiv__(lhs)

    def __truediv__(s, rhs):
        return s.__div__(rhs)

    def __floordiv__(s, rhs):
        return s.__div__(rhs)

    def __div__(s, rhs):
        from ecpy.utils import egcd

        x = s.x
        y = s._lift(rhs)
        return s.__class__(s.ring, x * egcd(y, s.ring.modulo)[1])

    def __rdiv__(s, lhs):
        from ecpy.utils import egcd

        x = s._lift(lhs)
        y = s.x
        return s.__class__(s.ring, x * egcd(y, s.ring.modulo)[1])

    def __neg__(s):
        return s.__class__(s.ring, s.x)

    def __eq__(s, rhs):
        x = s.x
        y = s._lift(rhs)
        return x == y

    def __str__(s):
        return str(s.x)

    def __iter__(s):
        if hasattr(s.x, "__iter__"):
            return iter(s.x)
        return RingElement.__iter__(s)
