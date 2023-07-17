class Ring(object):
    def __init__(s, element_class):
        s.element_class = element_class

    def __repr__(s):
        return "%s()" % s.__class__.__name__

    def __str__(s):
        return s.__class__.__name__

    def __call__(s, *x):
        return s.element_class(s, *x)

    def _add(s, a, b):
        raise NotImplementedError()

    def _mul(s, a, b):
        raise NotImplementedError()

    def _neg(s, a, b):
        raise NotImplementedError()

    def _equ(s, a, b):
        raise NotImplementedError()

    def order(s):
        return 0

    def _ord(s, a):
        return 0

    def degree(s):
        return 1


class RingElement(object):
    def __init__(s, ring, x):
        s.ring = ring
        s.x = x

    def change_ring(s, _ring):
        return s.__class__(_ring, *tuple(s))

    def __mul__(s, rhs, mod=None):
        return s.ring._mul(tuple(s), s._to_tuple(rhs))

    def __rmul__(s, lhs):
        return s.ring._mul(s._to_tuple(lhs), tuple(s))

    def __getitem__(s, idx):
        return s._to_tuple(s)[idx]

    def __len__(s):
        return 1

    def order(s):
        return s.ring._ord(tuple(s))

    def int(s):
        return int(s.x)

    def __add__(s, rhs):
        return s.ring._add(tuple(s), s._to_tuple(rhs))

    def __sub__(s, rhs):
        return s.ring._add(tuple(s), s._to_tuple(-rhs))

    def __neg__(s):
        return s.ring._neg(tuple(s))

    def __radd__(s, lhs):
        return s.ring._add(s._to_tuple(lhs), tuple(s))

    def __rsub__(s, lhs):
        return s.ring._add(s._to_tuple(-lhs), tuple(s))

    def __ne__(s, rhs):
        return not (s == rhs)

    def __eq__(s, rhs):
        return s.ring._equ(tuple(s), s._to_tuple(rhs))

    def __repr__(s):
        return "%r(%s)" % (s.ring, s.x)

    def __str__(s):
        return "%s" % s.x

    def _to_tuple(s, d):
        if isinstance(d, s.__class__):
            return tuple(d)
        elif isinstance(d, tuple):
            return d
        else:
            return (d,)

    def __iter__(s):
        return (s.x,).__iter__()

    def __int__(s):
        return s.int()

    def __hash__(s):
        return hash(str(s.x) + str(s.ring))
