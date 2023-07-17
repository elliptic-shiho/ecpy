from ..rings.CommutativeRing import CommutativeRing, CommutativeRingElement


class Field(CommutativeRing):
    def __init__(s, element_class):
        CommutativeRing.__init__(s, element_class)

    def _inv(s, a):
        raise NotImplementedError()

    def _mod(s, a, b):
        raise NotImplementedError()

    def _div(s, a, b):
        return s._mul(a, (s._inv(b)))


class FieldElement(CommutativeRingElement):
    def __init__(s, field, x):
        CommutativeRingElement.__init__(s, field, x)
        s.field = field
        if isinstance(x, FieldElement):
            x = x.x
        s.x = x

    def change_field(s, _field):
        return s.__class__(_field, *tuple(s))

    def __mul__(s, rhs):
        return s.field._mul(tuple(s), s._to_tuple(rhs))

    def __div__(s, rhs):
        return s.field._div(tuple(s), s._to_tuple(rhs))

    def __rdiv__(s, lhs):
        return s.field._div(s._to_tuple(lhs), tuple(s))

    def __truediv__(s, rhs):
        return s.field._div(tuple(s), s._to_tuple(rhs))

    def __rtruediv__(s, lhs):
        return s.field._div(s._to_tuple(lhs), tuple(s))

    def __floordiv__(s, rhs):
        if hasattr(s.field, "_fdiv"):
            return s.field._fdiv(tuple(s), s._to_tuple(rhs))
        else:
            return s.field._div(tuple(s), s._to_tuple(rhs))

    def __rfloordiv__(s, lhs):
        if hasattr(s.field, "_fdiv"):
            return s.field._fdiv(s._to_tuple(lhs), tuple(s))
        else:
            return s.field._div(s._to_tuple(lhs), tuple(s))

    def __pow__(s, rhs, mod=None):
        from six.moves import map

        if rhs == 0:
            return s.__class__(s.field, 1)
        d = int(rhs)
        if d < 0:
            x = 1 / s
            d = -d
        else:
            x = s
        bits = list(map(int, bin(d)[2:]))[::-1]
        if bits[0]:
            res = x
        else:
            res = s.field(1)
        for cur in bits[1:]:
            x *= x
            if cur:
                res *= x
        return res

    def __mod__(s, rhs):
        return s.field._mod(tuple(s), s._to_tuple(rhs))

    def __rmod__(s, lhs):
        return s.field._mod(s._to_tuple(lhs), tuple(s))

    def __rmul__(s, lhs):
        return s * lhs

    def int(s):
        return int(s.x)

    def __ne__(s, rhs):
        return not (s == rhs)

    def __eq__(s, rhs):
        return s.field._equ(tuple(s), s._to_tuple(rhs))

    def __iter__(s):
        return (s.x,).__iter__()

    def __int__(s):
        return s.int()

    def __hash__(s):
        return s.x
