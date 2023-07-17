from .CommutativeRing import CommutativeRing, CommutativeRingElement


class Integer(CommutativeRing):
    def __init__(s):
        CommutativeRing.__init__(s, IntegerElement)

    def _add(s, a, b):
        return a[0] + b[0]

    def _neg(s, a):
        return -a[0]

    def _equ(s, a, b):
        return len(a) == len(b) == 1 and a[0] == b[0]


class IntegerElement(CommutativeRingElement):
    def __init__(s, ring, x):
        CommutativeRingElement.__init__(s, ring, x)

    def __mul__(s, rhs, mod=None):
        import six

        if isinstance(rhs, IntegerElement):
            return s.__class__(s.ring, s.x * rhs.x)
        elif isinstance(rhs, six.integer_types):
            return s.__class__(s.ring, s.x * rhs)
        else:
            raise ValueError("Invalid Value: %r" % rhs)

    def __floordiv__(s, rhs):
        return s.__div__(rhs)

    def __truediv__(s, rhs):
        return s.__div__(rhs)

    def __rfloordiv__(s, lhs):
        return s.__rdiv__(lhs)

    def __rtruediv__(s, lhs):
        return s.__rdiv__(lhs)

    def __div__(s, rhs):
        import six

        if rhs == 0:
            raise ZeroDivisionError()
        if isinstance(rhs, IntegerElement):
            return s.__class__(s.ring, s.x // rhs.x)
        elif isinstance(rhs, six.integer_types):
            return s.__class__(s.ring, s.x // rhs)
        elif hasattr(rhs, "__rfloordiv__"):
            return rhs.__rfloordiv__(s.x)
        else:
            raise ValueError("Invalid Value: %r" % rhs)

    def __rdiv__(s, lhs):
        import six

        if lhs == 0:
            raise ZeroDivisionError()
        if isinstance(lhs, IntegerElement):
            return s.__class__(s.ring, lhs.x // s.x)
        elif isinstance(lhs, six.integer_types):
            return s.__class__(s.ring, lhs // s.x)
        elif hasattr(lhs, "__floordiv__"):
            return lhs.__floordiv__(s.x)
        else:
            raise ValueError("Invalid Value: %r" % rhs)


ZZ = Integer()
