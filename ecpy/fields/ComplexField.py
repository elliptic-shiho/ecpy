from .Field import Field, FieldElement


class ComplexField(Field):
    def __init__(s):
        super(ComplexField, s).__init__(ComplexFieldElement)

    def _add(s, a, b):
        return s.element_class(s, a[0] + b[0], a[1] + b[1])

    def _neg(s, a):
        return s.element_class(s, -a[0], -a[1])

    def _mul(s, a, b):
        return s.element_class(s, a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

    def _equ(s, a, b):
        return a[0] == b[0] and a[1] == b[1]

    def _inv(s, a):
        a = s.element_class(s, a[0], a[1])
        x = a.absolute_value()
        abar = a.conjugate()
        return s.element_class(s, abar.x / x, abar.y / x)

    def _div(s, z, w):
        a, b, c, d = z[0], z[1], w[0], w[1]
        u = 1.0 / (c**2 + d**2)
        return s.element_class(s, (a * c + b * d) * u, (b * c - a * d) * u)


class ComplexFieldElement(FieldElement):
    def __init__(s, field, x, y=0):
        super(ComplexFieldElement, s).__init__(field, x)
        s.y = y

    def norm(s):
        """
        Calculate Norm
        Returns:
          ||self||
        """
        return s.x**2 - s.y**2

    def absolute_value(s):
        """
        Calculate Absolute Value
        Returns:
          |self|
        """
        import math

        return math.sqrt(s.x**2 + s.y**2)

    def conjugate(s):
        """
        Calculate conjugate complex
        Returns:
          \overline{self}
        """
        return s.__class__(s.field, s.x, -s.y)

    def __repr__(s):
        return "%r(%r, %r)" % (s.field, s.x, s.y)

    def __str__(s):
        return "%r + %ri" % (s.x, s.y)

    def __iter__(s):
        return (s.x, s.y).__iter__()

    def _to_tuple(s, d):
        if isinstance(d, s.__class__):
            return tuple(d)
        elif isinstance(d, tuple):
            return d
        else:
            return (d, 0)


CC = ComplexField()
