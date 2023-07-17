from ecpy.fields import FiniteField, ExtendedFiniteField
from ecpy.utils.util import is_enable_native, _native
from ecpy.utils import modular_square_root
from random import randint


def EllipticCurve(field, *args, **kwargs):
    """
    Return Elliptic Curve Instance.
    """
    if isinstance(field, FiniteField):
        return FiniteFieldEllipticCurve(field, *args, **kwargs)
    else:
        return GenericEllipticCurve(field, *args, **kwargs)


def EC_from_j_invariant(field, j0):
    """
    Return Elliptic Curve Instance which has j-invariant is `j0`.
    """
    assert isinstance(field, FiniteField)
    R = randint(1, field.n)
    if j0 == 0:
        return FiniteFieldEllipticCurve(field, 0, R)
    if j0 == 1728:
        return FiniteFieldEllipticCurve(field, R, 0)
    ADR = 3 * j0 * R**2
    BDR = 2 * j0 * R**3
    u = 1728 - j0
    uinv = 1 / field(u)
    return FiniteFieldEllipticCurve(field, ADR * uinv, BDR * uinv)


class GenericEllipticCurve(object):
    """
    Elliptic Curve on General Field
    """

    def __init__(s, field, a, b):
        """
        Constructor of Elliptic Curve.
          y^2 = x^3 + `a`x+ `b` on `field`
        """
        s.element_class = GenericEllipticCurvePoint
        s.field = field
        s.a = a
        s.b = b
        s.O = s.element_class(s, 0, 1, 0)

    def is_on_curve(s, point):
        """
        Is on curve `point`?
        """
        return s._is_on_curve(point.x, point.y, point.z)

    def _is_on_curve(s, x, y, z=1):
        """
        Is on curve (`x`, `y`, `z`)?
        """
        x = s.field(x)
        y = s.field(y)
        z = s.field(z)
        return y * y * z == x * x * x + s.a * x * z * z + s.b * z * z * z

    def determinant(s):
        """
        Calculate Determinant of Curve.
        """
        return -16 * (4 * s.a**3 + 27 * s.b**2)

    def j_invariant(s):
        """
        Calculate j-Invariant of Curve.
        """
        return -1728 * ((4 * s.a**3) / s.determinant())

    def __repr__(s):
        return "EllipticCurve(%r, %r, %r)" % (s.field, s.a, s.b)

    def __call__(s, *x):
        return s.element_class(s, *x)

    def __str__(s):
        res = "Elliptic Curve y^2 = x^3"
        if s.a != 0:
            if s.a == 1:
                res += " + x"
            else:
                res += " + %rx" % s.a
        if s.b != 0:
            res += " + %r" % s.b
        res += " over %r" % s.field
        return res

    def _add(s, P, Q):
        from six import integer_types

        """
    Add Operation on Perspective Coordinate
    P : tuple (x, y, z)
    Q : tuple (u, v, w)
    return: R = P + Q
    """
        Px, Py, Pz = P
        Qx, Qy, Qz = Q
        Rx, Ry, Rz = s.O
        if s._equ(P, Q):
            X, Y, Z = Px, Py, Pz
            u = 3 * X * X + s.a * Z * Z
            v = Y * Z
            a = Y * v
            w = u * u - 8 * X * a
            Rx = 2 * v * w
            Ry = u * (4 * X * a - w) - 8 * a * a
            Rz = 8 * v * v * v
        else:
            u = Qy * Pz - Py * Qz
            v = Qx * Pz - Px * Qz
            v2 = v * v
            v3 = v2 * v
            w = u * u * Pz * Qz - v3 - 2 * v2 * Px * Qz
            Rx = v * w
            Ry = u * (v2 * Px * Qz - w) - v3 * Py * Qz
            Rz = v3 * Pz * Qz
        if isinstance(Rz, integer_types):
            z = 1 // s.field(Rz)
        else:
            z = 1 // Rz
        if z == 0:
            return s.O
        return s.element_class(s, Rx * z, Ry * z, 1)

    def _equ(s, P, Q):
        """
        Is P equals to Q?
        """
        return P[0] * Q[1] == P[1] * Q[0]

    def _neg(s, P):
        """
        return -P
        """
        if P == (0, 1, 0):
            return s.O
        return s.element_class(s, P[0], -P[1])


class GenericEllipticCurvePoint(object):
    """
    Elliptic Curve Point on General Field
    """

    def __init__(s, group, x, y, z=1):
        def F(x):
            if isinstance(x, tuple):
                return group.field(*x)
            return group.field(x)

        s.group = group
        s.x = F(x)
        s.y = F(y)
        s.z = F(z)
        s.inf = s.x == 0 and s.y == 1 and s.z == 0
        if not (s.inf or s.group.is_on_curve(s)):
            raise ArithmeticError("Invalid Point: (%s, %s, %s)" % (s.x, s.y, s.z))

    def is_infinity(s):
        """
        Returns:
          Is self equals to O?
        """
        return s.inf

    def change_group(s, _group):
        return s.__class__(_group, *tuple(s))

    def line_coeff(s, Q):
        from six.moves import map

        """
    Calculate Line Coefficient of Line self to Q
    """
        P = s
        x1, y1, z1 = map(s.group.field, P)
        x2, y2, z2 = map(s.group.field, Q)
        assert z1 == z2 == 1  # is normalized?
        if x1 == x2:
            l = (3 * x1 * x1 + s.group.a) // (2 * y1)
        else:
            l = (y2 - y1) // (x2 - x1)
        return l

    def __add__(s, rhs):
        if isinstance(rhs, GenericEllipticCurvePoint) and rhs.is_infinity():
            return s
        d = s._to_tuple(rhs)
        if s.is_infinity():
            return s.__class__(s.group, d[0], d[1])
        else:
            return s.group._add(tuple(s), d)

    def __sub__(s, rhs):
        return s + (-rhs)

    def __mul__(s, rhs):
        """
        Multiplication Operation
        """
        from six.moves import map

        d = rhs
        if d < 0:
            b = -1
            d = -d
        else:
            b = 1
        if d == 0:
            return s.group.O
        bits = list(map(int, bin(d)[2:]))[::-1]
        x = s
        if bits[0]:
            res = x
        else:
            res = s.group.O
        for cur in bits[1:]:
            x += x
            if cur:
                res += x
        if b == -1:
            res = -res
        return res

    def __neg__(s):
        """
        Negate Operation Wrapper
        """
        return s.group._neg(tuple(s))

    def __rmul__(s, lhs):
        return s * lhs

    def __eq__(s, rhs):
        """
        Is self equals to rhs?
        """
        if rhs == None:
            return False
        return s.group._equ(tuple(s), s._to_tuple(rhs))

    def _to_tuple(s, d):
        if isinstance(d, GenericEllipticCurvePoint):
            return tuple(d)
        elif isinstance(d, tuple):
            return d
        else:
            raise ArithmeticError("Invalid Parameter: %r" % d)

    def __iter__(s):
        return (s.x, s.y, s.z).__iter__()

    def __repr__(s):
        if s.is_infinity():
            return "%r.O" % s.group
        return "%r(%r, %r, %r)" % (s.group, s.x, s.y, s.z)

    def __str__(s):
        if s.is_infinity():
            return "Infinity Point (0 : 1 : 0) on %s" % s.group
        return "Point (%s : %s : %s) on %s" % (s.x, s.y, s.z, s.group)


class FiniteFieldEllipticCurve(GenericEllipticCurve):
    """
    Elliptic Curve on Finite Field or Extended Finite Field
    """

    def __init__(s, field, a, b):
        s.field = field
        s.a = a
        s.b = b
        s.element_class = FiniteFieldEllipticCurvePoint
        s.O = s.element_class(s, 0, 1, 0)
        if is_enable_native:
            if isinstance(field, ExtendedFiniteField):
                cond = {1: "x^2+1", 2: "x^2+x+1"}
                poly = cond[field.t]
                s.base = _native.EF(field.p, poly)
            elif isinstance(field, FiniteField):
                s.base = _native.FF(field.p)
            s.ec = _native.EC(s.base, a, b)
            s._add = s.__add_native
        else:
            s.add_func = s.__add

    def get_corresponding_y(s, x):
        """
        Calculate `y` coordinate corresponding to given x.
        """
        x = s.field(x)
        y_square = x * x * x + s.a * x + s.b
        for y in modular_square_root(y_square, s.field.p ** s.field.degree()):
            if pow(y, 2, s.field.p ** s.field.degree()) == y_square:
                return y
        return None

    def embedding_degree(s, m):
        """
        Calculate Embedding Degree.
          <=> minimum `k` satisfy m | p^k - 1
        """
        k = 1
        while True:
            if (s.field.p ** (k * s.field.degree()) - 1) % m == 0:
                return k
            k += 1

    def random_point(s):
        """
        return random point on this curve.
        """
        rnd = [randint(0, s.field.order())] * s.field.degree()
        x = s.field(*rnd)
        while True:
            y = s.get_corresponding_y(x)
            if y != None:
                if s._is_on_curve(x, y):
                    return s.element_class(s, x, y)
            x += 1

    def __add(s, P, Q):
        return super(FiniteFieldEllipticCurve, s)._add(P, Q)

    def __add_native(s, P, Q):
        R = _native.EC_elem(s.ec, 0, 0)
        P = _native.EC_elem(s.ec, tuple(P[0]), tuple(P[1]), tuple(P[2]))
        Q = _native.EC_elem(s.ec, tuple(Q[0]), tuple(Q[1]), tuple(Q[2]))
        s.ec.add(R, P, Q)
        R = FiniteFieldEllipticCurvePoint(s, *R.to_python(), normalize=True)
        return R


class FiniteFieldEllipticCurvePoint(GenericEllipticCurvePoint):
    def __init__(s, group, x, y, z=1, normalize=False):
        def F(x):
            if type(x) == tuple:
                return group.field(*x)
            return group.field(x)

        s.group = group
        s.x = F(x)
        s.y = F(y)
        s.z = F(z)
        if normalize and s.z != 0:
            s.x = s.x / s.z
            s.y = s.y / s.z
            s.z = s.z / s.z
        s.inf = s.x == 0 and s.y == 1 and s.z == 0
        if not (s.inf or s.group.is_on_curve(s)):
            raise ArithmeticError("Invalid Point: (%s, %s, %s)" % (s.x, s.y, s.z))
        if is_enable_native:
            s.__mul_method__ = s._mul_native
        else:
            s.__mul_method__ = super(FiniteFieldEllipticCurvePoint, s).__mul__

    def _mul_native(s, rhs):
        P = tuple(s)
        R = _native.EC_elem(s.group.ec, 0, 1, 0)
        P = _native.EC_elem(s.group.ec, tuple(P[0]), tuple(P[1]), tuple(P[2]))
        m = rhs
        s.group.ec.mul(R, P, m)
        R = FiniteFieldEllipticCurvePoint(s.group, *R.to_python(), normalize=True)
        return R

    def __mul__(s, rhs):
        return s.__mul_method__(rhs)

    def distortion_map(s):
        """
        IMPORTANT: If you want to use this function,
                    definition field should be Extended Finite Field.
        return \phi(self), \phi is Distortion map
        Polynomial: x^2+1 or x^2+x+1
        """

        def to_tuple(x):
            from six import integer_types

            if type(x) in integer_types:
                return (x, 0)
            return tuple(x)

        x = to_tuple(s.x)
        y = to_tuple(s.y)
        if s.group.field.t == 1:
            x = (-x[0], -x[1])
            y = (y[1], y[0])
        elif s.group.field.t == 2:
            x = (x[1], x[0])
            y = (y[0], y[1])
        return s.__class__(s.group, x, y)

    def order(s):
        """
        return order of self
        """
        r = s.change_group(s.group)
        i = 2
        while True:
            if r.is_infinity():
                return i
            r += s
            i += 1
