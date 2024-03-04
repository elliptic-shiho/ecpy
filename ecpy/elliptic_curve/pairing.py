def miller(E, P, Q, m):
    from six.moves import map

    """
  Calculate Divisor by Miller's Algorithm
  Args:
    E: The Elliptic Curve
    P: A point over E which has order m
    Q: A point over E which has order m to apply function f_P
    m: The order of P, Q on E
  Returns:
    f_P(Q)
  """

    def h(P, Q, R):
        # if \lambda is infinity
        if (P == Q and P.y == 0) or (P != Q and P.x == Q.x):
            return R.x - P.x
        L = P.line_coeff(Q)
        p = R.y - P.y - L * (R.x - P.x)
        q = R.x + P.x + Q.x - L * L
        return p / q

    if P == Q:
        return 1
    b = map(int, bin(m)[2:])
    next(b)
    f = 1
    T = P
    for i in b:
        f = f * f * h(T, T, Q)
        T = T + T
        if i:
            f = f * h(T, P, Q)
            T = T + P
    return f


def weil_pairing(E, P, Q, m, S=None):
    """
    Calculate Weil Pairing
    Args:
      E: The Elliptic Curve
      P: A point over E which has order m
      Q: A point over E which has order m
      m: The order of P, Q on E
      S: [Optional] A random point on E
    Returns:
      e_m(P, Q)
    """
    if S is None:
        S = E.random_point()
    from ecpy.utils.util import is_enable_native, _native
    from ecpy.fields.ExtendedFiniteField import ExtendedFiniteFieldElement

    if is_enable_native:
        P = _native.EC_elem(E.ec, tuple(P.x), tuple(P.y), tuple(P.z))
        Q = _native.EC_elem(E.ec, tuple(Q.x), tuple(Q.y), tuple(Q.z))
        S = _native.EC_elem(E.ec, tuple(S.x), tuple(S.y), tuple(S.z))
        if E.ec.type == 1:
            t = _native.FF_elem(0)
        elif E.ec.type == 2:
            t = _native.EF_elem(0, 0)
        _native.weil_pairing(t, E.ec, P, Q, S, m)
        if E.ec.type == 1:
            return t.to_python()
        elif E.ec.type == 2:
            t = t.to_python()
            return ExtendedFiniteFieldElement(E.field, t[0], t[1])
    else:
        fpqs = miller(E, P, Q + S, m)
        fps = miller(E, P, S, m)
        fqps = miller(E, Q, P - S, m)
        fqs = miller(E, Q, -S, m)
        return E.field._inv(fps * fqps) * fpqs * fqs


def tate_pairing(E, P, Q, m, k=2):
    """
    Calculate Tate Pairing
    Args:
      E: The Elliptic Curve
      P: A point over E which has order m
      Q: A point over E which has order m
      m: The order of P, Q on E
      k: [Optional] The Embedding Degree of m on E
    """
    from ecpy.utils.util import is_enable_native, _native

    if is_enable_native:
        P = _native.EC_elem(E.ec, tuple(P.x), tuple(P.y), tuple(P.z))
        Q = _native.EC_elem(E.ec, tuple(Q.x), tuple(Q.y), tuple(Q.z))
        if E.ec.type == 1:
            t = _native.FF_elem(0)
        elif E.ec.type == 2:
            t = _native.EF_elem(0, 0)
        _native.tate_pairing(t, E.ec, P, Q, m, k)
        if E.ec.type == 1:
            from ecpy.fields.Zmod import ZmodElement

            return ZmodElement(E.field, t.to_python())
        elif E.ec.type == 2:
            from ecpy.fields.ExtendedFiniteField import ExtendedFiniteFieldElement

            t = t.to_python()
            return ExtendedFiniteFieldElement(E.field, t[0], t[1])
    else:
        f = miller(E, P, Q, m)
        return f ** (((E.field.p**k) - 1) // m)


def MapToPoint(E, y):
    """
    MapToPoint Function: Given by Boneh-Durfee's ID-based Encryption Paper.
    Args:
      E: The Elliptic Curve
      y: Any Value (should be E.field element)

    Returns:
      Correspond point of y on E
    """
    from ecpy.utils import cubic_root

    x = cubic_root(y**2 - 1)
    Q = E(x, y)
    return 6 * Q


def gen_supersingular_ec(bits=70):
    """
    Generate Super-Singluar Elliptic Curve
    Args:
      bits: The Security Parameter: log_2 p = bits

    Returns:
      A (Super Singular) Elliptic Curve, Extended Finite Field, l
      l is need to calculate Pairing
    """
    from ecpy.fields import ExtendedFiniteField
    from .EllipticCurve import EllipticCurve

    def _next_prime(n):
        from ecpy import is_prime

        """
    return next prime of n
    """
        while not is_prime(n):
            n += 1
        return n

    """
  If you have gmpy, use gmpy.next_prime
  in other hand, use slow function
  """
    try:
        from gmpy import next_prime
    except:
        next_prime = _next_prime

    def gen_prime():
        from ecpy import is_prime
        from random import randint

        while True:
            p = int(next_prime(randint(2 ** (bits - 1), 2**bits)))
            if is_prime(p * 6 - 1):
                break
        return p * 6 - 1, p

    p, l = gen_prime()
    F = ExtendedFiniteField(p, "x^2+x+1")
    return EllipticCurve(F, 0, 1), F, l


def find_point_by_order(E, l):
    """
    Find a Elliptic Curve Point P which has order l.
    Args:
      E: The Elliptic Curve
      l: Order of Point on E

    Returns:
      Point on E which has order l.
    """
    i = 3
    while True:
        r = E.get_corresponding_y(i)
        if r != None:
            P = E(i, r)
            if (P * l).is_infinity():
                return P
        i += 1


def symmetric_weil_pairing(E, P, Q, m):
    """
    Symmetric Weil Pairing
    \hat{e}(P, Q) = e(P, \phi(Q)) (\phi is Distortion Map)
    Args:
      E: The Elliptic Curve
      P: A point on E which has order m
      Q: A point on E which has order m
      m: The order of P, Q
    """
    return weil_pairing(E, P, Q.distortion_map(), m)


def symmetric_tate_pairing(E, P, Q, m, k=2):
    """
    Symmetric Tate Pairing
    \hat{e}(P, Q) = e(P, \phi(Q)) (\phi is Distortion Map)
    Args:
      E: The Elliptic Curve
      P: A point on E which has order m
      Q: A point on E which has order m
      m: The order of P, Q
      k: [Optional] The Embedding Degree of m on E
    """
    return tate_pairing(E, P, Q.distortion_map(), m)
