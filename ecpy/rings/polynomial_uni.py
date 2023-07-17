from copy import deepcopy
from ecpy.utils import memoize
from .Ring import Ring, RingElement
from six.moves import map, xrange, zip_longest


class UnivariatePolynomialRing(Ring):
    """
    Univariate Polynomial Ring
    """

    def __init__(s, field, gen="x"):
        """
        Args:
          field : A base field class
        """
        Ring.__init__(s, UnivariatePolynomialElement)
        s.field = field
        s.gen_name = gen

    def gen(s):
        return s.element_class(s, 0, 1)  # x

    def __str__(s):
        return "Univariate Polynomial Ring over %s" % s.field

    def _add(s, A, B):
        if len(A) == 1:
            return s.element_class(s, [A[0] + B[0]] + list(B[1:]))
        elif len(B) == 1:
            return s.element_class(s, [A[0] + B[0]] + list(A[1:]))
        ret = []
        for x, y in zip_longest(A, B, fillvalue=0):
            ret += [(x + y)]
        return s.element_class(s, ret)

    def _mul(s, A, B):
        if len(A) == 1:
            return s.element_class(s, map(lambda x: A[0] * x, B))
        elif len(B) == 1:
            return s.element_class(s, map(lambda x: x * B[0], A))
        ret = [0] * (len(A) + len(B) - 1)
        for x, y in enumerate(A):
            for u, v in enumerate(B):
                ret[x + u] += y * v
        return s.element_class(s, ret)

    def _neg(s, A):
        return s.element_class(s, map(lambda x: -x, A))

    def _equ(s, A, B):
        from six.moves import zip_longest

        if len(A) == 1 and len(B) == 1:
            return A[0] == B[0]
        return all([x == y for x, y in zip_longest(A, B, fillvalue=0)])


class UnivariatePolynomialElement(RingElement):
    def __init__(s, poly_ring, *args):
        """
        coeffs = [x^0, x^1, x^2, ... degree's coefficient]
        """
        RingElement.__init__(s, poly_ring, args)
        if len(args) == 1:
            v = args[0]
            if isinstance(v, UnivariatePolynomialElement):
                s.coeffs = v.coeffs
            elif hasattr(v, "__iter__"):
                s.coeffs = list(v)
            else:
                s.coeffs = [v]
        else:
            s.coeffs = args
        s.trim()
        s.coeffs = list(map(s.ring.field, s.coeffs))

    def __rdivmod__(s, lhs):
        if not isinstance(lhs, UnivariatePolynomialElement):
            lhs = s.__class__(s.ring, lhs)
        return divmod(lhs, s)

    def __divmod__(s, rhs):
        assert rhs != 0
        if isinstance(rhs, UnivariatePolynomialElement):
            if len(rhs) == 1:
                q = UnivariatePolynomialElement(
                    s.ring, map(lambda x: x / rhs[0], s.coeffs)
                )
                r = UnivariatePolynomialElement(
                    s.ring, map(lambda x: x % rhs[0], s.coeffs)
                )
                return q, r
            q = UnivariatePolynomialElement(s.ring, 0)
            r = s
            d = rhs.degree()
            c = rhs[-1]
            while r.degree() >= d:
                t = UnivariatePolynomialElement(s.ring, [r[-1] / c]).shift(
                    r.degree() - d
                )
                q = q + t
                r = r - t * rhs
            return q, r
        else:
            q = UnivariatePolynomialElement(s.ring, map(lambda x: x / rhs, s.coeffs))
            r = UnivariatePolynomialElement(s.ring, map(lambda x: x % rhs, s.coeffs))
            return q, r

    def __div__(s, rhs):
        return divmod(s, rhs)[0]

    def __truediv__(s, rhs):
        return s.__div__(rhs)

    def __floordiv__(s, rhs):
        return s.__div__(rhs)

    def __getitem__(s, idx):
        return s.coeffs[idx]

    def __iter__(s):
        return iter(s.coeffs)

    def __len__(s):
        return len(s.coeffs)

    def __mod__(s, rhs):
        return divmod(s, rhs)[1]

    def __rmod__(s, lhs):
        return divmod(lhs, s)[1]

    def __pow__(s, rhs, mod=0):
        if rhs == 0:
            return UnivariatePolynomialElement(s.ring, 1)
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
            res = UnivariatePolynomialElement(s.ring, 1)
        b = 0
        for cur in bits[1:]:
            b += 1
            x *= x
            if mod > 0:
                x %= mod
                x.trim()
            if cur:
                res *= x
                if mod > 0:
                    res %= mod
                    res.trim()
        return res

    def __str__(s):
        if len(s.coeffs) == 0:
            return "0"
        res = []
        deg = 0
        for x in s:
            r = ""
            if x != 0:
                if (x != -1 and x != 1) or deg == 0:
                    r += str(x)
                elif x == -1:
                    r += "-"
                if deg > 0:
                    r += s.ring.gen_name
                    if deg > 1:
                        r += "^%d" % deg
            if r != "":
                res += [r]
            deg += 1
        if len(res) == 0:
            return "0"
        return "+".join(res[::-1]).replace("+-", "-")

    @memoize
    def apply(s, x0):
        """
        Using Horner's method
        """
        ret = 0
        for x in s.coeffs[::-1]:
            ret = ret * x0 + x
        return ret

    def degree(s):
        return len(s.coeffs) - 1

    def monic(s):
        return s / s[-1]

    def shift(s, i):
        if i > 0:
            return UnivariatePolynomialElement(s.ring, [0] * abs(i) + s.coeffs)
        else:
            return UnivariatePolynomialElement(s.ring, s.coeffs[abs(i) :])

    def trim(s):
        i = 0
        for x in s.coeffs[::-1]:
            if x != 0:
                break
            i += 1
        if 0 < i < len(s.coeffs):
            s.coeffs = s.coeffs[:-i]
        elif i == len(s.coeffs):
            s.coeffs = [0]

    def __hash__(s):
        return hash("".join(map(str, s.coeffs)) + str(s.ring))
