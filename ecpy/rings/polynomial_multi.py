from .polynomial_uni import UnivariatePolynomialElement
from .polynomial_uni import UnivariatePolynomialRing
from .Ring import Ring, RingElement
from six.moves import map, xrange, zip_longest
from ecpy.utils import memoize
from copy import deepcopy
import itertools


class BivariatePolynomialRing(Ring):
    """
    Bivariate Polynomial Ring
    """

    def __init__(s, field, gens):
        """
        Args:
          field : A base field class
          gens : generator names

        Example:
          PR = BivariatePolynomialRing(RR, ['x', 'y'])
        """
        Ring.__init__(s, BivariatePolynomialElement)
        assert len(gens) == 2
        s.field = field
        s.gen_names = gens

    def gens(s):
        ret = []
        ret += [s.element_class(s, [[0, 1]])]
        ret += [s.element_class(s, [[0], [1]])]
        return ret

    def _add(s, A, B):
        # A + B
        if len(A) == 1 and len(A[0]) == 1:
            A = A[0][0]
            return s.element_class(s, [[A + B[0][0]] + list(B[0][1:])] + list(B[1:]))
        elif len(B) == 1 and len(B[0]) == 1:
            B = B[0][0]
            return s.element_class(s, [[A[0][0] + B] + list(A[0][1:])] + list(A[1:]))
        ret = []
        for x, y in zip_longest(A, B, fillvalue=[0]):
            t = []
            for xs, ys in zip_longest(x, y, fillvalue=0):
                t += [xs + ys]
            ret += [t]
        return s.element_class(s, ret)

    def _mul(s, A, B):
        if len(A) == 1 and len(A[0]) == 1:
            A = A[0][0]
            return s.element_class(s, map(lambda y: map(lambda x: A * x, y), B))
        elif len(B) == 1 and len(B[0]) == 1:
            B = B[0][0]
            return s.element_class(s, map(lambda y: map(lambda x: x * B, y), A))
        deg_total_1 = max([len(X) + len(Y) - 1 for X, Y in itertools.product(A, B)])
        deg_total_2 = len(list(itertools.product(A, B)))
        ret = [[0] * deg_total_1 for _ in xrange(deg_total_2)]
        deg1 = 0
        for X in A:
            deg2 = 0
            for Y in B:
                for x, y in enumerate(X):
                    for u, v in enumerate(Y):
                        ret[deg1 + deg2][x + u] += y * v
                deg2 += 1
            deg1 += 1
        return s.element_class(s, ret)

    def _neg(s, A):
        return s.element_class(s, map(lambda y: map(lambda x: -x, y), A))

    def _equ(s, A, B):
        if len(A) == 1 and len(A[0]) == 1 and len(B) == 1 and len(B[0]) == 1:
            return A[0][0] == B[0][0]
        return all(
            [
                all([s == t for s, t in zip_longest(x, y, fillvalue=0)])
                for x, y in zip(A, B)
            ]
        )

    def __repr__(s):
        return "%s(%r, %r)" % (s.__class__.__name__, s.field, s.gen_names)

    def __str__(s):
        return "Bivariate Polynomial Ring over %s" % s.field


class BivariatePolynomialElement(RingElement):
    def __init__(s, poly_ring, args):
        RingElement.__init__(s, poly_ring, args)
        if isinstance(args, BivariatePolynomialElement):
            s.coeffs = args.coeffs
        elif hasattr(args, "__iter__"):
            s.coeffs = list(args)
        else:
            s.coeffs = [[args]]
        s.trim()
        s.coeffs = list(map(lambda y: list(map(s.ring.field, y)), s.coeffs))

    def __getitem__(s, idx):
        return s.coeffs[idx]

    def __iter__(s):
        return iter(s.coeffs)

    def __len__(s):
        return len(s.coeffs)

    def __pow__(s, rhs, mod=0):
        if rhs == 0:
            return s.__class__(s.ring, [[1]])
        d = rhs
        if d < 0:
            x = 1 / s
            d = -d
        else:
            x = s
        bits = list(map(int, bin(d)[2:]))[::-1]
        if bits[0]:
            res = x
        else:
            res = s.__class__(s.ring, [[1]])
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
        deg1 = 0
        for y in s.coeffs:
            deg2 = 0
            res2 = []
            for x in y:
                r = ""
                if x != 0:
                    if (x != -1 and x != 1) or deg2 == 0:
                        r += str(x)
                    elif x == -1:
                        r += "-"
                    if deg2 > 0:
                        r += s.ring.gen_names[0]
                        if deg2 > 1:
                            r += "^%d" % deg2
                if r != "":
                    res2 += [r]
                deg2 += 1
            r = "+".join(res2[::-1]).replace("+-", "-")
            if len(r) != 0:
                if deg1 > 0:
                    if deg2 > 2:
                        r = "(" + r + ")"
                    if deg1 == 1:
                        r += s.ring.gen_names[1]
                    else:
                        r += "%s^%d" % (s.ring.gen_names[1], deg1)
                    if r[:2] == "1y":
                        r = r.lstrip("1")
            if r != "":
                res += [r]
            deg1 += 1
        if len(res) == 0:
            return "0"
        return "+".join(res).replace("+-", "-")

    @memoize
    def apply(s, x0, y0):
        """
        Using Horner's method
        """
        ret = 0
        for y in s.coeffs[::-1]:
            t = 0
            for x in y[::-1]:
                t = t * x0 + x
            ret = ret * y0 + t
        return ret

    def degree(s):
        return max([len(x) + i - 1 for i, x in enumerate(s.coeffs)])

    def monic(s):
        return s / s[-1]

    def trim(s):
        res = []
        for y in s.coeffs:
            i = 0
            r = list(y)
            for x in r[::-1]:
                if x != 0:
                    break
                i += 1
            if 0 < i < len(r):
                r = r[:-i]
            elif i == len(r):
                r = [0]
            res += [r]
        i = 0
        for r in res[::-1]:
            if not (len(r) == 1 and r[0] == 0):
                break
            i += 1
        if i < len(res):
            if 0 < i:
                s.coeffs = res[:-i]
            else:
                s.coeffs = res
        elif i == len(res):
            s.coeffs = [[0]]

    def __hash__(s):
        return hash("".join(map(str, s.coeffs)) + str(s.ring))

    def _to_tuple(s, d):
        if isinstance(d, s.__class__):
            return tuple(d)
        elif isinstance(d, tuple):
            return d
        else:
            return ((d,),)
