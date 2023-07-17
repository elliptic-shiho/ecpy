from __future__ import print_function
from ecpy import *
from random import randint
from six.moves import xrange
import sys

ac_count = 0
wa_count = 0


def _assert(a, b, msg, cond):
    global ac_count, wa_count
    msg = msg.ljust(16)
    print(("[+] %s..." % (msg)).ljust(30), end=" ")
    var = {"a": a, "b": b}
    if eval("a %s b" % cond, var):
        print("\x1b[33m[  OK  ]\x1b[0m %r" % (b,))
        ac_count += 1
    else:
        print("\x1b[31m[ Fail ]\x1b[0m Expected: %r, Result: %r" % (b, a))
        wa_count += 1


def assert_neq(a, b, m):
    _assert(a, b, m, "!=")


def assert_eq(a, b, m):
    _assert(a, b, m, "==")


def test():
    F = FiniteField(101)
    x = F(2)
    y = F(203)
    assert_neq(x, 1, "x != 1")
    assert_eq(x, 2, "x == 2")

    F = FiniteField(5)
    x = F(3)
    y = F(7)  # = 2
    print(F)
    print("[+] x, y = %s, %s" % (x, y))
    assert_eq(x + y, F(0), "x+y == F(0)")
    assert_eq(x + y, 0, "x+y == 0")
    assert_eq(x - y, 1, "x-y == 1")
    assert_eq(x * y, 1, "x*y == 1")
    x = F(2)
    y = F(3)
    print("[+] x, y = %s, %s" % (x, y))
    # commutive!
    assert_eq(1 / x, y, "1/x == y")
    assert_eq(util.modinv(x.x, F.p), y, "modinv(x) == y")
    assert_eq(1 / y, x, "1/y == x")

    assert_eq(x**3, y, "x^3 == y")

    assert_eq(util.crt([3, 4], [4, 9]), 31, "CRT Test")
    # assert_eq(util.crt([7, 13], [12, 18]), 31, "CRT Test 2")

    assert_eq(F.order(), 4, "|F| = 4")
    assert_eq(x.order(), 4, "|x| = 4")

    F = FiniteField(17)
    E = EllipticCurve(F, 1, 0)

    P = E(1, 6)
    Q = E(11, 4)
    print("P, Q = %r, %r" % (P, Q))
    assert_eq(P + Q, E(3, 8), "P+Q")
    assert_eq(P + P, E(0, 0), "P+P")
    assert_eq(P * 2, E(0, 0), "P*2")
    assert_eq(2 * P, E(0, 0), "2*P")
    assert_eq(P.order(), 4, "|P| = 4")

    print("Random Test: ")
    i = 0
    while i < 10:
        while True:
            r = randint(-50, 50)
            if r != 0:
                break
        print("[+] random 1 = %d" % r)
        assert_eq((util.modinv(r, 101) * r) % 101, 1, "modinv")
        while True:
            q = randint(-50, 50)
            if q != 0:
                break
        print("[+] random 2 = %d" % q)
        assert_eq(r * (q * P), q * (r * P), "ECDH test")
        i += 1

    # The arithmetic of elliptic curves: p.397 example of miller algorithm
    F = FiniteField(631)
    E = EllipticCurve(F, 30, 34)
    m = 5

    P = E(36, 60)
    Q = E(121, 387)
    S = E(0, 36)

    print("P, Q, S = %r, %r, %r" % (P, Q, S))
    assert_eq(E.embedding_degree(m), 1, "embedding degree")
    assert_eq(miller(E, P, Q + S, m), 103, "miller(P, Q+S)")
    assert_eq(miller(E, P, S, m), 219, "miller(P, S)")
    assert_eq(miller(E, Q, P - S, m), 284, "miller(Q, P-S)")
    assert_eq(miller(E, Q, -S, m), 204, "miller(Q, -S)")
    assert_eq(weil_pairing(E, P, Q, m, S), 242, "weil_pairing")
    assert_eq(tate_pairing(E, P, Q, m, 1), 279, "tate_pairing")
    g = tate_pairing(E, P, Q, m)
    print("[+] g = %s" % g)
    assert_eq(tate_pairing(E, 2 * P, Q, m), g**2, "e(2P, Q) == g^2")
    assert_eq(tate_pairing(E, P, 2 * Q, m), g**2, "e(P, 2Q) == g^2")
    assert_eq(tate_pairing(E, P, Q, m) ** 2, g**2, "e(P, Q)^2 == g^2")

    print("[+] SSSA-Attack")
    F = FiniteField(16857450949524777441941817393974784044780411511252189319)

    A = 16857450949524777441941817393974784044780411507861094535
    B = 77986137112576

    E = EllipticCurve(F, A, B)

    """
  print "Random Point"
  for x in xrange(10):
    print E.random_point()
  """

    P = E(
        5732560139258194764535999929325388041568732716579308775,
        14532336890195013837874850588152996214121327870156054248,
    )
    Q = E(
        2609506039090139098835068603396546214836589143940493046,
        8637771092812212464887027788957801177574860926032421582,
    )

    assert_eq(
        SSSA_Attack(F, E, P, Q),
        6418297401790414611703852603267852625498215178707956450,
        "SSSA-Attack",
    )

    p = 0xD3CEEC4C84AF8FA5F3E9AF91E00CABACAAAECEC3DA619400E29A25ABECECFDC9BD678E2708A58ACB1BD15370ACC39C596807DAB6229DCA11FD3A217510258D1B
    A = 0x95FC77EB3119991A0022168C83EEE7178E6C3EEAF75E0FDF1853B8EF4CB97A9058C271EE193B8B27938A07052F918C35ECCB027B0B168B4E2566B247B91DC07
    B = 0x926B0E42376D112CA971569A8D3B3EDA12172DFB4929AEA13DA7F10FB81F3B96BF1E28B4A396A1FCF38D80B463582E45D06A548E0DC0D567FC668BD119C346B2
    Gx = 0xCF634030986CF41C1ADD87E71D638B9CC723C764059CF4C9B8ED2A0AAF5D51DC770372503EBFAAD746AB9220E992C09822916978226465AD31D354A3EFEE51DA
    Gy = 0x65EAAD8848B2787103FCE02358B45D8A61420031989EB6B4B70D82FE20D85583AE542EB8F76749DC640B0F13F682228819B8B2F04BD7A5A17A4C675540FE1C90
    Px = 10150325274093651859575658519947563789222194633356867789068177057343771571940302488270622886585658965620106459791565259790154958179860547267338437952379763
    Py = 6795014289013853849339410895464797184780777251924203530417684718894057583288011725702609805686960505075072642102076744937056900144377846048950215257629102

    F = FiniteField(p)
    E = EllipticCurve(F, A, B)
    G = E(Gx, Gy)
    P = E(Px, Py)
    assert_eq(
        SSSA_Attack(F, E, G, P),
        0x746A6374667B6F6F70735F656C6C31707431635F6375727665735F525F683472647D,
        "SSSA-Attack (TJCTF 2016 Crypto 200: curvature2)\n",
    )
    z = CC(1, 2)  # 1+2i
    w = CC(5, 1)  # 5+i
    print("z, w = %r, %r" % (z, w))
    assert_eq(z + w, CC(6, 3), "z+w")
    assert_eq(z - w, CC(-4, 1), "z-w")
    assert_eq(z * w, CC(3, 11), "z*w")
    assert_eq(z / w, CC(0.2692307692307693, 0.34615384615384615), "z/w")

    F = ExtendedFiniteField(59)
    a = F(0, 1)
    E = EllipticCurve(F, 1, 0)
    P = E(25, 30)
    assert_eq(tuple(P), (25, 30, 1), "extended field EC")
    Q = P.distortion_map()
    assert_eq(tuple(Q), (F(34), F(0, 30), 1), "extended field EC 2")

    assert_eq(Q.distortion_map(), P, "distortion map")

    l = 56453
    m = l
    p = l * 6 - 1
    F = ExtendedFiniteField(p, "x^2+x+1")
    print("Random Test 2:")
    for x in xrange(10):
        r1 = randint(0, p)
        r2 = randint(0, p)
        r = F(r1, r2)
        print("[+] r = %s" % r)
        assert_eq(r ** (p**2), r, "r^(p^2) == r")

    E = EllipticCurve(F, 0, 1)
    P = E(3, 1164)
    print(P)
    print(P.distortion_map())

    g = symmetric_weil_pairing(E, P, P, m)
    print("[+] g = %s" % g)

    assert_eq(symmetric_weil_pairing(E, P, 2 * P, m), g**2, "e(P, 2P) == g^2")
    assert_eq(symmetric_weil_pairing(E, 2 * P, P, m), g**2, "e(2P, 2P) == g^2")
    assert_eq(symmetric_weil_pairing(E, P, P, m) ** 2, g**2, "e(P, P)^2 == g^2")

    g = symmetric_tate_pairing(E, P, P, m)
    print("[+] g = %s" % g)

    assert_eq(symmetric_tate_pairing(E, P, 2 * P, m), g**2, "e(P, 2P) == g^2")
    assert_eq(symmetric_tate_pairing(E, 2 * P, P, m), g**2, "e(2P, 2P) == g^2")

    assert_eq(
        F(53521, 219283) / F(297512, 101495), F(333099, 288028), "r/prev_r test: 0"
    )
    assert_eq(F(281317, 98371) / F(53521, 219283), F(323815, 46359), "r/prev_r test: 1")
    assert_eq(F(31851, 95658) / F(281317, 98371), F(5298, 9638), "r/prev_r test: 2")
    assert_eq(F(92937, 215632) / F(31851, 95658), F(278130, 175879), "r/prev_r test: 3")
    assert_eq(
        F(61703, 173508) / F(92937, 215632), F(189715, 176788), "r/prev_r test: 4"
    )
    assert_eq(F(80979, 72727) / F(61703, 173508), F(15407, 212022), "r/prev_r test: 5")
    assert_eq(F(311516, 184895) / F(80979, 72727), F(225531, 44087), "r/prev_r test: 6")
    assert_eq(
        F(326035, 114920) / F(311516, 184895), F(213234, 100495), "r/prev_r test: 7"
    )
    assert_eq(
        F(294922, 165746) / F(326035, 114920), F(113566, 200451), "r/prev_r test: 8"
    )
    assert_eq(
        F(73542, 195813) / F(294922, 165746), F(201397, 252614), "r/prev_r test: 9"
    )
    assert_eq(1 / F(338714, 3), F(37635, 188176), "division by b = -a")

    assert_eq(F(302128, 326350) * F(39563, 131552), F(151684, 28719), "multiple test")

    assert_eq(
        miller(E, P, P.distortion_map(), m), F(239139, 508), "miller function check"
    )

    a = F(234687, 190012)
    b = F(218932, 251221)

    print("[+] a = %s" % a)
    print("[+] b = %s" % b)

    assert_eq(a + b, F(114902, 102516), "a+b")
    assert_eq(a - b, F(15755, 277508), "a-b")
    assert_eq(a * b, F(217278, 89209), "a*b")
    assert_eq(a / b, F(167345, 81997), "a/b")
    assert_eq(a * a * a * a, F(31723, 160374), "a*a*a*a")
    assert_eq(a**4, F(31723, 160374), "a^4")
    assert_eq(a * a * b * b * a * a * b, b * a * a * b * b * a * a, "a^4")

    # E, F, l = gen_supersingular_ec()
    l = 670879118046860960563
    p = 4025274708281165763377
    F = ExtendedFiniteField(p, "x^2+x+1")
    E = EllipticCurve(F, 0, 1)
    # P = find_point_by_order(E, l)
    P = E(7, 658176732497617012595, 1)
    Q = P.distortion_map()
    g = tate_pairing(E, P, Q, l)
    print(E)
    for x in xrange(10):
        a = randint(2**15, 2**16)
        b = randint(2**15, 2**16)
        print("Random Pairing Test: a = %d, b = %d" % (a, b))
        gab = g ** (a * b)
        assert_eq(tate_pairing(E, a * P, b * Q, l), gab, "e(aP, bQ)")
        assert_eq(tate_pairing(E, a * P, Q, l) ** b, gab, "e(aP, Q)^b")
        assert_eq(tate_pairing(E, b * P, Q, l) ** a, gab, "e(bP, Q)^a")
        assert_eq(tate_pairing(E, P, a * Q, l) ** b, gab, "e(P, aQ)^b")
        assert_eq(tate_pairing(E, P, b * Q, l) ** a, gab, "e(P, bQ)^a")

    assert_eq(
        util.prime_factorization(12345678),
        {2: 1, 3: 2, 47: 1, 14593: 1},
        "prime factor1",
    )
    assert_eq(util.prime_factorization(12345), {3: 1, 5: 1, 823: 1}, "prime factor2")
    assert_eq(util.euler_phi(12345), 6576, "phi(12345)")

    for name in [
        "secp192k1",
        "secp192r1",
        "secp224k1",
        "secp224r1",
        "secp256k1",
        "secp256r1",
        "secp384r1",
        "secp521r1",
    ]:
        F, E, G, n = EllipticCurveRepository(name)
        assert_eq(G * n, E.O, name)

    # secp224k1 bug (Issue #3)
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFE56D
    a = 0
    b = 5
    F = FiniteField(p)
    E = EllipticCurve(F, a, b)
    Gx = 0xA1455B334DF099DF30FC28A169A467E9E47075A90F7E650EB6B7A45C
    Gy = 0x7E089FED7FBA344282CAFBD6F7E319F7C0B0BD59E2CA4BDB556D61A5
    G = E(Gx, Gy)
    n = 0x010000000000000000000000000001DCE8D2EC6184CAF0A971769FB1F7
    assert_eq(n * G, E.O, "Issue #3")

    PR = UnivariatePolynomialRing(RR, "X")
    X = PR.gen()
    assert_eq(str(X), "X", "Polynomial ring generator name")

    poly_test = lambda p, q: all(map(lambda x: x[0] == x[1], zip(p, q)))

    assert_eq(poly_test(X + 1, [1, 1]), True, "x+1")
    assert_eq(poly_test(X - 1, [-1, 1]), True, "x-1")
    assert_eq(poly_test(X * X, [0, 0, 1]), True, "x^2")
    assert_eq(poly_test((X**2 + 1) - X, [1, -1, 1]), True, "(x^2+1) - x")
    assert_eq(
        poly_test((X**3 + 8 * X**2 + x + 8) / (X**2 + 1), [8, 1]),
        True,
        "(x^3+8x^2+x+8)/(x^2+1)",
    )
    assert_eq(
        poly_test((X**2 + 2 * X - 2) % (X**2 + 1), [-3, 2]),
        True,
        "(x^2+2x-2) % (x^2+1)",
    )

    print("=== Test for Univariate Polynomial Ring over Z/5Z ===")

    PR = UnivariatePolynomialRing(Zmod(5))
    x = PR.gen()
    P = x**2 + 1
    R = x**3 - 3
    assert_eq(poly_test(P * R, [2, 0, 2, 1, 0, 1]), True, "P * R")
    assert_eq(poly_test((P * R) / R, P), True, "(P * R) / R")
    assert_eq(
        poly_test((P * R) / (3 * x**2 + 3), [4, 0, 0, 2]), True, "(P * R) / 3x^2 + 3"
    )

    PR = BivariatePolynomialRing(RR, ["x", "y"])
    print(PR)
    x, y = PR.gens()
    assert_eq(str(x), "x", "str(x)")
    assert_eq(str(y), "y", "str(y)")
    assert_eq(str(x + y), "x+y", "str(x+y)")
    assert_eq(str(y * (x + y)), "xy+y^2", "str(y*(x+y))")
    print(y**2 - x**3 + 3 * x - 1)

    x = ZZ(0xCAFEBABE)
    y = ZZ(0xDEADBEEF)
    assert_eq(x + y, 7141620141, "Integer Addition")
    assert_eq(x - y, -330236977, "Integer Subtract")
    assert_eq(x * y, 12723420444339690338, "Integer Multiplication")

    QR = QuotientRing(ZZ, 5)
    print(QR)
    x = QR(3)
    y = QR(7)
    assert_eq(x, 3, "3 mod 5")
    assert_eq(y, 2, "7 mod 5")
    assert_eq(x + y, 0, "x+y")
    assert_eq(x - y, 1, "x-y")
    assert_eq(x * y, 1, "x*y")
    assert_eq(x / y, 4, "x/y")

    PR = UnivariatePolynomialRing(ZZ)
    x = PR.gen()
    QR = QuotientRing(PR, x**2 + 1)
    print(QR)
    a = QR(1)
    b = QR(x)
    assert_eq(a, 1, "a")
    assert_eq(b, x, "b")
    assert_eq(a + b, x + 1, "a+b")
    assert_eq(a * b, x, "a*b")
    assert_eq(b * b, -1, "b^2")
    assert_eq(b * b + a, 0, "b^2 + a")

    print(
        "[+] %d Test(s) finished. %d Test(s) success, %d Test(s) fail."
        % (ac_count + wa_count, ac_count, wa_count)
    )
    sys.exit(wa_count)


if __name__ == "__main__":
    test()
