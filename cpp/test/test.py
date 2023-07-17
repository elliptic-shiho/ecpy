from __future__ import print_function
from ecpy.native import *
import time
import sys

ac_count = 0
wa_count = 0


def _assert(a, b, msg, cond):
    global ac_count, wa_count
    msg = msg.ljust(16)
    print(("[+] %s..." % (msg)).ljust(30), end=" ")
    var = {"a": a, "b": b}
    if eval("a %s b" % cond, var):
        r = repr(b)[:64]
        if repr(b) != r:
            r += "..."
        print("\x1b[33m[  OK  ]\x1b[0m %s" % (r,))
        ac_count += 1
    else:
        print("\x1b[31m[ Fail ]\x1b[0m Expected: %r, Result: %r" % (b, a))
        wa_count += 1


def assert_neq(a, b, m):
    _assert(a, b, m, "!=")


def assert_eq(a, b, m):
    _assert(a, b, m, "==")


def main():
    F = FF(31)
    x = FF_elem(25)
    y = FF_elem(10)
    t = FF_elem(0)
    F.add(t, x, y)
    assert_eq(t.to_python(), 4, "x+y")
    F.sub(t, x, y)
    assert_eq(t.to_python(), 15, "x+y")
    F.mul(t, x, y)
    assert_eq(t.to_python(), 2, "x*y")
    F.div(t, x, y)
    assert_eq(t.to_python(), 18, "x/y")
    F.pow(t, x, y)
    assert_eq(t.to_python(), 25, "x^y")

    F = EF(41, "x^2+x+1")
    x = EF_elem(15, 25)
    y = EF_elem(39, 10)
    t = EF_elem(0, 0)

    F.add(t, x, y)
    assert_eq(t.to_python(), (13, 35), "x+y=13+35w")
    F.sub(t, x, y)
    assert_eq(t.to_python(), (17, 15), "x-y=17+15w")
    F.mul(t, x, y)
    assert_eq(t.to_python(), (7, 14), "x+y=7+14w")
    F.div(t, x, y)
    assert_eq(t.to_python(), (29, 5), "x+y=29+5w")
    print(
        "[+] %d Test(s) finished. %d Test(s) success, %d Test(s) fail."
        % (ac_count + wa_count, ac_count, wa_count)
    )

    p = 17
    F = FF(p)
    E = EC(F, 0, 1)
    P = EC_elem(E, 1, 6)
    Q = EC_elem(E, 14, 5)
    t = EC_elem(E, 0, 0)

    E.add(t, P, Q)
    assert_eq(t.to_python(), (4, 10, 4), "P+Q")
    E.add(t, P, P)
    assert_eq(t.to_python(), (1, 13, 11), "P+P")
    E.mul(t, P, 2)
    assert_eq(t.to_python(), (1, 13, 11), "2*P")
    E.mul(t, P, 4)
    assert_eq(t.to_python(), (15, 3, 7), "4*P")
    E.mul(t, P, 9)
    assert_eq(t.to_python(), (0, 1, 0), "#P = 9")

    p = 631
    m = 5
    F = FF(p)
    E = EC(F, 30, 34)
    P = EC_elem(E, 36, 60)
    Q = EC_elem(E, 121, 387)
    S = EC_elem(E, 0, 36)
    t = FF_elem(0)
    weil_pairing(t, E, P, Q, S, m)
    assert_eq(t.to_python(), 242, "weil_pairing")
    tate_pairing(t, E, P, Q, m, 1)
    assert_eq(t.to_python(), 279, "tate_pairing")

    sys.exit(wa_count)


if __name__ == "__main__":
    main()
