from __future__ import print_function
from six.moves import xrange
from ecpy import *
from timeit import timeit
from random import randint

test_param = []


def init_test(t, *args):
    global test_param
    cond = {"weil": weil_pairing, "tate": tate_pairing}
    test_param = [cond[t], args[0]]


def do_test():
    global test_param
    test_param[0](*test_param[1])


def show_results(name, result, count):
    per_pass = 1000000 * (result / count)
    print("%s: %.2f usec/pass" % (name, per_pass))


def main():
    l = 2201426263
    p = 6 * l - 1
    F = ExtendedFiniteField(p, "x^2+x+1")
    E = EllipticCurve(F, 0, 1)
    i = 3
    while True:
        y = E.get_corresponding_y(i)
        if y != None:
            P = E(i, y)
            if (l * P).is_infinity():
                break
        i += 1
    print(P)

    rand = [randint(2**31, 2**32) for _ in xrange(10)]
    count = 20

    weil_time = []
    print("[+] Weil Pairing: ")
    for x in rand:
        r = timeit(
            "bench_pairing.do_test()",
            setup="import bench_pairing; from ecpy import EllipticCurve, ExtendedFiniteField; bench_pairing.init_test('weil', [%r, %r, %r.distortion_map(), %r])"
            % (E, P, x * P, l),
            number=count,
        )
        weil_time += [r]
        show_results("weil", r, count)

    print("[+] Tate Pairing: ")
    tate_time = []
    for x in rand:
        r = timeit(
            "bench_pairing.do_test()",
            setup="import bench_pairing; from ecpy import EllipticCurve, ExtendedFiniteField; bench_pairing.init_test('tate', [%r, %r, %r.distortion_map(), %r, 2])"
            % (E, P, x * P, l),
            number=count,
        )
        tate_time += [r]
        show_results("tate", r, count)

    print("=" * 64)
    show_results("weil", sum(weil_time), count * 10)
    show_results("tate", sum(tate_time), count * 10)


if __name__ == "__main__":
    main()
