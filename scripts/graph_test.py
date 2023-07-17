from ecpy import ExtendedFiniteField, EllipticCurve, weil_pairing
from random import randint


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

    rand = [randint(2**31, 2**32) for _ in range(10)]

    print("[+] Weil Pairing: ")
    for x in rand:
        print(weil_pairing(E, P, (x * P).distortion_map(), l))


if __name__ == "__main__":
    main()
