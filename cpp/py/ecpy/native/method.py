from .library import *


def miller(ret, E, P, Q, m):
    cond = {1: lib.EC_FF_miller, 2: lib.EC_EF_miller}
    cond[E.type](ret.ptr, E.ptr, P.ptr, Q.ptr, to_char_ptr(str(m)))


def weil_pairing(ret, E, P, Q, S, m):
    cond = {1: lib.EC_FF_weil_pairing, 2: lib.EC_EF_weil_pairing}
    cond[E.type](ret.ptr, E.ptr, P.ptr, Q.ptr, S.ptr, to_char_ptr(str(m)))


def tate_pairing(ret, E, P, Q, m, k):
    cond = {1: lib.EC_FF_tate_pairing, 2: lib.EC_EF_tate_pairing}
    cond[E.type](ret.ptr, E.ptr, P.ptr, Q.ptr, to_char_ptr(str(m)), int(k))
