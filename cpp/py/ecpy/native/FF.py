from .library import lib, NativeProxy, c_char_p
from ctypes import c_void_p


class FF(NativeProxy):
    _FF_create = lib.FF_create
    _FF_create.restype = c_void_p

    def __init__(s, p):
        s.p = p
        NativeProxy.__init__(
            s,
            FF._FF_create(c_char_p(str(p).encode("us-ascii"))),
            lib.FF_to_string,
            lib.FF_delete,
        )

    def add(s, ret, a, b):
        assert isinstance(ret, FF_elem)
        assert isinstance(a, FF_elem)
        assert isinstance(b, FF_elem)
        lib.FF_add(s.ptr, ret.ptr, a.ptr, b.ptr)

    def sub(s, ret, a, b):
        assert isinstance(ret, FF_elem)
        assert isinstance(a, FF_elem)
        assert isinstance(b, FF_elem)
        lib.FF_sub(s.ptr, ret.ptr, a.ptr, b.ptr)

    def mul(s, ret, a, b):
        assert isinstance(ret, FF_elem)
        assert isinstance(a, FF_elem)
        assert isinstance(b, FF_elem)
        lib.FF_mul(s.ptr, ret.ptr, a.ptr, b.ptr)

    def div(s, ret, a, b):
        assert isinstance(ret, FF_elem)
        assert isinstance(a, FF_elem)
        assert isinstance(b, FF_elem)
        lib.FF_div(s.ptr, ret.ptr, a.ptr, b.ptr)

    def pow(s, ret, a, b):
        assert isinstance(ret, FF_elem) and isinstance(a, FF_elem)
        lib.FF_pow(s.ptr, ret.ptr, a.ptr, c_char_p(str(b).encode("us-ascii")))


class FF_elem(NativeProxy):
    _FF_elem_create = lib.FF_elem_create
    _FF_elem_create.restype = c_void_p

    def __init__(s, v):
        s.v = v
        NativeProxy.__init__(
            s,
            FF_elem._FF_elem_create(c_char_p(str(v).encode("us-ascii"))),
            lib.FF_elem_to_string,
            lib.FF_elem_delete,
        )

    def to_python(s):
        return int(str(s))
