from .library import *
import ast


class EF(NativeProxy):
    def __init__(s, p, poly):
        s.p = p
        s.poly = poly
        NativeProxy.__init__(
            s,
            lib.EF_create(to_char_ptr(str(p)), to_char_ptr(poly)),
            lib.EF_to_string,
            lib.EF_delete,
        )

    def add(s, ret, a, b):
        assert (
            isinstance(ret, EF_elem)
            and isinstance(a, EF_elem)
            and isinstance(b, EF_elem)
        )
        lib.EF_add(s.ptr, ret.ptr, a.ptr, b.ptr)

    def sub(s, ret, a, b):
        assert (
            isinstance(ret, EF_elem)
            and isinstance(a, EF_elem)
            and isinstance(b, EF_elem)
        )
        lib.EF_sub(s.ptr, ret.ptr, a.ptr, b.ptr)

    def mul(s, ret, a, b):
        assert (
            isinstance(ret, EF_elem)
            and isinstance(a, EF_elem)
            and isinstance(b, EF_elem)
        )
        lib.EF_mul(s.ptr, ret.ptr, a.ptr, b.ptr)

    def div(s, ret, a, b):
        assert (
            isinstance(ret, EF_elem)
            and isinstance(a, EF_elem)
            and isinstance(b, EF_elem)
        )
        lib.EF_div(s.ptr, ret.ptr, a.ptr, b.ptr)

    def pow(s, ret, a, b):
        assert isinstance(ret, EF_elem) and isinstance(a, EF_elem)
        lib.EF_pow(s.ptr, ret.ptr, a.ptr, to_char_ptr(str(b)))


class EF_elem(NativeProxy):
    def __init__(s, u, v):
        NativeProxy.__init__(
            s,
            lib.EF_elem_create(to_char_ptr(str(u)), to_char_ptr(str(v))),
            lib.EF_elem_to_string,
            lib.EF_elem_delete,
        )

    def to_python(s):
        r = str(s).lstrip("EF_elem")
        return tuple(ast.literal_eval(r))
