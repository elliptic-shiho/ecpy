from .library import *
from .FF import FF, FF_elem
from .EF import EF, EF_elem
import ast


class EC(NativeProxy):
    def __init__(s, base, a, b):
        s.base = base
        s.a = a
        s.b = b
        if isinstance(base, FF):
            ptr = lib.EC_FF_create(to_char_ptr(str(a)), to_char_ptr(str(b)), base.ptr)
            tostring_func = lib.EC_FF_to_string
            del_func = lib.EC_FF_delete
            s.add_func = lib.EC_FF_add
            s.sub_func = lib.EC_FF_sub
            s.mul_func = lib.EC_FF_mul
            s.type = 1
        elif isinstance(base, EF):
            ptr = lib.EC_EF_create(
                to_char_ptr(str(a)),
                to_char_ptr(str(b)),
                base.ptr,
                to_char_ptr(base.poly),
            )
            tostring_func = lib.EC_EF_to_string
            del_func = lib.EC_EF_delete
            s.add_func = lib.EC_EF_add
            s.sub_func = lib.EC_EF_sub
            s.mul_func = lib.EC_EF_mul
            s.type = 2
        NativeProxy.__init__(s, ptr, tostring_func, del_func)

    def add(s, ret, a, b):
        assert (
            isinstance(ret, EC_elem)
            and isinstance(a, EC_elem)
            and isinstance(b, EC_elem)
        )
        s.add_func(s.ptr, ret.ptr, a.ptr, b.ptr)

    def sub(s, ret, a, b):
        assert (
            isinstance(ret, EC_elem)
            and isinstance(a, EC_elem)
            and isinstance(b, EC_elem)
        )
        s.sub_func(s.ptr, ret.ptr, a.ptr, b.ptr)

    def mul(s, ret, a, b):
        assert isinstance(ret, EC_elem) and isinstance(a, EC_elem)
        s.mul_func(s.ptr, ret.ptr, a.ptr, to_char_ptr(str(b)))

    def div(s, ret, a, b):
        raise NotImplementedError()

    def pow(s, ret, a, b):
        raise NotImplementedError()


class EC_elem(NativeProxy):
    def __init__(s, curve, x, y, z=1):
        from six import integer_types

        def conv(x):
            if s.curve.type == 1:
                if isinstance(x, tuple):
                    return FF_elem(x[0])
                else:
                    return FF_elem(x)
            elif s.curve.type == 2:
                if isinstance(x, tuple):
                    return EF_elem(x[0], x[1])
                else:
                    return EF_elem(x, 0)

        assert isinstance(curve, EC)
        s.x = x
        s.y = y
        s.z = z
        s.curve = curve
        s.base = curve.base
        if isinstance(x, integer_types + (tuple,)):
            x = conv(x)
        if isinstance(y, integer_types + (tuple,)):
            y = conv(y)
        if isinstance(z, integer_types + (tuple,)):
            z = conv(z)

        if s.curve.type == 1:
            ptr = lib.EC_elem_FF_create(x.ptr, y.ptr, z.ptr)
            tostring_func = lib.EC_elem_FF_to_string
            del_func = lib.EC_elem_FF_delete
        elif s.curve.type == 2:
            ptr = lib.EC_elem_EF_create(x.ptr, y.ptr, z.ptr)
            tostring_func = lib.EC_elem_EF_to_string
            del_func = lib.EC_elem_EF_delete

        NativeProxy.__init__(s, ptr, tostring_func, del_func)

    def to_python(s):
        r = str(s).lstrip("EC_elem").replace("EF_elem", "").replace("FF_elem", "")
        return tuple(ast.literal_eval(r))
