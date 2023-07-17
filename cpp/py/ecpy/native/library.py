from ctypes import *

lib = cdll.LoadLibrary("libecpy_native.so")


def to_char_ptr(string):
    return c_char_p(string.encode("us-ascii"))


class NativeProxy(object):
    def __init__(s, object_ptr, tostring_func, del_func):
        s.ptr = c_void_p(object_ptr)
        s.tostring_func = tostring_func
        s.del_func = del_func

    def __to_string(s, bufsize):
        b = create_string_buffer(bufsize)
        s.tostring_func(s.ptr, b, bufsize)
        b = b.value
        if len(b) == 0:  # not enough buffer size
            return s.__to_string(2 * bufsize)
        return b

    def __str__(s):
        return str(s.__to_string(1024).decode("us-ascii"))

    def __del__(s):
        s.del_func(s.ptr)
