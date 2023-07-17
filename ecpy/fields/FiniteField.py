from .Zmod import Zmod, ZmodElement


class FiniteField(Zmod):
    """
    Finite Field Class
    """

    def __init__(s, p):
        """
        Constructor of FiniteField
        p should be prime
        """
        s.n = s.p = p
        s.element_class = ZmodElement

    def __str__(s):
        return Zmod.__str__(s, "p")
