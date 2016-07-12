from Zmod import Zmod


class FiniteField(Zmod):
  """
  Finite Field Class
  """
  def __init__(s, p):
    """
    Constructor of FiniteField
    p should be prime
    """
    Zmod.__init__(s, p)
    s.p = p

  def __str__(s):
    return Zmod.__str__(s, "p")
