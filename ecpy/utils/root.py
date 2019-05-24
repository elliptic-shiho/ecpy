from __future__ import division
import random


def _find_power_divisor(base, x, modulo=None):
  k = 0
  m = base
  while x % m == 0:
    k += 1
    m = pow(m * base, 1, modulo)
  return k


def _find_power(power_base, x, crib, modulo=None):
  k = 1
  r = power_base
  while pow(x, r, modulo) != crib:
    k += 1
    r *= power_base
  return k


def modular_square_root(a, m, force=False):
  """
  Calculate Quadratic Residue
  Args:
    a: The Target
    m: Modulus
  Returns:
    A Quadratic Residue of a modulo m
  """
  from ecpy.fields.ExtendedFiniteField import ExtendedFiniteFieldElement
  if isinstance(a, ExtendedFiniteFieldElement) and not force:
    return modular_square_root_extended(a)
  else:
    return __modular_square_root(int(a), m)


def __modular_square_root(a, m):
  from ecpy.utils import is_prime, legendre_symbol, modinv
  if is_prime(m):
    if legendre_symbol(a, m) == -1:
      return []
    # Tonelli-Shanks Algorithm
    if m % 4 == 3:
      r = pow(a, (m + 1) // 4, m)
      return [r, m - r]
    s = _find_power_divisor(2, m - 1)
    q = (m - 1) // 2**s
    z = 0
    while legendre_symbol(z, m) != -1:
      z = random.randint(1, m)
    c = pow(z, q, m)
    r = pow(a, (q + 1) // 2, m)
    t = pow(a, q, m)
    l = s
    while True:
      if t % m == 1:
        assert (r ** 2) % m == a
        return [r, m - r]
      i = _find_power(2, t, 1, m)
      power = l - i - 1
      if power < 0:
        power = modinv(2**-power, m)
      else:
        power = 2**power
      b = pow(c, power, m)
      r = (r * b) % m
      t = (t * (b**2)) % m
      c = pow(b, 2, m)
      l = i
  if m == 2:
    return a
  if m % 4 == 3:
    r = pow(a, (m + 1) // 4, m)
    return [r, m - r]
  if m % 8 == 5:
    v = pow(2 * a, (m - 5) // 8, m)
    i = pow(2 * a * v, 2, m)
    r = a * v * (i - 1) % m
    return [r, m - r]
  if m % 8 == 1:
    e = _find_power_divisor(2, m - 1)
    q = (m - 1) // 2**e
    z = 1
    while pow(z, 2**(e - 1), m) == 1:
      x = random.randint(1, m)
      z = pow(x, q, m)
    y = z
    r = e
    x = pow(a, (q - 1) // 2, m)
    v = a * x % m
    w = v * x % m
    while True:
      if w == 1:
        return [v, m - v]
      k = _find_power(2, w, 1, m)
      d = pow(y, 2**(r - k - 1), m)
      y = pow(d, 2, m)
      r = k
      v = d * v % m
      w = w * y % m


def extended_legendre_symbol(a):
  """
  Legendre Symbol on the Extended Field
  Args:
    a: The Target
  Returns:
    Legendre Symbol of a
  """
  from six.moves import xrange
  from ecpy.utils import legendre_symbol
  m = a.field.degree()
  p = a.field.p
  b = pow(a, sum([p**i for i in xrange(0, m)]), p)
  return legendre_symbol(b, p)


def modular_square_root_extended(x):
  """
  Calculate Quadratic Residue on Extended Field
  Args:
    x: The Target
  Returns:
    A square root of x
  """
  if extended_legendre_symbol(x) != 1:
    return []
  a = x

  m = a.field.degree()
  p = a.field.p
  q = p**(m // 2)
  if m % 2 == 0:
    if pow(q, m // 2, 4) == 1:
      c0 = 1
      while c0 == 1:
        c = x.field(random.randint(0, q**2), random.randint(0, q**2))
        c0 = extended_legendre_symbol(c)
      d = pow(c, (q - 1) // 2)
      e = 1 // (c * d)
      f = (c * d)**2
      b = pow(a, (q - 1) // 4)
      b2 = b**2
      a0 = (b2) ** q * b2
      if a0 == -1:
        return []
      R = b2 * a
      S = b**q
      if b**(q + 1) != 1:
        S *= e
        R *= f
      x0 = modular_square_root(R, q, force=True)[0]
      return [x0 * S]
    elif pow(q, m // 2, 4) == 3:
      pass
    else:
      pass
  else:
    raise NotImplementedError("If you want to use function, please implement "
                              "Shanks, Atkin, Kong et al, Mueller, Tonelli-"
                              "Shanks, ... algorithm.")


def cubic_root(x):
  """
  Calculate Cubic Residue
  Args:
    x: The Target
  Returns:
    A Cubic Residue of x
  """
  F = x.field
  p = F.p
  m = F.degree()
  rho = 1
  pm = p**m
  r = (pm - 1) // 3
  rho = 1
  while rho ** r == 1:
    rho = F(random.randint(1, p - 1), random.randint(1, p - 1))
  t = 1
  while True:
    r = 3**t
    if (pm - 1) % r == 0:
      t += 1
    else:
      t -= 1
      s = (pm - 1) // 3**t
      if (s + 1) % 3 == 0:
        l = (s + 1) // 3
        break
      elif (s - 1) % 3 == 0:
        l = (s - 1) // 3
        break
  a = rho ** s
  a_ = rho ** (3**(t - 1) * s)
  b = x ** s
  h = 1
  i = 1
  while i < t:
    d = b ** (3**(t - 1 - i))
    if d == 1:
      k = 0
    elif d == a_:
      k = 2
    else:
      k = 1
    b = b * (a * a * a)**k
    h = h * a**k
    a = a * a * a
    i += 1
  r = x**l * h
  if s == 3 * l + 1:
    r = 1 // r
  return r
