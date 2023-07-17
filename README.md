# ecpy
[![Code Climate](https://codeclimate.com/github/elliptic-shiho/ecpy/badges/gpa.svg)](https://codeclimate.com/github/elliptic-shiho/ecpy)
[![Build Status](http://travis-ci.org/elliptic-shiho/ecpy.svg?branch=master)](https://travis-ci.org/elliptic-shiho/ecpy)

Elliptic-Curve Cryptography Library (Implemented by Python)

# Features
* Calculation between elliptic curve points
  - Elliptic curve\: `y^2 = x^3 + ax + b` over Field `K` (`char(K) != 2, 3`).
  - You can choose K from below:
    + Complex Field, Rational Field, FiniteField and Extended Finite Field(Irreducible polynomial: `x^2 + 1` or `x^2 + x + 1`)
* Weil/Tate Pairing
  - Distortion Map, symmetric-pairing functions
* SSSA-Attack Implementation
  - [Sample](https://gist.github.com/elliptic-shiho/40d42dbab87065e06d6c473ef93e244e) (CTF Challenge Writeup)
* Example Applications
  - `example/id_based_encryption.py`: Implementation of Boneh-Franklin's Identity-based Encryption Scheme
  - `example/boneh_lynn_shacham_short_signature.py`: Implementation of Boneh-Lynn-Shacham's Short Signature Scheme
  - `example/ecdsa.py`: Implementation of ECDSA Signature Scheme

# Setup

```
> git clone https://github.com/elliptic-shiho/ecpy && cd ecpy
> pip install --upgrade .
```

## Require modules:
* six

If installed following modules, `ecpy` uses that.

* `gmpy2`
* `primefac` or [primefac-fork](https://github.com/elliptic-shiho/primefac-fork/)

# Setup with Native modules
`ecpy` has a native module. If you want to speed-up ecpy, you can use that.

```
> cd ecpy/cpp
> make
> sudo make install
```

The description is here: [cpp/README.md](cpp/README.md)

# Speed
See `bench_pairing.py` .

```
Sat Jan 28 03:55:47 JST 2017 ~/prog/lab/mathematics/ecpy 100%
> time python scripts/bench_pairing.py
[+] found gmpy! use gmpy.is_prime
[+] Found native module!
Point (25 : 5635866983L + 11271733966Lw : 1) on Elliptic Curve y^2 = x^3 + 1 over ExtendedFiniteField(13208557577, "x^2+x+1")
[+] Weil Pairing: 
weil: 14744.50 usec/pass
weil: 14842.80 usec/pass
weil: 14626.05 usec/pass
weil: 14773.50 usec/pass
weil: 14582.25 usec/pass
weil: 14560.95 usec/pass
weil: 14622.40 usec/pass
weil: 15073.19 usec/pass
weil: 14874.60 usec/pass
weil: 14695.20 usec/pass
[+] Tate Pairing: 
tate: 3079.89 usec/pass
tate: 3083.21 usec/pass
tate: 3079.10 usec/pass
tate: 3072.70 usec/pass
tate: 3069.81 usec/pass
tate: 3089.30 usec/pass
tate: 3077.90 usec/pass
tate: 3116.00 usec/pass
tate: 3095.35 usec/pass
tate: 3193.65 usec/pass
================================================================
weil: 14739.54 usec/pass
tate: 3095.69 usec/pass

real    0m4.293s
user    0m4.256s
sys     0m0.012s
```

# Special Thanks
Cybozu Labs Youth (5th, 6th) - http://labs.cybozu.co.jp/youth.html

# License
This software released under MIT License. Please see /LICENSE.

