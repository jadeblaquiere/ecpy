# Copyright (c) 2016, Joseph deBlaquiere <jadeblaquiere@yahoo.com>
# All rights reserved
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of ecpy nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from ecpy.point import Point
import ecpy.curves as curves
from Crypto.Random import random
import hashlib

_curve = curves.curve_secp256k1

_default_hashalg = hashlib.sha256


def _modinv(a, m):
    lastr, r, x, lastx = a, m, 0, 1
    while r:
        lastr, (q, r) = r, divmod(lastr, r)
        x, lastx = lastx - q*x, x
    return lastx % m


class ECDSA (object):
    p = _curve['p']
    n = _curve['n']
    a = _curve['a']
    bits = _curve['bits']
    G = Point(_curve['G'][0], _curve['G'][1])

    def __init__(self, hashalg=_default_hashalg):
        self.hashalg = hashalg

    @classmethod
    def set_curve(cls, c):
        cls.p = c['p']
        cls.n = c['n']
        cls.a = c['a']
        cls.bits = c['bits']

    @classmethod
    def set_generator(cls, G):
        cls.G = G

    def sign(self, privkey, ciphertext, additional=None):
        hstr = self.hashalg(ciphertext)
        if additional is not None:
            hstr.update(additional)
        e = int(hstr.hexdigest(), 16)
        s = 0
        while s == 0:
            r = 0
            k = 0
            while r == 0:
                k = random.randint(1, ECDSA.n - 1)
                R = ECDSA.G * k
                r = R.affine()[0]
            kinv = _modinv(k, ECDSA.n)
            s = ((e + privkey * r) * kinv) % ECDSA.n
        return (r, s)

    def verify(self, pubkey, signature, ciphertext, additional=None):
        r = signature[0]
        s = signature[1]
        if r < 1 or r > ECDSA.n:
            return False
        if s < 1 or s > ECDSA.n:
            return False
        hstr = self.hashalg(ciphertext)
        if additional is not None:
            hstr.update(additional)
        e = int(hstr.hexdigest(), 16)
        w = _modinv(s, ECDSA.n)
        u1 = (e * w) % ECDSA.n
        u2 = (r * w) % ECDSA.n
        P = (u1 * ECDSA.G) + (u2 * pubkey)
        return (P.affine()[0] % ECDSA.n) == (r % ECDSA.n)
