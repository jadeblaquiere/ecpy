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

import ecpy.curves as curves
import abc
import sys

if sys.version_info < (3,):
    integer_types = (int, long,)
else:
    integer_types = (int,)


class _abstractstaticmethod(staticmethod):
    __slots__ = ()

    def __init__(self, function):
        super(_abstractstaticmethod, self).__init__(function)
        function.__isabstractmethod__ = True
    __isabstractmethod__ = True


_curve = curves.curve_secp256k1
# _curve = curves.curve_secp384r1
# _curve = curves.curve_bauer9

_generator_LUT_bits = 8
_generator_list = []


def _modinv(a, m):
    lastr, r, x, lastx = a, m, 0, 1
    while r:
        lastr, (q, r) = r, divmod(lastr, r)
        x, lastx = lastx - q*x, x
    return lastx % m


"""ecpy.point: Native python implementation of elliptic curve point math"""


class PointBase (object):
    """Abstract base clase for elliptic curve points in finite fields"""
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def compress(self):
        """Return a string representing the compressed form of the point"""
        pass

    @_abstractstaticmethod
    def decompress(textrep):
        """Construct a point from a string representing the compressed form"""
        pass

    @abc.abstractmethod
    def __eq__(self, q):
        pass

    @abc.abstractmethod
    def __ne__(self, q):
        pass

    @abc.abstractmethod
    def __add__(self, q):
        """Operator for adding one point to another or itself (doubling). """
        """Points must be from the same curve"""
        pass

    @abc.abstractmethod
    def __mul__(self, n):
        """Operator for multiplication of Point by scalar value n"""
        pass

    @abc.abstractmethod
    def __rmul__(self, n):
        pass

    @abc.abstractmethod
    def affine(self):
        """Return a tuple (x,y) of the the affine coordinates of the point"""
        pass

    @abc.abstractmethod
    def __str__(self):
        pass

    @abc.abstractmethod
    def __repr__(self):
        pass


class Point (PointBase):
    """Point class for short Weierstrass curve algebra. Short Weierstrass
    format curves are defined by the equation y**2 = x**3 + a*x + b.

    Point is derived from abstract base class PointBase.
    """
    p = _curve['p']
    n = _curve['n']
    a = _curve['a']
    b = _curve['b']
    bits = _curve['bits']

    def __init__(self, x=None, y=None, z=None, infinity=False):
        self.x = x
        self.y = y
        self.z = z if z is not None else 1
        self.is_infinite = True if x is None else infinity
        self.dcache = None

    @classmethod
    def set_curve(cls, c):
        """Sets curve parameters. Takes a dictionary as input with keys:
            p - defines prime field Fp. Coordinate values are modulo p
            G - generator point for cyclic subgroup of points on cruve
            n - order of G, also size of the subgroup
            a, b - curve equation coefficients
            bits - bitsize of prime field, i.e. log2(p)
        """
        cls.p = c['p']
        cls.n = c['n']
        cls.a = c['a']
        cls.b = c['b']
        cls.bits = c['bits']

    def _from_jacobian(self):
        if not self.is_infinite:
            if self.z != 1:
                zinv = _modinv(self.z, Point.p)
                self.x = (self.x * zinv ** 2) % Point.p
                self.y = (self.y * zinv ** 3) % Point.p
                self.z = 1

    def compress(self):
        """Return a string representing the compressed form of the point"""
        if self.is_infinite:
            return b'infinity'
        P = self.affine()
        pfmt = b'%%0%dx' % (int((Point.bits + 7) / 8) * 2)
        return (b'03' if (P[1] % 2) else b'02') + (pfmt % P[0])

    @staticmethod
    def decompress(textrep):
        """Construct a point from a string representing the compressed form"""
        if textrep == b'infinity':
            return Point(infinity=True)
        P = [0, 0]
        P[0] = x = int(textrep[2:], 16)
        sign = int(textrep[:2], 16) & 1
        beta = pow(int(x * x * x + Point.a * x + Point.b),
                   int((Point.p + 1) // 4), Point.p)
        P[1] = (Point.p - beta) if ((beta + sign) & 1) else beta
        return Point(P[0], P[1])

    def _copy(self):
        return Point(self.x, self.y, self.z, self.is_infinite)

    def _double(self):
        if self.is_infinite or self.y == 0:
            return Point(infinity=True)
        ysq = (self.y * self.y) % Point.p
        S = (4 * self.x * ysq) % Point.p
        M = (3 * self.x * self.x + Point.a * self.z ** 4) % Point.p
        nx = (M * M - 2 * S) % self.p
        ny = (M * (S - nx) - 8 * ysq * ysq) % Point.p
        nz = (2 * self.y * self.z) % Point.p
        return Point(nx, ny, nz)

    def __eq__(self, q):
        assert isinstance(q, Point)
        if self.is_infinite:
            if q.is_infinite:
                return True
            else:
                return False
        U1 = (self.x * q.z * q.z) % Point.p
        U2 = (q.x * self.z * self.z) % Point.p
        S1 = (self.y * q.z ** 3) % Point.p
        S2 = (q.y * self.z ** 3) % Point.p
        if U1 == U2:
            if S1 == S2:
                return True
        return False

    def __ne__(self, q):
        return not (self == q)

    def __add__(self, q):
        """Operator for adding one point to another or itself (doubling).
        Points must be from the same curve"""
        assert isinstance(q, Point)
        if self.is_infinite:
            return q._copy()
        if q.is_infinite:
            return self._copy()
        U1 = (self.x * q.z * q.z) % Point.p
        U2 = (q.x * self.z * self.z) % Point.p
        S1 = (self.y * q.z ** 3) % Point.p
        S2 = (q.y * self.z ** 3) % Point.p
        if U1 == U2:
            if S1 != S2:
                return Point(infinity=True)
            return self._double()
        H = U2 - U1
        R = S2 - S1
        H2 = (H * H) % Point.p
        H3 = (H * H2) % Point.p
        U1H2 = (U1 * H2) % Point.p
        nx = (R * R - H3 - 2 * U1H2) % Point.p
        ny = (R * (U1H2 - nx) - S1 * H3) % Point.p
        nz = (H * self.z * q.z) % Point.p
        return Point(nx, ny, nz)

    def __mul__(self, n):
        """Operator for multiplication of Point by scalar value n"""
        assert isinstance(n, integer_types)
        if self.is_infinite or n == 0:
            return Point(infinity=True)
        nmod = n % Point.n
        if self.dcache is None:
            d = self._copy()
            self.dcache = []
            accum = Point(infinity=True)
            cn = self.bits
            while cn > 0:
                dval = (d.x, d.y, d.z, d.is_infinite)
                self.dcache.append(dval)
                if (nmod % 2) == 1:
                    accum += (d)
                if cn == 1:
                    return accum
                d = d._double()
                nmod = nmod // 2
                cn -= 1
        else:
            accum = Point(infinity=True)
            idx = 0
            while nmod > 0:
                if (nmod % 2) == 1:
                    dval = self.dcache[idx]
                    accum += Point(dval[0], dval[1], dval[2], dval[3])
                if nmod == 1:
                    return accum
                nmod = nmod // 2
                idx += 1

    def __rmul__(self, n):
        return self.__mul__(n)

    def affine(self):
        """Return a tuple (x,y) of the the affine coordinates of the point"""
        self._from_jacobian()
        if self.is_infinite:
            return b'infinity'
        return (self.x, self.y)

    def __str__(self):
        return self.compress()

    def __repr__(self):
        self._from_jacobian()
        return "Point(" + str(self.x) + ", " + str(self.y) + ")"


class Generator (Point):
    """Optimized subclass for elliptic curve generator points. Look up tables
    are generated in constructor which accererate scalar multiplication.
    """

    def __init__(self, x, y, z=None):
        super(self.__class__, self).__init__(x, y, z)
        self._recalculate_tables()

    @classmethod
    def set_curve(cls, c):
        Point.set_curve(c)

    @staticmethod
    def init(x, y, z=None):
        """Will find a matching generator or construct one if required.
        Creating generators is computationally expensive so this method is
        preferred to calling normal constructor."""
        P = Point(x, y, z)
        for gen in _generator_list:
            if P == gen:
                return gen
        gen = Generator(x, y, z)
        _generator_list.append(gen)
        return gen

    def _recalculate_tables(self):
        self.ntables = ((Generator.bits - 1) // _generator_LUT_bits) + 1
        self.tablesize = 2 ** _generator_LUT_bits
        self.maskbits = self.tablesize - 1
        super(self.__class__, self).__mul__(1)
        bitbase = 0
        self.lut = []
        while bitbase < Generator.bits:
            nlut = []
            for i in range(self.tablesize):
                accum = Point(infinity=True)
                for j in range(_generator_LUT_bits):
                    if bitbase+j < Generator.bits:
                        if (i & (0x01 << j)) != 0:
                            dval = self.dcache[bitbase + j]
                            accum = accum + Point(
                                dval[0], dval[1], dval[2], dval[3])
                nlut.append(accum)
            accum._from_jacobian()
            self.lut.append(nlut)
            bitbase += _generator_LUT_bits

    def __mul__(self, n):
        nshift = n % Generator.n
        iidx = 0
        accum = Point(infinity=True)
        while nshift != 0:
            idx = nshift & self.maskbits
            accum = accum + self.lut[iidx][idx]
            iidx += 1
            nshift >>= _generator_LUT_bits
        return accum

    def __rmul__(self, n):
        return self.__mul__(n)
