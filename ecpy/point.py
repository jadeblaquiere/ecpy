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

    def __init__(self, x=None, y=None, z=None, infinity=False, curve=None):
        if curve is None:
            self.p = Point.p
            self.n = Point.n
            self.a = Point.a
            self.b = Point.b
            self.bits = Point.bits
        else:
            self.p = curve['p']
            self.n = curve['n']
            self.a = curve['a']
            self.b = curve['b']
            self.bits = curve['bits']
        self.x = x % self.p if x is not None else None
        self.y = y % self.p if y is not None else None
        self.z = z % self.p if z is not None else 1
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
            ------------------------------------------------------------
            The implementation presumes a Short Weierstrass format curves
            y^2 = x^3 + ax + b. The internal calculations are performed
            in Jacobian (X = x/z^2, Y=y/z^3). no restriction is placed
            on the values of a, b
        """
        cls.p = c['p']
        cls.n = c['n']
        cls.a = c['a']
        cls.b = c['b']
        cls.bits = c['bits']

    def curve(self):
        """Returns curve parameters as a dictionary"""
        return {'p': self.p,
                'n': self.n,
                'a': self.a,
                'b': self.b,
                'bits': self.bits}

    def _from_jacobian(self):
        if not self.is_infinite:
            if self.z != 1:
                zinv = _modinv(self.z, self.p)
                self.x = (self.x * zinv ** 2) % self.p
                self.y = (self.y * zinv ** 3) % self.p
                self.z = 1

    def uncompressed_format(self):
        """Return a string of the uncompressed form (x and y) of the point"""
        if self.is_infinite:
            return b'infinity'
        P = self.affine()
        pfmt = '%%0%dx' % (int((self.bits + 7) / 8) * 2)
        return (b'04') + (pfmt % P[0]).encode() + (pfmt % P[1]).encode()

    def compress(self):
        """Return a string of the compressed form of the point"""
        if self.is_infinite:
            return b'infinity'
        P = self.affine()
        pfmt = '%%0%dx' % (int((self.bits + 7) / 8) * 2)
        return (b'03' if (P[1] % 2) else b'02') + (pfmt % P[0]).encode()

    @staticmethod
    def decompress(textrep, curve=None):
        """Construct a point from a string representing the compressed form"""
        if isinstance(textrep, str):
            textrep = textrep.encode()
        if curve is None:
            curve = {'p': Point.p,
                     'n': Point.n,
                     'a': Point.a,
                     'b': Point.b,
                     'bits': Point.bits}
        if textrep == b'infinity':
            return Point(infinity=True, curve=curve)
        P = [0, 0]
        sign = int(textrep[:2], 16) & 1
        if int(textrep[:2], 16) & 4 != 0:
            bytelen = (len(textrep) - 2) // 2
            P[0] = x = int(textrep[2:2+bytelen], 16)
            P[1] = y = int(textrep[2+bytelen:], 16)
        else:
            P[0] = x = int(textrep[2:], 16)
            beta = pow(int(x * x * x + curve['a'] * x + curve['b']),
                       int((curve['p'] + 1) // 4), curve['p'])
            P[1] = (curve['p'] - beta) if ((beta + sign) & 1) else beta
        return Point(P[0], P[1], curve=curve)

    def is_valid(self):
        """Validate the the Point is on the curve"""
        if self.is_infinite:
            return True
        self._from_jacobian()
        ysq = (self.y * self.y) % self.p
        xcu = (self.x * self.x * self.x) % self.p
        ax = (self.a * self.x) % self.p
        right = (xcu + ax + self.b) % self.p
        if ysq == right:
            return True
        return False

    def _copy(self):
        return Point(self.x, self.y, self.z, self.is_infinite,
                     curve=self.curve())

    def _double(self):
        # double uses Bernstein-Lange 2007 formula
        # https://hyperelliptic.org/EFD/g1p/data/shortw/jacobian/doubling/dbl-2007-bl
        if self.is_infinite or self.y == 0:
            return Point(infinity=True, curve=self.curve())
        ysq = (self.y * self.y) % self.p
        S = (4 * self.x * ysq) % self.p
        M = (3 * self.x * self.x + self.a * self.z ** 4) % self.p
        nx = (M * M - 2 * S) % self.p
        ny = (M * (S - nx) - 8 * ysq * ysq) % self.p
        nz = (2 * self.y * self.z) % self.p
        return Point(nx, ny, nz, curve=self.curve())

    def _same_curve(self, q):
        assert isinstance(q, Point)
        if self.p != q.p:
            return False
        if self.n != q.n:
            return False
        if self.a != q.a:
            return False
        if self.b != q.b:
            return False
        if self.bits != q.bits:
            return False
        return True

    def __eq__(self, q):
        assert isinstance(q, Point)
        assert self._same_curve(q) is True
        if self.is_infinite:
            if q.is_infinite:
                return True
            else:
                return False
        U1 = (self.x * q.z * q.z) % self.p
        U2 = (q.x * self.z * self.z) % self.p
        S1 = (self.y * q.z ** 3) % self.p
        S2 = (q.y * self.z ** 3) % self.p
        if U1 == U2:
            if S1 == S2:
                return True
        return False

    def __ne__(self, q):
        return not (self == q)

    def __add__(self, q):
        """Operator for adding one point to another or itself (doubling).
        Points must be from the same curve"""
        # add uses Cohen, Miyaji, Ono formula
        # https://hyperelliptic.org/EFD/g1p/data/shortw/jacobian/addition/add-1998-cmo-2
        assert isinstance(q, Point)
        assert self._same_curve(q) is True
        if self.is_infinite:
            return q._copy()
        if q.is_infinite:
            return self._copy()
        U1 = (self.x * q.z * q.z) % self.p
        U2 = (q.x * self.z * self.z) % self.p
        S1 = (self.y * q.z ** 3) % self.p
        S2 = (q.y * self.z ** 3) % self.p
        if U1 == U2:
            if S1 != S2:
                return Point(infinity=True, curve=self.curve())
            return self._double()
        H = U2 - U1
        R = S2 - S1
        H2 = (H * H) % self.p
        H3 = (H * H2) % self.p
        U1H2 = (U1 * H2) % self.p
        nx = (R * R - H3 - 2 * U1H2) % self.p
        ny = (R * (U1H2 - nx) - S1 * H3) % self.p
        nz = (H * self.z * q.z) % self.p
        return Point(nx, ny, nz, curve=self.curve())

    def __mul__(self, n):
        """Operator for multiplication of Point by scalar value n"""
        assert isinstance(n, integer_types)
        if self.is_infinite or n == 0:
            return Point(infinity=True, curve=self.curve())
        nmod = n % self.n
        if self.dcache is None:
            d = self._copy()
            self.dcache = []
            accum = Point(infinity=True, curve=self.curve())
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
            accum = Point(infinity=True, curve=self.curve())
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
        return self.compress().decode()

    def __repr__(self):
        self._from_jacobian()
        return "self.decompress(" + str(self.x) + ", " + str(self.y) + ")"


class Generator (Point):
    """Optimized subclass for elliptic curve generator points. Look up tables
    are generated in constructor which accererate scalar multiplication.
    """

    def __init__(self, x, y, z=None, curve=None):
        super(self.__class__, self).__init__(x, y, z, curve=curve)
        self._recalculate_tables()

    @classmethod
    def set_curve(cls, c):
        Point.set_curve(c)

    @staticmethod
    def init(x, y, z=None, curve=None):
        """Will find a matching generator or construct one if required.
        Creating generators is computationally expensive so this method is
        preferred to calling normal constructor."""
        P = Point(x, y, z, curve)
        for gen in _generator_list:
            if P._same_curve(gen):
                if P == gen:
                    return gen
        gen = Generator(x, y, z, curve)
        assert gen.is_valid()
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
                accum = Point(infinity=True, curve=self.curve())
                for j in range(_generator_LUT_bits):
                    if bitbase+j < Generator.bits:
                        if (i & (0x01 << j)) != 0:
                            dval = self.dcache[bitbase + j]
                            accum = accum + Point(dval[0], dval[1], dval[2],
                                                  dval[3], self.curve())
                nlut.append(accum)
            accum._from_jacobian()
            self.lut.append(nlut)
            bitbase += _generator_LUT_bits

    def __mul__(self, n):
        nshift = n % Generator.n
        iidx = 0
        accum = Point(infinity=True, curve=self.curve())
        while nshift != 0:
            idx = nshift & self.maskbits
            accum = accum + self.lut[iidx][idx]
            iidx += 1
            nshift >>= _generator_LUT_bits
        return accum

    def __rmul__(self, n):
        return self.__mul__(n)
