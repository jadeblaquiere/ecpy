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

import curves
import abc


class _abstractstaticmethod(staticmethod):
    __slots__ = ()
    def __init__(self, function):
        super(_abstractstaticmethod, self).__init__(function)
        function.__isabstractmethod__ = True
    __isabstractmethod__ = True


_curve = curves.curve_secp256k1
# _curve = curves.curve_secp384r1
# _curve = curves.curve_bauer9
_ed_curve = curves.ed_curve_curve41417

_generator_LUT_bits = 8


def _modinv(a, m):
    lastr, r, x, lastx = a, m, 0, 1
    while r:
        lastr, (q, r) = r, divmod(lastr, r)
        x, lastx = lastx - q*x, x
    return lastx % m

"""ecpy.point: Native python implementation of elliptic curve point math"""

class PointBase (object):
    """Abstract base clase for elliptic curve points in finite integer fields"""
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
        # todo: validate point is on curve
        self.x = x
        self.y = y
        self.z = z if z is not None else 1
        self.is_infinite = True if x is None else infinity
        self.dcache = None

    @classmethod
    def set_curve(cls, curve):
        """Sets curve parameters. Takes a dictionary as input with keys:
            p - defines prime field Fp. Coordinate values are modulo p
            G - generator point for cyclic subgroup of points on cruve
            n - order of G, also size of the subgroup
            a, b - curve equation coefficients
            bits - bitsize of prime field, i.e. log2(p)
        """
        cls.p = curve['p']
        cls.n = curve['n']
        cls.a = curve['a']
        cls.b = curve['b']
        cls.bits = curve['bits']

    def _from_jacobian(self):
        if not self.is_infinite:
            if self.z != 1:
                zinv = _modinv(self.z, Point.p)
                self.x = (self.x * zinv ** 2) % Point.p
                self.y = (self.y * zinv ** 3) % Point.p
                self.z = 1

    def compress(self):
        """Return a string representing the compressed form of the point"""
        P = self.affine()
        pfmt = '%%0%dx' % (int((Point.bits + 7) / 8) * 2)
        return ('03' if (P[1] % 2) else '02') + (pfmt % P[0])

    @staticmethod
    def decompress(textrep):
        """Construct a point from a string representing the compressed form"""
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
        assert isinstance(q, Point)
        return not (self == q)

    def __add__(self, q):
        """Operator for adding one point to another or itself (doubling). """
        """Points must be from the same curve"""
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
        assert isinstance(n, (int, long))
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
                    accum += d
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
            return 'infinity'
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

    def __init__(self, x, y, z=None, infinity=False):
        super(self.__class__, self).__init__(x, y)
        self._recalculate_tables()

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


class EdPoint (PointBase):
    """Point class for Edwards curve algebra. Edwards format curves are
    defined by the equation X**2 + y**2 = c**2 ( 1 + d * x**2 * y**2 ). 
    
    Point is derived from abstract base class PointBase.
    
    """
    p = _ed_curve['p']
    l = _ed_curve['l']
    h = _ed_curve['h']
    c = _ed_curve['c']
    d = _ed_curve['d']
    bits = _ed_curve['bits']

    def __init__(self, x=None, y=None, z=None, t=None, infinity=False):
        # todo: validate point is on curve
        self.dcache = None
        if infinity:
            self.is_infinte = True
            self.x, self.y, self.t, self.z = None, None, None. None
            return
        self.x = x
        self.y = y
        if t is None:
            assert self.x is not None
            assert self.y is not None
            self.t = (self.x * self.y) * EdPoint.p
        else: 
            self.t = t
        self.z = z if z is not None else 1

    @classmethod
    def set_curve(cls, curve):
        """Sets curve parameters. Takes a dictionary as input with keys:
            p - defines prime field Fp. Coordinate values are modulo p
            G - generator point for cyclic subgroup of points on cruve
            l - l*h = order of G, also size of the subgroup, l prime
            h - cofactor of order
            a, b - curve equation coefficients
            bits - bitsize of prime field, i.e. log2(p)
        """
        cls.p = curve['p']
        cls.l = curve['l']
        cls.h = curve['h']
        cls.c = curve['c']
        cls.d = curve['d']
        cls.ccd2 = (2 * cls.c * cls.c * cls.d) % cls.p
        cls.bits = curve['bits']

    def _from_extended_projective(self):
        if not self.is_infinite:
            if self.z != 1:
                zinv = _modinv(self.z, EdPoint.p)
                self.x = (self.x * zinv) % EdPoint.p
                self.y = (self.y * zinv) % EdPoint.p
                self.t = (self.x * self.y) % EdPoint.p
                self.z = 1

    def compress(self):
        """Return a string representing the compressed form of the point"""
        # todo: actually implement compression
        # per https://eprint.iacr.org/2015/673.pdf
        P = self.affine()
        pfmt = '%%0%dx' % (int((EdPoint.bits + 7) / 8) * 2)
        return '04' + (pfmt % P[0]) + (pfmt %P[1])

    @staticmethod
    def decompress(textrep):
        """Construct a point from a string representing the compressed form"""
        psz = int((EdPoint.bits + 7) / 8) * 2
        x = int(textrep[2:2+psz], 16)
        y = int(textrep[2+psz:], 16)
        return EdPoint(P[0], P[1])

    def _copy(self):
        return Point(self.x, self.y, self.z, self.t, self.is_infinite)

    def _double(self):
        if self.is_infinite:
            return Point(infinity=True)
        # https://hyperelliptic.org/EFD/g1p/auto-edwards-projective.html
        A = (self.x + self.y) % EdPoint.p
        C = (self.x * self.x) % EdPoint.p
        D = (self.y * self.y) % EdPoint.p
        E = (C + D) % EdPoint.p
        if EdPoint.c == 1:
            H = (self.z * self.z) % EdPoint.p
            J = (E - (2 * H)) % EdPoint.p
            nx = ((B - E) * J) % EdPoint.p
            ny = (E * (C - D)) % EdPoint.p
            nz = (E * J) % EdPoint.p
        else:
            H = (EdPoint.c * self.z * self.z) % EdPoint.p
            J = (E - (2 * H)) % EdPoint.p
            nx = (EdPoint.c * (B - E) * J) % EdPoint.p
            ny = (EdPoint.c * E * (C - D)) % EdPoint.p
            nz = (E * J) % EdPoint.p
        return EdPoint(nx, ny, nz)

    def __eq__(self, q):
        assert isinstance(q, EdPoint)
        if self.is_infinite:
            if q.is_infinite:
                return True
            else:
                return False
        U1 = (self.x * q.z) % EdPoint.p
        U2 = (q.x * self.z) % EdPoint.p
        S1 = (self.y * q.z) % EdPoint.p
        S2 = (q.y * self.z) % EdPoint.p
        if U1 == U2:
            if S1 == S2:
                return True
        return False

    def __ne__(self, q):
        assert isinstance(q, EdPoint)
        return not (self == q)

    def __add__(self, q):
        """Operator for adding one point to another or itself (doubling). """
        """Points must be from the same curve"""
        assert isinstance(q, EdPoint)
        if self.is_infinite:
            return q._copy()
        if q.is_infinite:
            return self._copy()
        A = (self.z * q.z) % EdPoint.p
        B = (A * A) % EdPoint.p
        C = (self.x * q.x) % EdPoint.p
        C = (self.y * q.y) % EdPoint.p
        E = (EdPoint.d * C * D) % EdPoint.p
        F = (B - E) % EdPoint.p
        G = (B + E) % EdPoint.p
        tmp = (self.x + self.y) * (q.x + q.y) % EdPoint.p
        nx = (A * F * (tmp - C - D)) % EdPoint.p
        ny = (A * G * (D - C)) % EdPoint.p
        nz = (EdPoint.c * F * G) % EdPoint.p
        return EdPoint(nx, ny, nz)

    def __mul__(self, n):
        """Operator for multiplication of Point by scalar value n"""
        assert isinstance(n, (int, long))
        if self.is_infinite or n == 0:
            return EdPoint(infinity=True)
        nmod = n % EdPoint.n
        if self.dcache is None:
            d = self._copy()
            self.dcache = []
            accum = EdPoint(infinity=True)
            cn = self.bits
            while cn > 0:
                dval = (d.x, d.y, d.z, d.t, d.is_infinite)
                self.dcache.append(dval)
                if (nmod % 2) == 1:
                    accum += d
                if cn == 1:
                    return accum
                d = d._double()
                nmod = nmod // 2
                cn -= 1
        else:
            accum = EdPoint(infinity=True)
            idx = 0
            while nmod > 0:
                if (nmod % 2) == 1:
                    dval = self.dcache[idx]
                    accum += EdPoint(dval[0], dval[1], dval[2], dval[3], dval[4])
                if nmod == 1:
                    return accum
                nmod = nmod // 2
                idx += 1

    def __rmul__(self, n):
        return self.__mul__(n)

    def affine(self):
        """Return a tuple (x,y) of the the affine coordinates of the point"""
        self._from_extended_projective()
        if self.is_infinite:
            return 'infinity'
        return (self.x, self.y)

    def __str__(self):
        return self.compress()

    def __repr__(self):
        self._from_jacobian()
        return "EdPoint(" + str(self.x) + ", " + str(self.y) + ")"


