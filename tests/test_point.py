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

from ecpy.point import Point, Generator
import ecpy.curves as curves
import time
from Crypto.Random import random

#_curve = curves.curve_secp112r1
_curve = curves.curve_secp256k1
#_curve = curves.curve_secp384r1
#_curve = curves.curve_bauer9
P = _curve['p']
N = _curve['n']
A = _curve['a']

Generator.set_curve(_curve)

def modinv(a, m): 
    lastremainder, remainder, x, lastx = a, m, 0, 1
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(lastremainder, remainder)
        x, lastx = lastx - quotient*x, x
    return lastx % m

def inv(a, n):
    if a == 0:
        return 0
    lm, hm = 1, 0
    low, high = a % n, n
    while low > 1:
        r = high // low
        nm, new = hm - lm * r, high - low * r
        lm, low, hm, high = nm, new, lm, low
    return lm % n

#def inv(a,n):
#    return pow(a,n-2,n)



def isinf(p):
    return p[0] == 0 and p[1] == 0


def to_jacobian(p):
    o = (p[0], p[1], 1)
    return o


def jacobian_double(p):
    if not p[1]:
        return (0, 0, 0)
    ysq = (p[1] ** 2) % P
    S = (4 * p[0] * ysq) % P
    M = (3 * p[0] ** 2 + A * p[2] ** 4) % P
    nx = (M**2 - 2 * S) % P
    ny = (M * (S - nx) - 8 * ysq ** 2) % P
    nz = (2 * p[1] * p[2]) % P
    return (nx, ny, nz)


def jacobian_add(p, q):
    if not p[1] and not p[0]:
        return q
    if not q[1] and not q[0]:
        return p
    #print 'p', str(p)
    #print 'q', str(q)
    U1 = (p[0] * q[2] ** 2) % P
    U2 = (q[0] * p[2] ** 2) % P
    S1 = (p[1] * q[2] ** 3) % P
    S2 = (q[1] * p[2] ** 3) % P
    #print 't', U1, U2, S1, S2
    if U1 == U2:
        if S1 != S2:
            return (0, 0, 1)
        return jacobian_double(p)
    H = U2 - U1
    R = S2 - S1
    H2 = (H * H) % P
    H3 = (H * H2) % P
    U1H2 = (U1 * H2) % P
    #print 'b', H, R, H2, H3, U1H2
    nx = (R ** 2 - H3 - 2 * U1H2) % P
    ny = (R * (U1H2 - nx) - S1 * H3) % P
    nz = (H * p[2] * q[2]) % P
    return (nx, ny, nz)


def from_jacobian(p):
    z = inv(p[2], P)
    return ((p[0] * z**2) % P, (p[1] * z**3) % P)


def jacobian_multiply(a, n):
    if a[1] == 0 or n == 0:
        return (0, 0, 1)
    if n == 1:
        return a
    if n < 0 or n >= N:
        return jacobian_multiply(a, n % N)
    if (n % 2) == 0:
        return jacobian_double(jacobian_multiply(a, n//2))
    if (n % 2) == 1:
        return jacobian_add(jacobian_double(jacobian_multiply(a, n//2)), a)


def fast_multiply(a, n):
    return from_jacobian(jacobian_multiply(to_jacobian(a), n))


def fast_add(a, b):
    return from_jacobian(jacobian_add(to_jacobian(a), to_jacobian(b)))

if __name__ == '__main__':
    p = _curve['p']
    if False:
        for x in range(min(512,p)):
            xi = inv(x, p)
            test = (x*xi) % p
            #print x, xi, x*xi, test
            if x != 0:
                assert test == 1
    
    print('Setting up Generators')
    G = _curve['G']
    print('ref Gen ' + str(G))
    Gpt = Point(G[0],G[1])
    print('Point Gen ' + str(Gpt))
    GenG = Generator.init(G[0],G[1])
    print('Generator Gen ' + str(GenG))
    GenG2 = Generator.init(G[0],G[1])
    print('Generator Gen ' + str(GenG2))
    assert GenG2 is GenG
    assert Gpt.is_valid()
    assert GenG.is_valid()
    print('Generator uncompressed ' + str(GenG.uncompressed_format()))
    InfP = Point()
    ex = InfP.compress()
    exraw = InfP.uncompressed_format()
    InfP2 = Point.decompress(ex)
    InfP3 = Point.decompress(exraw)
    InfP4 = Point.decompress(ex.decode())
    assert ex == exraw
    assert InfP == InfP2
    assert InfP == InfP3
    assert InfP == InfP4
    print('Point @ Infinity')

    if True:
        limit = min(32, _curve['n']) // 2 
        for ax in range(1-limit,limit+1):
            print('ax=', ax)
            x = ax % _curve['n']
            for ay in range(1-limit,limit+1):
                y = ay % _curve['n']
                X = fast_multiply(G, x)
                Xpt = Point(X[0], X[1]) if x != 0 else Point(infinity=True)
                assert Xpt.is_valid()
                Xa = Xpt.affine()
                if x != 0:
                    Xinv = Point(Xa[0],(Xa[1] + 1) % Point.p)
                    assert Xinv.is_valid() != True
                Y = fast_multiply(G, y)
                Ypt = Point(Y[0], Y[1]) if y != 0 else Point(infinity=True)
                assert Ypt.is_valid()
                Ya = Ypt.affine()
                if y != 0:
                    Yinv = Point(Ya[0],(Ya[1] + 1) % Point.p)
                    assert Yinv.is_valid() != True
                 #print 'fd'
                XX = fast_add(X, X)
                XXpt = Xpt._double()
                XXa = XXpt.affine()
                #print 'fa'
                XY = fast_add(X, Y)
                XYpt = Xpt + Ypt
                XYa = XYpt.affine()
                #print str(x), str(y), str(Xa), str(GXa), str(X), str(Ya), str(GYa), str(Y), str(XX), str(XXa), str(GXXa), str(XY), str(XYa), str(GXYa)
                if XXpt.is_infinite:
                    assert XX[0] == 0 and XX[1] == 0
                else:
                    assert XXa[0] == XX[0] and XXa[1] == XX[1]
                if XYpt.is_infinite:
                    assert XY[0] == 0 and XY[1] == 0
                else:
                    assert XYa[0] == XY[0] and XYa[1] == XY[1]
                xy = (x+y) % _curve['n']
                Fxy = fast_multiply(G,xy)
                Pxy = Gpt * xy
                Pxya = Pxy.affine()
                if not Pxy.is_infinite:
                    Pxyc = Pxy.compress()
                    Pxycd = Point.decompress(Pxyc)
                    Pxycds = Point.decompress(Pxyc.decode())
                    Pxyuc = Pxy.uncompressed_format()
                    Pxyucd = Point.decompress(Pxyuc)
                if Pxy.is_infinite:
                    assert Fxy[0] == 0 and Fxy[1] == 0
                    assert XYpt.is_infinite
                else:
                    assert Pxya[0] == Fxy[0] and Pxya[1] == Fxy[1]
                    assert Pxya[0] == XYa[0] and Pxya[1] == XYa[1]
                    assert Pxy == Pxycd
                    assert Pxy == Pxycds
                    assert Pxy == Pxyucd
                    assert Pxyc != Pxyuc
                #print str(xy) + ':', _curve['n'], Fxy, Pxya, Pxyc
                #print Fxy, '%x' % Fxy[0], '%x' % Fxy[1], Pxya, Pxyc
                #print Pxyc
                #print
                if not Xpt.is_infinite and not Ypt.is_infinite:
                    Pecxy = Xpt * y
                    Pecyx = Ypt * x
                    assert Pecxy == Pecyx
    
    print('validating curve parameter is observed')
    G256 = Generator.init(curves.curve_secp256k1['G'][0], curves.curve_secp256k1['G'][1], curve=curves.curve_secp256k1)
    G384 = Generator.init(curves.curve_secp384r1['G'][0], curves.curve_secp384r1['G'][1], curve=curves.curve_secp384r1)
    A256 = 3 * G256
    B256 = 5 * G256
    C256 = 8 * G256
    D256 = 2 * B256
    C256a = C256.affine()
    E256 = Point(C256a[0], C256a[1], curve=curves.curve_secp256k1)
    assert A256 + B256 == C256
    assert A256 + D256 == B256 + E256
    A384 = 3 * G384
    B384 = 5 * G384
    C384 = 8 * G384
    D384 = 2 * B384
    C384a = C384.affine()
    E384 = Point(C384a[0], C384a[1], curve=curves.curve_secp384r1)
    assert A384 + B384 == C384
    assert A384 + D384 == B384 + E384
    fail = False
    try:
        F = A256 + A384
    except AssertionError:
        fail = True
    assert fail == True
    fail = False
    try:
        F = A256 == A384
    except AssertionError:
        fail = True
    assert fail == True

    # generate small curve - not byte aligned, bits != 0 (mod 8)
    curve = {'b': 0, 'bits': 11, 'G': (1364, 917), 'a': 1, 'n': 31, 'p': 1487}
    Gsmall = Generator.init(curve['G'][0], curve['G'][1], curve=curve)

    print("benchmarking :")
    testset = []
    mpztestset = []
    for x in range(250):
        r = random.randint(0, _curve['n'] - 1)
        testset.append(r)
    
    tval = []
    for x in testset:
        tval.append(fast_multiply(G,x))
    brpd_time = time.time()
    for X in tval:
        for y in range(100):
            Pdoub = fast_add(X, X)
    arpd_time = time.time()
    
    tval = []
    for x in testset:
        Gx = Gpt * x
        #Gx._from_jacobian()
        tval.append(Gx)
    bpd_time = time.time()
    for X in tval:
        for y in range(100):
            Gdoub = X._double()
            #Gda = Gdoub.affine()
    apd_time = time.time()
    
    ref_time = arpd_time - brpd_time
    orig_time = apd_time - bpd_time
    
    print('doubling:')
    print('ref =', ref_time)
    print('ecp =', orig_time)
    
    testset = []
    mpztestset = []
    for x in range(250):
        r = random.randint(0, _curve['n'] - 1)
        s = random.randint(0, _curve['n'] - 1)
        testset.append((r, s))
    
    tval = []
    for x in testset:
        tval.append((fast_multiply(G,x[0]), fast_multiply(G,x[1])))
    brpd_time = time.time()
    for X in tval:
        for y in range(100):
            Padd = fast_add(X[0], X[1])
    arpd_time = time.time()
    
    tval = []
    for x in testset:
        Gx0 = Gpt * x[0]
        #Gx0._from_jacobian()
        Gx1 = Gpt * x[1]
        #Gx1._from_jacobian()
        tval.append((Gx0, Gx1))
    bpd_time = time.time()
    for X in tval:
        for y in range(100):
            Gadd = X[0] + X[1]
            #Gaa = Gadd.affine()
    apd_time = time.time()
    
    ref_time = arpd_time - brpd_time
    orig_time = apd_time - bpd_time
    
    print('addition')
    print('ref =', ref_time)
    print('ecp =', orig_time)

    testset = []
    mpztestset = []
    for x in range(25):
        r = random.randint(0, _curve['n'] - 1)
        s = random.randint(0, _curve['n'] - 1)
        testset.append((r, s))
    
    tval = []
    brpd_time = time.time()
    for x in testset:
        X = fast_multiply(G,x[0])
        for y in range(100):
            Pmul = fast_multiply(X, x[1])
    arpd_time = time.time()
    
    tval = []
    bpd_time = time.time()
    for x in testset:
        X = Gpt * x[0]
        #X._from_jacobian()
        for y in range(100):
            Gmul = X * x[1]
            #Gmula = Gmul.affine()
    apd_time = time.time()
    
    ref_time = arpd_time - brpd_time
    orig_time = apd_time - bpd_time
    
    print('scalar mulitplication:')
    print('ref =', ref_time)
    print('ecp =', orig_time)

    testset = []
    mpztestset = []
    for x in range(25):
        r = random.randint(0, _curve['n'] - 1)
        testset.append(r)
    
    tval = []
    brpd_time = time.time()
    for x in testset:
        for y in range(100):
            X = fast_multiply(G,x)
    arpd_time = time.time()
    
    tval = []
    bpd_time = time.time()
    for x in testset:
        for y in range(100):
            X = Gpt * x
    apd_time = time.time()
    
    tval = []
    bgpdb_time = time.time()
    for x in testset:
        for y in range(100):
            X = GenG * x
    agpdb_time = time.time()
    
    ref_time = arpd_time - brpd_time
    orig_time = apd_time - bpd_time
    gnew_time = agpdb_time - bgpdb_time

    print('generation:')
    print('ref =', ref_time)
    print('ecp =', orig_time)
    print('gen =', gnew_time)
