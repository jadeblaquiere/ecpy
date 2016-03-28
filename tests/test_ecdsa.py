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

from ecpy import Point, Generator
from ecpy import curves
from ecpy import ECDSA
from Crypto.Random import random
from Crypto.Cipher import AES
from Crypto.Util import Counter
import hashlib
import binascii
import base64

_curve = curves.curve_secp112r1
#_curve = curves.curve_secp256k1
#_curve = curves.curve_secp384r1

Point.set_curve(_curve)
_G = Generator(_curve['G'][0], _curve['G'][1])
ECDSA.set_curve(_curve)
ECDSA.set_generator(_G)

ecdsa = ECDSA()
recdsa = ECDSA()

friends = ['Alice', 'Bob', 'Carl', 'Donna', 'Eve']

data = []
for f in friends:
    datum = {}
    datum['name'] = f
    p = random.randint(1,_curve['n']-1)
    P = p * _G
    datum['privkey'] = p
    datum['pubkey'] = P
    data.append(datum)

print data

for i in range(100000):
    for d1 in data:
        for d2 in data:
            if d1 == d2:
                print 'skipping %s to %s' % (d1['name'], d2['name'])
            else:
                print 'sending %s to %s' % (d1['name'], d2['name'])
                msglen = random.randint(1,1023)
                msg=''
                for j in range(msglen):
                    if random.randint(0,7) == 0:
                        msg += ' '
                    else:
                        msg += chr(ord('a') + random.randint(0,25))
                adlen = random.randint(1,255)
                ad=''
                for j in range(adlen):
                    if random.randint(0,7) == 0:
                        ad += ' '
                    else:
                        ad += chr(ord('a') + random.randint(0,25))
                send = {}
                send['ad'] = ad
                ecdhkey = d2['pubkey'] * d1['privkey']
                bs = AES.block_size
                iv = random.randint(0, (1 << (8*bs))-1)
                fmt = '%%0%dx' % (bs * 2)
                ivhex = fmt % iv
                ivbin = binascii.unhexlify(ivhex)
                counter = Counter.new(AES.block_size * 8, initial_value=iv)
                key = hashlib.sha256(ecdhkey.compress()).digest()
                cryptor = AES.new(key, AES.MODE_CTR, counter = counter)
                ciphertext = ivbin + cryptor.encrypt(msg)
                b64cipher = base64.b64encode(ciphertext)
                send['b64cipher'] = b64cipher
                sig = ecdsa.sign(d1['privkey'], ciphertext, ad)
                send['sig'] = sig
                print send
                print
                print
                recv = send
                recdhkey = d1['pubkey'] * d2['privkey']
                assert recdhkey == ecdhkey
                rciphertext = base64.b64decode(recv['b64cipher'])
                assert rciphertext == ciphertext
                rve = recdsa.verify(d1['pubkey'], recv['sig'], rciphertext, recv['ad'])
                assert rve == True
                rbs = AES.block_size
                assert rbs == bs
                rivbin = rciphertext[:rbs]
                assert rivbin == ivbin
                riv = int(binascii.hexlify(rivbin),16)
                rcounter = Counter.new(AES.block_size * 8, initial_value=riv)
                rkey = hashlib.sha256(recdhkey.compress()).digest()
                assert rkey == key
                rcryptor = AES.new(rkey, AES.MODE_CTR, counter=rcounter)
                plaintext = rcryptor.decrypt(rciphertext[rbs:])
                assert plaintext == msg
            
            
