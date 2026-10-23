*! sobol_mata 0.4.2 - 23oct2026
/* Local MT19937, init_genrand (2002), uint32 output, no global RNG calls.
   Adapted from the algorithm by Makoto Matsumoto and Takuji Nishimura.

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
version 16.0
mata:

/* Low 32 bits of a product; all intermediate integers remain below 2^53. */
real scalar _sobol_mul32(real scalar a, b)
{
    real scalar a0, a1, b0, b1
    a0=mod(a,65536); a1=floor(a/65536)
    b0=mod(b,65536); b1=floor(b/65536)
    return(mod(a0*b0 + 65536*mod(a0*b1+a1*b0,65536),2^32))
}

real rowvector _sobol_mt_words(real scalar seed, n, real rowvector table)
{
    real rowvector mt, out, indices, mixed, y, z
    real scalar i, first, last, take, word, prior
    _sobol_integer(seed,0,2^32-1,"seed")
    _sobol_integer(n,0,32768,"number of scramble words")
    out=J(1,n,0)
    if (!n) return(out)
    mt=J(1,624,0)
    mt[1]=seed
    for (i=2; i<=624; i++) {
        prior=mt[i-1]
        /* prior>>30 is at most 3: only the low byte changes under XOR. */
        word=256*floor(prior/256)+table[1,mod(prior,256)+256*floor(prior/2^30)+1]
        mt[i]=mod(_sobol_mul32(1812433253,word)+i-1,2^32)
    }
    first=1
    while (first<=n) {
        /* In-place dependencies permit batches 1:227, 228:454, 455:623.
           The last word uses the already updated first word separately. */
        for (i=1;i<=455;i=i+227) {
            last=min((i+226,623))
            indices=i..last
            mixed=floor(mt[1,indices]/2^31)*2^31+mod(mt[1,indices:+1],2^31)
            y=_sobol_xor32(mt[1,mod(indices:-1:+397,624):+1],floor(mixed/2),table)
            mt[1,indices]=_sobol_xor32(y,mod(mixed,2)*2567483615,table)
        }
        word=floor(mt[624]/2^31)*2^31+mod(mt[1],2^31)
        y=_sobol_xor32(mt[397],floor(word/2),table)
        mt[624]=_sobol_xor32(y,mod(word,2)*2567483615,table)
        take=min((624,n-first+1))
        out[1,first..(first+take-1)]=mt[1,1..take]
        first=first+take
    }
    /* Temper the collected words together; all intermediates are exact doubles. */
    y=_sobol_xor32(out,floor(out/2^11),table)
    z=mod(y*2^7,2^32)
    z=((z:+2636928640)-_sobol_xor32(z,J(1,n,2636928640),table))/2
    y=_sobol_xor32(y,z,table)
    z=mod(y*2^15,2^32)
    z=((z:+4022730752)-_sobol_xor32(z,J(1,n,4022730752),table))/2
    y=_sobol_xor32(y,z,table)
    return(_sobol_xor32(y,floor(y/2^18),table))
}

end
