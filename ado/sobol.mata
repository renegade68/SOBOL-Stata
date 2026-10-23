*! sobol_mata 0.4.2 - 23oct2026
/* Pure Mata, Gray-code order, 32-bit integer numerators stored exactly
   as doubles. Direction data: Joe--Kuo D(6), first 1024 dimensions.
   Algorithm: Joe and Kuo (2008), Notes on generating Sobol sequences.
   Optional LMS+digital shift, local MT19937; never touches Stata's RNG. */
version 16.0
mata:

struct sobol_state {
    real scalar dimension, index, seed, scrambled
    real matrix directions, lms
    real rowvector numerator, xor_table, shift
    string scalar prng
}

void _sobol_integer(real scalar x, lo, hi, string scalar label)
{
    if (missing(x) || x != floor(x) || x < lo || x > hi) {
        errprintf("%s must be an integer in [%18.0f, %18.0f]\n", label, lo, hi)
        _error(3300)
    }
}

/* An 8-bit XOR lookup table, flattened columnwise. No global cache or RNG. */
real rowvector _sobol_xor_table()
{
    real scalar k
    real matrix t, high
    t = 0
    for (k=0; k<8; k++) {
        /* XOR's high bit differs only in the off-diagonal quadrants. */
        high = t :+ 2^k
        t = (t,high \ high,t)
    }
    return(vec(t)')
}

/* Inputs are equally sized rowvectors (scalars included), in [0,2^32). */
real rowvector _sobol_xor32(real rowvector a, b, table)
{
    real scalar k, weight
    real rowvector out, aa, bb
    aa = a
    bb = b
    out = J(1,cols(a),0)
    weight = 1
    for (k=1; k<=4; k++) {
        out = out + weight * table[1,mod(aa,256):+256*mod(bb,256):+1]
        aa = floor(aa:/256)
        bb = floor(bb:/256)
        weight = weight * 256
    }
    return(out)
}

real matrix _sobol_directions_scalar(real scalar d, real rowvector table)
{
    real matrix p, v
    real scalar j, k, r, s, a, x
    v = J(32,d,0)
    v[,1] = 2:^(31::0)
    if (d==1) return(v)
    p = _sobol_direction_data()
    for (j=2; j<=d; j++) {
        s = p[j-1,1]
        a = p[j-1,2]
        for (k=1; k<=s; k++) v[k,j] = p[j-1,k+2] * 2^(32-k)
        for (k=s+1; k<=32; k++) {
            x = _sobol_xor32(v[k-s,j],floor(v[k-s,j]/2^s),table)
            for (r=1; r<s; r++) {
                if (mod(floor(a/2^(s-1-r)),2)) {
                    x = _sobol_xor32(x,v[k-r,j],table)
                }
            }
            v[k,j] = x
        }
    }
    return(v)
}


/* Small dimensions retain v0.4.1 arithmetic to avoid grouping overhead. */
real matrix _sobol_directions(real scalar d, real rowvector table)
{
    real matrix p, v
    real scalar k, r, s
    real rowvector j, a, x, mask
    if (d<=16) return(_sobol_directions_scalar(d,table))
    v = J(32,d,0)
    v[,1] = 2:^(31::0)
    if (d==1) return(v)
    p = _sobol_direction_data()
    p = p[|1,1 \ d-1,cols(p)|]
    /* Equal polynomial degrees share the recurrence's row offsets. Each
       coordinate still uses its own initial values and coefficient mask. */
    for (s=1; s<=max(p[,1]); s++) {
        j = selectindex(p[,1]':==s):+1
        if (!cols(j)) continue
        a = p[j:-1,2]'
        for (k=1; k<=s; k++) v[k,j] = p[j:-1,k+2]' * 2^(32-k)
        for (k=s+1; k<=32; k++) {
            x = _sobol_xor32(v[k-s,j],floor(v[k-s,j]/2^s),table)
            for (r=1; r<s; r++) {
                mask = mod(floor(a/2^(s-1-r)),2)
                x = _sobol_xor32(x,v[k-r,j]:*mask,table)
            }
            v[k,j] = x
        }
    }
    return(v)
}

void _sobol_check_state(struct sobol_state scalar state)
{
    _sobol_integer(state.dimension,1,1024,"state dimension")
    _sobol_integer(state.index,0,2^32,"state index")
    if (rows(state.directions)!=32 || cols(state.directions)!=state.dimension ||
        cols(state.numerator)!=state.dimension || cols(state.shift)!=state.dimension ||
        cols(state.xor_table)!=65536) {
        errprintf("invalid Sobol state; create it with sobol_init()\n")
        _error(3200)
    }
}

/* start is zero based. 2^32 is an exhausted-state sentinel. */
void sobol_seek(struct sobol_state scalar state, real scalar start)
{
    real scalar g, k
    real rowvector x
    _sobol_check_state(state)
    _sobol_integer(start,0,2^32,"start")
    x = state.shift
    if (start < 2^32) {
        g = _sobol_xor32(start,floor(start/2),state.xor_table)
        k = 1
        while (g>0) {
            if (mod(g,2)) x = _sobol_xor32(x,state.directions[k,.],state.xor_table)
            g = floor(g/2)
            k++
        }
    }
    state.numerator = x
    state.index = start
}

/* Validate packed LMS rows: row r contains bits 1..r, with diagonal one.
   Bit 1 is the MOST significant fractional bit (weight 2^31). */
void _sobol_validate_lms(real matrix lms, real rowvector shift)
{
    real scalar r, weight, d
    d=cols(lms)
    _sobol_integer(d,1,1024,"LMS dimension")
    if (rows(lms)!=32 || cols(shift)!=d) {
        errprintf("LMS must be 32 by d and shift must be 1 by d\n")
        _error(3200)
    }
    if (any(lms:>=.) || any(lms:<0) || any(lms:>=2^32) || any(lms:!=floor(lms)) ||
        any(shift:>=.) || any(shift:<0) || any(shift:>=2^32) || any(shift:!=floor(shift))) {
        errprintf("LMS masks and digital shifts must be uint32 integers\n")
        _error(3300)
    }
    for (r=1; r<=32; r++) {
        weight=2^(32-r)
        if (any(mod(lms[r,.],weight):!=0) || any(mod(floor(lms[r,.]/weight),2):!=1)) {
            errprintf("LMS matrices must be lower triangular with unit diagonal (MSB first)\n")
            _error(3300)
        }
    }
}

void _sobol_apply_lms(struct sobol_state scalar state, real matrix lms,
                      real rowvector shift)
{
    real scalar j
    real matrix a, b, v
    real rowvector inverseweights
    real colvector weights
    _sobol_validate_lms(lms,shift)
    weights=2:^(31::0)
    inverseweights=(1:/weights)'
    v=J(32,state.dimension,0)
    for (j=1; j<=state.dimension; j++) {
        a=mod(floor(lms[,j]*inverseweights),2)
        b=mod(floor(state.directions[,j]*inverseweights),2)
        v[,j]=mod(a*b',2)'*weights
    }
    state.directions=v
    state.lms=lms
    state.shift=shift
    state.scrambled=1
}

/* Exactly 32 MT uint32 words per coordinate: shift, then rows 2..32.
   Each row uses the high r-1 random bits plus its fixed diagonal bit. */
void _sobol_seeded_lms(struct sobol_state scalar state, real scalar seed)
{
    real rowvector words, shift
    real matrix lms
    real scalar r, weight
    real rowvector offsets
    words=_sobol_mt_words(seed,32*state.dimension,state.xor_table)
    offsets=32*(0..(state.dimension-1))
    shift=words[1,offsets:+1]
    lms=J(32,state.dimension,0)
    lms[1,.]=J(1,state.dimension,2^31)
    for (r=2; r<=32; r++) {
        weight=2^(33-r)
        lms[r,.]=floor(words[1,offsets:+r]/weight)*weight:+weight/2
    }
    _sobol_apply_lms(state,lms,shift)
    state.seed=seed
    state.prng="MT19937-init_genrand-2002"
}

struct sobol_state scalar sobol_init(real scalar d, | real scalar start, seed)
{
    struct sobol_state scalar state
    if (args()==1) start=0
    _sobol_integer(d,1,1024,"dimension")
    _sobol_integer(start,0,2^32,"start")
    if (args()==3) _sobol_integer(seed,0,2^32-1,"seed")
    state.dimension = d
    state.index = 0
    state.numerator = J(1,d,0)
    state.shift = J(1,d,0)
    state.lms = J(0,d,0)
    state.seed = .
    state.scrambled = 0
    state.prng = "none"
    state.xor_table = _sobol_xor_table()
    state.directions = _sobol_directions(d,state.xor_table)
    if (args()==3) _sobol_seeded_lms(state,seed)
    sobol_seek(state,start)
    return(state)
}

/* Expert API: reproduce an explicit LMS configuration without any PRNG. */
struct sobol_state scalar sobol_init_lms(real matrix lms, real rowvector shift,
                                       | real scalar start)
{
    struct sobol_state scalar state
    if (args()==2) start=0
    _sobol_validate_lms(lms,shift)
    state=sobol_init(cols(lms),start)
    _sobol_apply_lms(state,lms,shift)
    state.prng="external-LMS"
    sobol_seek(state,start)
    return(state)
}

/* Broadcast a uint32 row over a matrix; bound workspace by the tile size. */
real matrix _sobol_xor_rows(real matrix a, real rowvector b, table)
{
    real matrix aa, out, indices
    real rowvector bb
    real scalar k, weight
    aa=a
    bb=b
    out=J(rows(a),cols(a),0)
    weight=1
    for (k=1; k<=4; k++) {
        indices=mod(aa,256):+256*mod(bb,256):+1
        out=out+weight*colshape(table[1,vec(indices)'],rows(a))'
        aa=floor(aa:/256)
        bb=floor(bb:/256)
        weight=weight*256
    }
    return(out)
}

/* Gray codes reflect on doubling. For B=2^k and 0<=r<B,
   g(q*B+r) XOR g(q*B) = g(r). This remains true after any LMS.
   A <=1024-row tile avoids one Mata iteration per point. */
real matrix _sobol_draw_tiled(struct sobol_state scalar state, real scalar n)
{
    struct sobol_state scalar cursor
    real matrix tile, out
    real scalar bits, B, k, first, offset, take, index
    bits=min((10,ceil(ln(n)/ln(2)),floor(ln(262144/state.dimension)/ln(2))))
    B=2^bits
    tile=J(1,state.dimension,0)
    for (k=1; k<=bits; k++) {
        tile=tile \ _sobol_xor_rows(tile[rows(tile)::1,.],
                                          state.directions[k,.],state.xor_table)
    }
    out=J(n,state.dimension,0)
    cursor=state
    first=1
    index=state.index
    while (first<=n) {
        offset=mod(index,B)
        take=min((B-offset,n-first+1))
        sobol_seek(cursor,index-offset)
        out[|first,1 \ first+take-1,state.dimension|]=
            _sobol_xor_rows(tile[|offset+1,1 \ offset+take,state.dimension|],
                            cursor.numerator,state.xor_table)/2^32
        first=first+take
        index=index+take
    }
    /* Commit only after the complete output and next position are available. */
    sobol_seek(cursor,index)
    state.numerator=cursor.numerator
    state.index=index
    return(out)
}

real matrix sobol_draw(struct sobol_state scalar state, real scalar n)
{
    real matrix out
    real rowvector x
    real scalar i, index, q, c
    _sobol_check_state(state)
    _sobol_integer(n,0,2^32,"n")
    if (n > 2^32-state.index) {
        errprintf("Sobol stream exhausted: start + n must not exceed 2^32\n")
        _error(3300)
    }
    if (n>=64) return(_sobol_draw_tiled(state,n))
    out = J(n,state.dimension,0)
    if (!n) return(out)
    x = state.numerator
    index = state.index
    for (i=1; i<=n; i++) {
        out[i,.] = x / 2^32
        index++
        if (index < 2^32) {
            q = index
            c = 1
            while (mod(q,2)==0) {
                c++
                q = q/2
            }
            x = _sobol_xor32(x,state.directions[c,.],state.xor_table)
        }
        else x = J(1,state.dimension,0)
    }
    /* Commit only after successful completion; failed calls do not advance. */
    state.numerator = x
    state.index = index
    return(out)
}

real matrix sobol(real scalar n, d, | real scalar start, seed)
{
    struct sobol_state scalar state
    if (args()==2) start=0
    _sobol_integer(n,0,2^32,"n")
    _sobol_integer(start,0,2^32,"start")
    if (n > 2^32-start) {
        errprintf("start + n must not exceed 2^32\n")
        _error(3300)
    }
    if (args()==4) state = sobol_init(d,start,seed)
    else state = sobol_init(d,start)
    return(sobol_draw(state,n))
}

real matrix sobol_base2(real scalar m, d, | real scalar seed)
{
    _sobol_integer(m,0,32,"m")
    if (args()==3) return(sobol(2^m,d,0,seed))
    return(sobol(2^m,d))
}

/* Called by the ado wrapper, after it has created double variables. */
void _sobol_store(string rowvector names, real scalar start, | real scalar seed)
{
    struct sobol_state scalar state
    real scalar first, last, n
    real rowvector variables
    n = st_nobs()
    variables = st_varindex(names)
    if (args()==3) state = sobol_init(cols(names),start,seed)
    else state = sobol_init(cols(names),start)
    for (first=1; first<=n; first=first+4096) {
        last = min((n,first+4095))
        st_store((first,last),variables,sobol_draw(state,last-first+1))
    }
}

end
