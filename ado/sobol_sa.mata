*! sobol_mata 0.4.2 - 23oct2026
/* Independent uniform inputs; scalar deterministic output; S1 and ST. */
version 16.0
mata:

struct sobol_sa_design {
    real scalar m, n, d, reps, blocksize, level, maxbytes, n_evaluations
    real scalar output_bytes, external_bytes, input_block_bytes
    real colvector seeds
    real matrix bounds
    string scalar package_version, layout
}
struct sobol_sa_batch {
    real matrix keys, x
}
struct sobol_sa_result {
    real rowvector S1, ST, se, outside_unit
    real matrix replicates, covariance, ci
    real colvector scaled_variance, output_scale
    real scalar near_constant
    struct sobol_sa_design scalar design
    string scalar estimator
}

struct sobol_sa_design scalar sobol_sa_plan(real scalar m, real matrix bounds,
    real vector seeds, | real scalar blocksize, level, maxbytes)
{
    struct sobol_sa_design scalar p
    real colvector ordered
    real scalar R, d
    if (args()<4) blocksize=64
    if (args()<5) level=95
    if (args()<6) maxbytes=512*1024^2
    _sobol_integer(m,1,32,"sensitivity m")
    d=rows(bounds)
    _sobol_integer(d,1,512,"sensitivity dimension")
    if (cols(bounds)!=2) _error(3200)
    if (hasmissing(bounds) || any(bounds[,2]:<=bounds[,1]) ||
        hasmissing(bounds[,2]-bounds[,1])) {
        errprintf("bounds must be finite increasing lower/upper pairs\n")
        _error(3300)
    }
    R=length(seeds)
    if (rows(seeds)!=1 && cols(seeds)!=1) _error(3200)
    _sobol_integer(R,2,2^32,"replicate count")
    if (hasmissing(seeds) || any(seeds:<0) || any(seeds:>=2^32) ||
        any(seeds:!=floor(seeds))) _error(3300)
    ordered=sort(vec(seeds),1)
    if (any(ordered[|2\R|]:==ordered[|1\R-1|])) {
        errprintf("replicate seeds must be distinct\n")
        _error(3300)
    }
    _sobol_integer(blocksize,1,2^32,"blocksize")
    _sobol_integer(maxbytes,1,2^53-1,"maxbytes")
    if (missing(level) || level<=0 || level>=100) _error(3300)
    if (2^m*R*(d+2)>2^53-1) _error(3300)
    p.m=m; p.n=2^m; p.d=d; p.reps=R
    p.bounds=bounds; p.seeds=vec(seeds)
    p.blocksize=min((blocksize,p.n)); p.level=level; p.maxbytes=maxbytes
    p.n_evaluations=p.n*R*(d+2)
    p.output_bytes=8*p.n*(d+2)
    /* Conservative workspace estimates, excluding caller-owned data, Mata
       allocator overhead, callback workspace and the loaded generator. */
    p.external_bytes=8*p.n_evaluations+6*p.output_bytes+96*d^2+32*R*(d+1)
    p.input_block_bytes=24*p.blocksize*(d+2)*(d+4)+64*p.blocksize*d
    p.package_version=sobol_version()
    p.layout="sobol-sa-csv-1"
    return(p)
}

void _sobol_sa_check(struct sobol_sa_design scalar p)
{
    struct sobol_sa_design scalar q
    q=sobol_sa_plan(p.m,p.bounds,p.seeds,p.blocksize,p.level,p.maxbytes)
    if (p.n!=q.n || p.d!=q.d || p.reps!=q.reps ||
        p.n_evaluations!=q.n_evaluations || p.output_bytes!=q.output_bytes ||
        p.external_bytes!=q.external_bytes || p.input_block_bytes!=q.input_block_bytes ||
        p.layout!=q.layout || p.package_version!=q.package_version) {
        errprintf("invalid or incompatible sensitivity plan; recreate with sobol_sa_plan()\n")
        _error(3300)
    }
}

void _sobol_sa_memory(struct sobol_sa_design scalar p, real scalar bytes, nrows)
{
    if (bytes>p.maxbytes || nrows>2147483647) {
        errprintf("sensitivity workspace estimate %g bytes exceeds budget or matrix row limit\n",bytes)
        errprintf("inspect the plan; reduce blocksize/N/R or explicitly increase maxbytes\n")
        _error(3300)
    }
}

real matrix _sobol_sa_scale(struct sobol_sa_design scalar p, real matrix U)
{
    return(p.bounds[,1]':+U:*(p.bounds[,2]-p.bounds[,1])')
}

/* Role-major blocks: all A, then all AB_1,...,AB_D, then all B.
   Keys, not row order, define the exchange contract. */
struct sobol_sa_batch scalar _sobol_sa_expand(struct sobol_sa_design scalar p,
    real scalar r, start, real matrix U)
{
    struct sobol_sa_batch scalar b
    real matrix A, B, C
    real scalar k, n, first, last
    n=rows(U)
    A=_sobol_sa_scale(p,U[,1..p.d])
    B=_sobol_sa_scale(p,U[,p.d+1..2*p.d])
    b.x=J(n*(p.d+2),p.d,.)
    b.keys=J(n*(p.d+2),4,.)
    for (k=0;k<=p.d+1;k++) {
        C=A
        if (k>0 && k<=p.d) C[,k]=B[,k]
        if (k==p.d+1) C=B
        first=k*n+1; last=(k+1)*n
        if (n) {
            b.x[|first,1\last,p.d|]=C
            b.keys[|first,1\last,4|]=(J(n,1,r),start:+(0::n-1),J(n,1,k),J(n,1,p.seeds[r]))
        }
    }
    return(b)
}

struct sobol_sa_batch scalar sobol_sa_inputs(struct sobol_sa_design scalar p,
    real scalar r, start, count)
{
    struct sobol_state scalar s
    _sobol_sa_check(p)
    _sobol_integer(r,1,p.reps,"replicate")
    _sobol_integer(start,0,p.n,"base point start")
    _sobol_integer(count,0,p.n-start,"base point count")
    _sobol_sa_memory(p,24*count*(p.d+2)*(p.d+4)+64*count*p.d,count*(p.d+2))
    s=sobol_init(2*p.d,start,p.seeds[r])
    return(_sobol_sa_expand(p,r,start,sobol_draw(s,count)))
}

/* Return S1..., ST..., scaled variance, output scale. All Y are model outputs. */
real rowvector _sobol_sa_estimate(real matrix Y)
{
    real matrix C
    real colvector A, B, AB, pooled
    real rowvector out
    real scalar N, d, j, scale, V, center
    N=rows(Y); d=cols(Y)-2
    if (hasmissing(Y)) _error(3300)
    scale=max(abs(Y))
    if (scale==0) {
        errprintf("sensitivity indices undefined: zero output variance\n")
        _error(3300)
    }
    C=Y:/scale
    center=sum(quadcolsum(C))/(N*(d+2))
    C=C:-center
    A=C[,1]; B=C[,d+2]; pooled=A\B
    V=quadcrossdev(pooled,mean(pooled),pooled,mean(pooled))/(2*N)
    if (missing(V) || V<=0) {
        errprintf("sensitivity indices undefined: zero output variance\n")
        _error(3300)
    }
    out=J(1,2*d+2,.)
    for (j=1;j<=d;j++) {
        AB=C[,j+1]
        out[j]=quadcross(B,AB-A)/(N*V)
        out[d+j]=quadcross(A-AB,A-AB)/(2*N*V)
    }
    out[2*d+1]=V; out[2*d+2]=scale
    if (hasmissing(out)) _error(3300)
    return(out)
}

struct sobol_sa_result scalar _sobol_sa_summary(struct sobol_sa_design scalar p,
    real matrix estimates)
{
    struct sobol_sa_result scalar result
    real rowvector avg
    real scalar critical
    result.replicates=estimates[,1..2*p.d]
    avg=quadcolsum(result.replicates:/p.reps)
    result.S1=avg[1..p.d]; result.ST=avg[p.d+1..2*p.d]
    result.covariance=quadcrossdev(result.replicates,avg,result.replicates,avg)/(p.reps*(p.reps-1))
    result.se=sqrt(diagonal(result.covariance))'
    critical=invttail(p.reps-1,(100-p.level)/200)
    result.ci=(avg-critical*result.se\avg+critical*result.se)
    if (hasmissing(result.covariance) || hasmissing(result.ci)) _error(3300)
    result.scaled_variance=estimates[,2*p.d+1]
    result.output_scale=estimates[,2*p.d+2]
    result.near_constant=any(result.scaled_variance:<1e-12)
    result.outside_unit=(avg:<0):|(avg:>1)
    result.design=p
    result.estimator="saltelli2010-s1/jansen-st; centered; pooled-popvar; mean-of-replicate-ratios"
    return(result)
}

real scalar _sobol_sa_key_ok(struct sobol_sa_design scalar p, real rowvector key)
{
    if (hasmissing(key) || any(key:!=floor(key))) return(0)
    if (key[1]<1 || key[1]>p.reps || key[2]<0 || key[2]>=p.n ||
        key[3]<0 || key[3]>p.d+1) return(0)
    return(key[4]==p.seeds[key[1]])
}

struct sobol_sa_result scalar _sobol_sa_packed(struct sobol_sa_design scalar p,
    real matrix Y)
{
    real matrix estimates
    real scalar r
    estimates=J(p.reps,2*p.d+2,.)
    for (r=1;r<=p.reps;r++) estimates[r,.]=_sobol_sa_estimate(Y[|(r-1)*p.n+1,1\r*p.n,p.d+2|])
    return(_sobol_sa_summary(p,estimates))
}

struct sobol_sa_result scalar sobol_sa_analyze(struct sobol_sa_design scalar p,
    real matrix keys, real colvector values)
{
    real matrix Y
    real scalar i, row, col
    _sobol_sa_check(p)
    if (cols(keys)!=4 || rows(keys)!=p.n_evaluations || rows(values)!=rows(keys) || cols(values)!=1) _error(3200)
    _sobol_sa_memory(p,p.external_bytes,p.n*p.reps)
    if (hasmissing(values)) _error(3300)
    Y=J(p.n*p.reps,p.d+2,.)
    for (i=1;i<=rows(keys);i++) {
        if (!_sobol_sa_key_ok(p,keys[i,.])) _error(3300)
        row=(keys[i,1]-1)*p.n+keys[i,2]+1; col=keys[i,3]+1
        if (!missing(Y[row,col])) {
            errprintf("duplicate sensitivity output key\n")
            _error(3300)
        }
        Y[row,col]=values[i]
    }
    return(_sobol_sa_packed(p,Y))
}

struct sobol_sa_result scalar sobol_sensitivity(pointer(real matrix function) scalar f,
    real scalar m, real matrix bounds, real vector seeds, | real scalar blocksize, level, maxbytes)
{
    struct sobol_sa_design scalar p
    struct sobol_state scalar s
    real matrix Y, U, A, B, C, values, estimates
    real scalar r, first, take, j
    if (args()<5) blocksize=64
    if (args()<6) level=95
    if (args()<7) maxbytes=512*1024^2
    p=sobol_sa_plan(m,bounds,seeds,blocksize,level,maxbytes)
    if (f==NULL) _error(3300)
    _sobol_sa_memory(p,6*p.output_bytes+128*p.blocksize*p.d+96*p.d^2+32*p.reps*(p.d+1),p.n)
    estimates=J(p.reps,2*p.d+2,.)
    for (r=1;r<=p.reps;r++) {
        Y=J(p.n,p.d+2,.)
        s=sobol_init(2*p.d,0,p.seeds[r])
        for (first=1;first<=p.n;first=first+p.blocksize) {
            take=min((p.blocksize,p.n-first+1))
            U=sobol_draw(s,take)
            A=_sobol_sa_scale(p,U[,1..p.d]); B=_sobol_sa_scale(p,U[,p.d+1..2*p.d])
            for (j=0;j<=p.d+1;j++) {
                C=A
                if (j>0 && j<=p.d) C[,j]=B[,j]
                if (j==p.d+1) C=B
                values=(*f)(C)
                if (rows(values)!=take || cols(values)!=1) {
                    errprintf("sensitivity callback must return one scalar per input row\n")
                    _error(3200)
                }
                if (hasmissing(values)) _error(3300)
                Y[|first,j+1\first+take-1,j+1|]=values
            }
        }
        estimates[r,.]=_sobol_sa_estimate(Y)
    }
    return(_sobol_sa_summary(p,estimates))
}

string scalar _sobol_sa_csv(real rowvector v)
{
    return(invtokens(strtrim(strofreal(v,"%24.17g")),","))
}
string colvector _sobol_sa_header(struct sobol_sa_design scalar p)
{
    string colvector h
    real scalar j
    h=J(2+p.d+p.reps,1,"")
    h[1]="#sobol_sa,1,"+p.package_version
    h[2]="#plan,"+_sobol_sa_csv((p.m,p.d,p.reps,p.level))
    for (j=1;j<=p.d;j++) h[2+j]="#bounds,"+_sobol_sa_csv((j,p.bounds[j,.]))
    for (j=1;j<=p.reps;j++) h[2+p.d+j]="#seed,"+_sobol_sa_csv((j,p.seeds[j]))
    return(h)
}
void _sobol_sa_iofail(real scalar fh, string scalar message)
{
    fclose(fh)
    errprintf("%s\n",message)
    _error(3300)
}
void _sobol_sa_put(real scalar fh, string scalar line)
{
    if (_fput(fh,line)<0) _sobol_sa_iofail(fh,"sensitivity exchange write failed")
}
string scalar _sobol_sa_get(real scalar fh)
{
    string matrix line
    line=_fget(fh)
    if (rows(line)==0 || fstatus(fh)<0) _sobol_sa_iofail(fh,"incomplete sensitivity exchange file")
    return(line)
}

/* New file only. A partial file has no valid #end trailer and is rejected. */
void sobol_sa_export(struct sobol_sa_design scalar p, string scalar path)
{
    struct sobol_state scalar s
    struct sobol_sa_batch scalar b
    string colvector h
    string scalar columns
    real scalar fh, j, r, start, count, i
    _sobol_sa_check(p)
    _sobol_sa_memory(p,p.input_block_bytes,p.blocksize*(p.d+2))
    if (fileexists(path)) {
        errprintf("exchange destination already exists; use a new filename\n")
        _error(3300)
    }
    fh=fopen(path,"w")
    h=_sobol_sa_header(p)
    for (j=1;j<=rows(h);j++) _sobol_sa_put(fh,h[j])
    columns="replicate,index,role,seed"
    for (j=1;j<=p.d;j++) columns=columns+",x"+strtrim(strofreal(j,"%9.0f"))
    _sobol_sa_put(fh,columns)
    for (r=1;r<=p.reps;r++) {
        s=sobol_init(2*p.d,0,p.seeds[r])
        for (start=0;start<p.n;start=start+p.blocksize) {
            count=min((p.blocksize,p.n-start))
            b=_sobol_sa_expand(p,r,start,sobol_draw(s,count))
            for (i=1;i<=rows(b.x);i++) _sobol_sa_put(fh,_sobol_sa_csv((b.keys[i,.],b.x[i,.])))
        }
    }
    _sobol_sa_put(fh,"#end,"+_sobol_sa_csv(p.n_evaluations))
    fclose(fh)
}

real rowvector _sobol_sa_parse(string scalar line)
{
    real rowvector v
    real scalar j, comma
    string scalar field
    v=J(1,5,.)
    for (j=1;j<=5;j++) {
        comma=strpos(line,",")
        if ((j<5 && comma==0) || (j==5 && comma!=0)) return(J(1,5,.))
        if (j<5) {
            field=substr(line,1,comma-1)
            line=substr(line,comma+1,.)
        }
        else field=line
        v[j]=strtoreal(strtrim(field))
    }
    return(v)
}

struct sobol_sa_result scalar sobol_sa_import(struct sobol_sa_design scalar p,
    string scalar path)
{
    real matrix Y
    real rowvector v
    real scalar fh, j, row, col
    string colvector h
    string matrix line
    _sobol_sa_check(p)
    _sobol_sa_memory(p,p.external_bytes,p.n*p.reps)
    Y=J(p.n*p.reps,p.d+2,.)
    h=_sobol_sa_header(p)
    fh=fopen(path,"r")
    for (j=1;j<=rows(h);j++) {
        if (_sobol_sa_get(fh)!=h[j]) _sobol_sa_iofail(fh,"sensitivity design metadata mismatch")
    }
    if (_sobol_sa_get(fh)!="replicate,index,role,seed,y") _sobol_sa_iofail(fh,"invalid response column header")
    for (j=1;j<=p.n_evaluations;j++) {
        v=_sobol_sa_parse(_sobol_sa_get(fh))
        if (hasmissing(v)) _sobol_sa_iofail(fh,"invalid or missing response value")
        if (!_sobol_sa_key_ok(p,v[1..4])) _sobol_sa_iofail(fh,"response key or seed mismatch")
        row=(v[1]-1)*p.n+v[2]+1; col=v[3]+1
        if (!missing(Y[row,col])) _sobol_sa_iofail(fh,"duplicate response key")
        Y[row,col]=v[5]
    }
    if (_sobol_sa_get(fh)!="#end,"+_sobol_sa_csv(p.n_evaluations)) _sobol_sa_iofail(fh,"invalid response trailer/count")
    line=_fget(fh)
    if (rows(line)!=0 || fstatus(fh)!=-1) _sobol_sa_iofail(fh,"extra data or read error after response trailer")
    fclose(fh)
    return(_sobol_sa_packed(p,Y))
}
end
