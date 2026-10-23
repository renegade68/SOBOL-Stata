*! sobol_mata 0.4.2 - 23oct2026
/* Replicated randomized integration on [0,1]^d; no automatic transforms. */
version 16.0
mata:

struct sobol_integral {
    real scalar n, dimension, reps, level, blocksize
    real colvector seeds
    real matrix replicates, covariance, ci
    real rowvector estimate, se
    string scalar package_version
}

/* f(U) must return a finite real matrix, one row per input point and one
   column per integrand. It must evaluate rows independently of batch size. */
struct sobol_integral scalar sobol_integrate(pointer(real matrix function) scalar f,
    real scalar m, d, real vector seeds, | real scalar blocksize, level)
{
    struct sobol_integral scalar result
    struct sobol_state scalar state
    real colvector ordered
    real matrix values
    real rowvector total, correction, delta, updated
    real scalar R, r, N, first, take, K, critical
    if (args()<5) blocksize=4096
    if (args()<6) level=95
    if (f==NULL) {
        errprintf("integrand must be a function pointer\n")
        _error(3300)
    }
    _sobol_integer(m,0,32,"m")
    _sobol_integer(d,1,1024,"dimension")
    _sobol_integer(blocksize,1,2^32,"blocksize")
    if (missing(level) || level<=0 || level>=100) {
        errprintf("confidence level must lie strictly between 0 and 100\n")
        _error(3300)
    }
    R=length(seeds)
    _sobol_integer(R,2,2^32,"number of replicate seeds")
    if (any(seeds:>=.) || any(seeds:<0) || any(seeds:>=2^32) ||
        any(seeds:!=floor(seeds))) {
        errprintf("replicate seeds must be uint32 integers\n")
        _error(3300)
    }
    ordered=sort(vec(seeds),1)
    if (any(ordered[|2 \ R|]:==ordered[|1 \ R-1|])) {
        errprintf("replicate seeds must be distinct; duplicate scrambles cannot estimate error\n")
        _error(3300)
    }
    N=2^m
    blocksize=min((blocksize,N))
    K=0
    for (r=1; r<=R; r++) {
        state=sobol_init(d,0,seeds[r])
        for (first=1; first<=N; first=first+blocksize) {
            take=min((blocksize,N-first+1))
            values=(*f)(sobol_draw(state,take))
            if (rows(values)!=take || cols(values)<1) {
                errprintf("integrand must return one row per point and at least one column\n")
                _error(3200)
            }
            if (hasmissing(values)) {
                errprintf("integrand returned missing/nonfinite values in replicate %g\n",r)
                _error(3300)
            }
            if (K==0) {
                K=cols(values)
                result.replicates=J(R,K,0)
            }
            if (cols(values)!=K) {
                errprintf("integrand changed its number of columns between calls\n")
                _error(3200)
            }
            if (first==1) {
                total=J(1,K,0)
                correction=J(1,K,0)
            }
            /* Scale before quad summation, then compensate across blocks. */
            delta=quadcolsum(values:/N)-correction
            updated=total+delta
            correction=(updated-total)-delta
            total=updated
        }
        result.replicates[r,.]=total
    }
    result.estimate=quadcolsum(result.replicates:/R)
    result.covariance=quadcrossdev(result.replicates,result.estimate,
                                  result.replicates,result.estimate)/(R*(R-1))
    result.se=sqrt(diagonal(result.covariance))'
    critical=invttail(R-1,(100-level)/200)
    result.ci=(result.estimate-critical*result.se \ result.estimate+critical*result.se)
    if (hasmissing(result.replicates) || hasmissing(result.covariance) ||
        hasmissing(result.ci)) {
        errprintf("integration summary overflowed; rescale the integrand\n")
        _error(3300)
    }
    result.n=N
    result.dimension=d
    result.reps=R
    result.level=level
    result.blocksize=blocksize
    result.seeds=vec(seeds)
    result.package_version=sobol_version()
    return(result)
}

/* Load-completion marker: define after all runtime sources compile. */
string scalar sobol_version()
{
    return("0.4.2")
}
end
