{smcl}
{* *! sobol_mata 0.4.2 23oct2026}{...}
{title:Title}

{pstd}{bf:sobol_integrate()} {hline 2} Replicated scrambled Sobol integration

{title:Setup and syntax}

{phang2}{cmd:sobol_setup}

{p 8 12 2}{it:struct sobol_integral scalar}
{cmd:sobol_integrate(}{cmd:&}{it:f}{cmd:(),}{it:m}{cmd:,}{it:d}{cmd:,}{it:seeds}
[{cmd:,}{it:blocksize}[{cmd:,}{it:level}]]{cmd:)}

{title:Description}

{pstd}Integrates a user-defined function over [0,1]^d using R separate
LMS+shift scrambles, each with N=2^m points beginning at index zero. The function
f(U) must accept a real matrix of points and return a finite real matrix with
one row per point and K>=1 columns, one per integrand. Each row's evaluation
must be independent of batching; K must remain fixed across all calls.
Write an explicit wrapper for transformations or additional model parameters.
No points are dropped and no coordinate transformations are applied automatically.

{pstd}m is an integer from 0 through 32; d is from 1 through 1,024.
seeds is a row or column vector of R>=2 distinct uint32 integers, including zero.
Duplicate seeds are rejected because repeated identical scrambles do not provide
additional replications. Distinct pseudo-random seeds are a reproducible design
convention, not proof of mathematical independence for arbitrary seed choices.

{pstd}blocksize defaults to 4,096, must be an integer from 1 through 2^32,
and is capped at N. level defaults to 95 and must be strictly between 0 and 100.
These are logical limits; large N*R still requires corresponding computation.

{title:Returned structure}

{pstd}For R replicates and K integrands:
{p2colset 5 26 28 2}{...}
{p2col:estimate}1 by K average of the replicate integral estimates{p_end}
{p2col:replicates}R by K estimates, in the supplied seed order{p_end}
{p2col:covariance}K by K estimated covariance of estimate (replicate sample covariance / R){p_end}
{p2col:se}1 by K square roots of the covariance diagonal{p_end}
{p2col:ci}2 by K marginal approximate t intervals: lower row, upper row; df=R-1{p_end}
{p2col:seeds}R by 1 input seeds, preserving order{p_end}
{p2col:n, dimension, reps}N, d, R{p_end}
{p2col:blocksize, level}Effective block size and confidence level{p_end}
{p2col:package_version}Generator package version{p_end}
{p2colreset}{...}

{pstd}The intervals are marginal, not simultaneous. They describe estimated
numerical integration error, not parameter sampling uncertainty. A finite
32-bit grid can introduce error unaccounted for by between-scramble variation.
Coverage is not guaranteed for all integrands. A zero standard error for a
constant integrand correctly produces a degenerate interval; a zero estimated
error in general does not prove exact integration.

{title:Example}

{pstd}Run this block in a fresh Stata session after adding the repository ado
folder to the search path. It requires no dataset or external program.
{phang2}{cmd:sobol_setup}
{phang2}{cmd:mata:}
{phang2}{cmd:real matrix f_demo(real matrix u) return((exp(rowsum(u)), rowsum(u):<1))}
{phang2}{cmd:fit = sobol_integrate(&f_demo(),12,3,20261000:+(1..32))}
{phang2}{cmd:fit.estimate}
{phang2}{cmd:fit.se}
{phang2}{cmd:fit.ci}
{phang2}{cmd:end}

{pstd}The exact integrals are (exp(1)-1)^3 and 1/6. Finite-sample estimates vary.

{pstd}Within a compiled Mata function, declare the result before assignment:
{phang2}{cmd:struct sobol_integral scalar fit}

{title:Numerical behavior and errors}

{pstd}Each block uses quad-precision summation after division by N; block
contributions use compensated addition. Different block sizes preserve points
but may differ slightly due to floating-point summation. The full N by d matrix
is not stored. Memory includes point and callback-output blocks, the R by K
replicate matrix, and the K by K covariance matrix.

{pstd}Missing, overflowed, wrongly shaped, or changing-column callback outputs
are errors. Callback errors propagate to the caller. The generator and this
integrator do not modify Stata's RNG or dataset; arbitrary side effects of a
user-provided callback remain the user's responsibility. Use a deterministic,
side-effect-free rowwise callback for reproducible integration.

{pstd}For endpoint handling see {help sobol}; for seed rules see {help mf_sobol}.
The seed protocol matches v0.2 and does not promise the same
points as standard SciPy Sobol given only the same numeric seed.
