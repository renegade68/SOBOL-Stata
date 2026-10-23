{smcl}
{* *! sobol_mata 0.4.2 23oct2026}{...}
{title:Title}

{pstd}{bf:sobol()} {hline 2} Pure Mata Sobol matrices and stateful streams

{title:Setup}

{pstd}First add the complete repository ado folder to Stata's search path as
shown in {help sobol}. Keep the .mata implementation files with the .ado files.
{phang2}{cmd:sobol_setup}

{pstd}sobol_setup loads the functions and the sobol_state structure into Mata.
It is idempotent and can be called again after mata clear. It does not clear
other functions, matrices, the dataset, or the random-number state.
Start a new Stata session if an earlier version is already loaded;
sobol_setup rejects mixing versions.

{title:Syntax}

{p 8 12 2}{it:real matrix} {cmd:sobol(}{it:n}{cmd:,}{it:d}[{cmd:,}{it:start}[{cmd:,}{it:seed}]]{cmd:)}
{p 8 12 2}{it:real matrix} {cmd:sobol_base2(}{it:m}{cmd:,}{it:d}[{cmd:,}{it:seed}]{cmd:)}
{p 8 12 2}{it:struct sobol_state scalar} {cmd:sobol_init(}{it:d}[{cmd:,}{it:start}[{cmd:,}{it:seed}]]{cmd:)}
{p 8 12 2}{it:struct sobol_state scalar} {cmd:sobol_init_lms(}{it:lms}{cmd:,}{it:shift}[{cmd:,}{it:start}]{cmd:)}
{p 8 12 2}{it:real matrix} {cmd:sobol_draw(}{it:state}{cmd:,}{it:n}{cmd:)}
{p 8 12 2}{it:void} {cmd:sobol_seek(}{it:state}{cmd:,}{it:start}{cmd:)}
{p 8 12 2}{it:string scalar} {cmd:sobol_version()}

{title:Arguments and results}

{pstd}d is an integer from 1 through 1,024. All returned matrices are n by d
and contain double-precision values in [0,1). n is a nonnegative integer.
start is a zero-based integer index; its default is zero.
start+n cannot exceed 2^32. Missing and noninteger arguments are errors.
n=0 returns a 0 by d matrix, including at the exhausted position 2^32.

{pstd}An optional seed enables LMS plus digital shift. It must be an integer
from 0 through 2^32-1; an explicitly supplied missing seed is an error.
Omitting it gives the original unscrambled sequence. Supply start=0 when
passing a seed to sobol() or sobol_init(). No global RNG is read or changed.

{pstd}sobol_base2(m,d) returns the prefix with n=2^m points, starting at zero;
m must be an integer from zero through 32. These are logical bounds, not a
promise that a matrix of this size fits in memory. The returned matrix alone
needs 8*n*d bytes.

{pstd}sobol_init creates a state. sobol_draw returns the next n rows and advances
the state. Rejected draws do not advance it. sobol_seek changes its next index
directly, without iterating through skipped points, and returns nothing.
The state at index 2^32 is exhausted and accepts only an empty draw or a seek.
The generator never silently wraps around.

{pstd}state.index is the next zero-based index. Configuration fields
state.scrambled, state.seed, state.prng, state.lms, and state.shift can be
read to record a scramble. Do not modify structure fields directly.
The other fields are internal implementation details.
Ordinary structure assignment creates an independent copy of a stream.
Keep the state between blocks to avoid rebuilding direction numbers and the
8-bit XOR table. Names beginning with _sobol are internal.

{title:Examples}

{phang2}{cmd:mata: U = sobol_base2(10,3)}
{phang2}{cmd:mata: V = sobol(1024,3,1024)}

{pstd}Stateful generation at the interactive Mata prompt:
{phang2}{cmd:mata: s = sobol_init(3)}
{phang2}{cmd:mata: A = sobol_draw(s,512)}
{phang2}{cmd:mata: B = sobol_draw(s,512)}
{phang2}{cmd:mata: assert((A \ B)==sobol_base2(10,3))}
{phang2}{cmd:mata: sobol_seek(s,0)}

{pstd}Scrambled matrices and a stream with the same explicit seed:
{phang2}{cmd:mata: R = sobol_base2(10,3,20261023)}
{phang2}{cmd:mata: T = sobol(1024,3,1024,20261023)}
{phang2}{cmd:mata: s = sobol_init(3,0,20261023)}
{phang2}{cmd:mata: A = sobol_draw(s,512)}
{phang2}{cmd:mata: B = sobol_draw(s,512)}
{phang2}{cmd:mata: assert((A \ B)==R)}

{pstd}Inside a compiled Mata function, explicitly declare the state:
{phang2}{cmd:struct sobol_state scalar s}

{title:Algorithm and compatibility}

{pstd}The direction set is new-joe-kuo-6.21201, revision 16 September 2010,
restricted to its first 1,024 dimensions. The implementation computes 32
direction integers per coordinate and applies the Gray-code recurrence.
The numerators are integers below 2^32 and are exactly representable as doubles;
division by 2^32 is exact. XOR is implemented using an eight-bit lookup table.
No integer-64 emulation or floating-point bit casts are used.

{pstd}Version 0.3 batches draws of at least 64 points using reflected Gray codes
and matrix XOR operations. Temporary tiles have at most 1,024 rows and 262,144
elements. Smaller requests retain the scalar recurrence. Values, ordering, and
seed consumption are unchanged. The full requested output matrix is still
allocated; retain a stream and draw blocks to limit memory use.

{pstd}With a seed, a private MT19937 using init_genrand (2002) generates 32
uint32 words per dimension: the first is the digital shift; words 2 through 32
provide the high r-1 entries of LMS row r, whose diagonal entry is one.
Row 1 is fixed. Bits are ordered from most significant to least significant.
These unit-diagonal lower-triangular matrices multiply fractional-bit column
vectors over GF(2), then the shift is applied by XOR. Generation and seek use
the same fixed configuration; they consume no further random bits. Increasing
the dimension preserves the coordinate prefix for a given seed.

{pstd}The reference is SciPy Sobol with bits=32 and the same direction set,
ordering, and LMS matrices and shifts, without postprocessing or optimization.
The same numeric seed in standard SciPy Sobol is insufficient because its
random-bit generation and consumption differ. This is not a claim of bitwise
compatibility with an unidentified historical LPTAU implementation.

{title:Explicit LMS configurations}

{pstd}sobol_init_lms(lms,shift,start) accepts a 32 by d matrix of packed
unsigned 32-bit row masks and a 1 by d vector of unsigned 32-bit shifts.
In row r, weight 2^(32-r) is the unit diagonal; lower weights must be zero;
higher weights are the r-1 free entries. For example the identity LMS has
rows 2^31, 2^30, ..., 1 in every column. It returns a state with
state.seed missing, state.prng equal to external-LMS, and state.scrambled=1.
All integers are stored exactly in Mata doubles. Invalid shapes, fractional
values, missing values, or non-unit lower-triangular matrices are rejected.

{pstd}For example, to reconstruct a recorded configuration:
{phang2}{cmd:mata: t = sobol_init_lms(s.lms,s.shift,0)}
{phang2}{cmd:mata: assert(sobol_draw(t,1024)==R)}

{pstd}For each dimension, let w_r be its rth uint32 MT output. The shift is w_1;
the first packed LMS row is 2^31. For r=2,...,32 the packed row mask is
floor(w_r/2^(33-r))*2^(33-r)+2^(32-r). Exactly 32 words are allocated per
dimension, in coordinate order. No further random words are consumed by draw
or seek. Use stateful generation to amortize initialization cost.

{pstd}See {help sobol} for balance, endpoint, and double-storage considerations.
For replicated integration of a user-defined function, see {help mf_sobol_integrate}.
