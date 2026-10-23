{smcl}
{* *! sobol_mata 0.4.2 23oct2026}{...}
{title:Title}

{pstd}{bf:sobol_sensitivity()} {hline 2} First-order and total-effect Sobol indices

{title:Setup and scope}

{phang2}{cmd:sobol_setup}
{pstd}For independent uniform inputs, a deterministic scalar model, and S1/ST.
The dimension D is 1 through 512: a 2D-dimensional Sobol net supplies A and B.
Bounds is a D by 2 matrix of finite increasing lower/upper pairs. No S2,
correlated-input decomposition, stochastic-model noise decomposition, or
built-in nonuniform distributions are provided. This is experimental software.

{title:Syntax}

{phang2}{cmd:p = sobol_sa_plan(m,bounds,seeds[,blocksize[,level[,maxbytes]]])}
{phang2}{cmd:b = sobol_sa_inputs(p,replicate,start,count)}
{phang2}{cmd:fit = sobol_sa_analyze(p,keys,values)}
{phang2}{cmd:fit = sobol_sensitivity(&f(),m,bounds,seeds[,blocksize[,level[,maxbytes]]])}
{phang2}{cmd:sobol_sa_export(p,request_path)}
{phang2}{cmd:fit = sobol_sa_import(p,response_path)}

{pstd}Inside compiled functions, declare p as {cmd:struct sobol_sa_design scalar},
b as {cmd:struct sobol_sa_batch scalar}, and fit as {cmd:struct sobol_sa_result scalar}.
Plan fields are read-only by convention. Recreate a plan to change its settings.

{title:Plan and input batches}

{pstd}m is an integer 1..32 and N=2^m base points are used per replicate.
These are logical bounds, not a promise of feasible allocation. seeds is a
row or column vector of at least two distinct uint32 integers. The default
blocksize is 64 base points, level is 95, and maxbytes is 536870912 (512 MiB).
level is strictly between 0 and 100. maxbytes is a positive integer budget for
estimated internal workspace. There is no automatic point deletion or reduction of N.

{pstd}Plan fields include m, n, d, reps, bounds, seeds, blocksize, level,
maxbytes, n_evaluations, output_bytes, external_bytes, input_block_bytes,
package_version, and layout. n_evaluations=R*N*(D+2); output_bytes=8*N*(D+2)
is one replicate's output storage, not total working memory.

{pstd}sobol_sa_inputs returns b.x with count*(D+2) rows and D columns and
b.keys with four columns: replicate, index, role, seed. replicate is 1-based,
index and start are 0-based, and count counts base points, not model evaluations.
role=0 means A, role=j means A with column j taken from B, and role=D+1 means B.
All coordinates are mapped into bounds. Each requested block contains all roles.
The current row order is role-major; always retain the keys when evaluating a model.
Count zero returns empty matrices, including at start=N. start+count must not exceed N.

{pstd}Within each replicate, A is the first D columns and B is the last D
columns of a single scrambled 2D-dimensional Sobol net. Using two identical
D-dimensional streams as A and B would be incorrect. The sampler can seek to
any base-point block without altering the plan. No Stata RNG calls are made.

{title:Model and analysis}

{pstd}f(X) must return a finite real column vector with one model value per
input row. It must be deterministic and invariant to batch size. The wrapper
generates blocks and retains only one replicate's N by (D+2) output matrix.
Arbitrary callback RNG calls and other side effects remain the caller's responsibility.

{pstd}sobol_sa_analyze accepts keys as a total-evaluations by 4 matrix and
values as a column vector. Every planned key must occur exactly once. Records
may arrive in any order, including across replicate boundaries. Invalid keys,
seed mismatches, duplicates, missing records, and missing values are errors.
The caller is responsible for associating values with the correct model inputs.

{pstd}Output columns within each replicate are A, AB_1,...,AB_D,B. Center all
outputs by their common mean after scaling by their maximum absolute value.
With V the population variance of pooled A and B outputs, compute
S1_j=mean(B*(AB_j-A))/V and ST_j=mean((A-AB_j)^2)/(2*V).
These are the Saltelli first-order and Jansen-form total-effect estimators used
by SALib v1.5.1. Constant/zero-variance sampled outputs are errors, not zero indices.

{pstd}Compute ratios within each scramble and then average these ratios across
replicates. Finite-sample estimates may lie outside [0,1] or violate S1<=ST;
values are not clipped or reordered.

{title:Results}

{p2colset 5 27 29 2}{...}
{p2col:S1, ST}1 by D mean indices{p_end}
{p2col:replicates}R by 2D; S1 columns followed by ST columns{p_end}
{p2col:covariance}2D by 2D sample covariance of replicate ratios divided by R{p_end}
{p2col:se}1 by 2D square roots of covariance diagonal{p_end}
{p2col:ci}2 by 2D; lower and upper approximate marginal t intervals{p_end}
{p2col:scaled_variance}R by 1 pooled variance after output scaling{p_end}
{p2col:output_scale}R by 1 maximum absolute output before centering{p_end}
{p2col:near_constant}1 if any scaled variance is below 1e-12; diagnostic only{p_end}
{p2col:outside_unit}1 by 2D flags for reported mean indices outside [0,1]{p_end}
{p2col:design}The complete input plan{p_end}
{p2col:estimator}Text identifier of estimator and aggregation rules{p_end}
{p2colreset}{...}

{pstd}The t intervals use R-1 degrees of freedom and between-scramble
variation. Coverage, finite-N bias, and simultaneous coverage are not certified.
Do not treat rows within a net as independent observations for an iid standard
error or bootstrap. A zero estimated SE is not proof of an exact index.
near_constant uses variance relative to max(abs(output)) squared; a large
additive offset can trigger it. It is a scale diagnostic, not a hypothesis test.

{title:Example}

{pstd}Run this block in a fresh Stata session after adding the repository ado
folder to the search path. It requires no dataset or external program.
{phang2}{cmd:sobol_setup}
{phang2}{cmd:mata:}
{phang2}{cmd:real matrix sa_demo(real matrix X) return(X*(1\2\3))}
{phang2}{cmd:bounds = (0,1\0,1\0,1)}
{phang2}{cmd:fit = sobol_sensitivity(&sa_demo(),12,bounds,20261000:+(1..16))}
{phang2}{cmd:fit.S1}
{phang2}{cmd:fit.ST}
{phang2}{cmd:fit.se}
{phang2}{cmd:end}
{pstd}The exact S1 and ST are (1,4,9)/14. Finite-sample estimates vary.

{title:External models and memory}

{pstd}sobol_sa_export writes a new request file with a versioned metadata
preamble, keyed inputs, and a completion trailer. It refuses an existing path.
The model must copy the preamble and keys unchanged and return a scalar y for
each row. sobol_sa_import checks the full plan, keys, row count, and trailer
before producing results. It allows reordered records. The file contract is
given below and can be implemented in the language used by the external model.
An incomplete write may leave a partial file, which is rejected on import.

{pstd}The analyzer/importer stores all R replicates' scalar outputs to support
arbitrary record order. Its memory requirement can exceed that of the callback
wrapper. Requests are exported in blocks. Input batches expand to (D+2) model
inputs per base point; high D may require smaller count/blocksize.
Workspace estimates are checked before large internal allocations and matrix
row counts are limited to 2^31-1. Estimates exclude caller-owned matrices,
Mata allocator overhead, loaded generator state and callback workspace; maxbytes
is not a hard process memory limit. Review the plan before an expensive model run.

{pstd}Restart Stata when upgrading from an earlier package version. Old plans
and exchange files retain their package version; use that version to process
them. This interface does not identify or replace the historical LPTAU routine
or its ABI.

{title:External-model file contract}

{pstd}For example, create a plan and export a new request file in the current
working directory. This block creates request.csv and refuses an existing path.
Choose a new filename for another run.
{phang2}{cmd:sobol_setup}
{phang2}{cmd:mata: p=sobol_sa_plan(10,(0,1\0,1\0,1),20261000:+(1..16))}
{phang2}{cmd:mata: p.n_evaluations}
{phang2}{cmd:mata: p.external_bytes}
{phang2}{cmd:mata: sobol_sa_export(p,"request.csv")}

{pstd}The request and response have the same preamble. The following notation
shows field meanings; copy the actual exported preamble as text without editing:
{phang2}{cmd:#sobol_sa,1,0.4.2}
{phang2}{cmd:#plan,m,D,R,level}
{phang2}{cmd:#bounds,j,lower_j,upper_j} (one line for each dimension)
{phang2}{cmd:#seed,r,seed_r} (one line for each replicate)

{pstd}The request column header is replicate,index,role,seed,x1,...,xD.
The response header must be exactly {cmd:replicate,index,role,seed,y}.
For every request row evaluate the specified deterministic model on its x values
and return one finite scalar y with the same four keys. Records may be reordered,
but each key must occur exactly once. Do not separate values from their keys.
Preserve the request's completion trailer, {cmd:#end,total_evaluations}, where
total_evaluations=R*N*(D+2). The response ends immediately after this line;
a normal final newline is allowed.

{pstd}Use unquoted decimal numbers with a period as decimal separator and optional
E/e exponents. Do not use locale decimal commas, D exponents, NaN/Inf, missing
values, extra columns, blank lines, or comments after the column header.
Each response record has exactly five nonempty numeric fields. Keep integer keys
exact. Preamble numbers were exported at 17 significant digits; retain them as text.
The importer binds m, bounds, seeds, level, layout and package version to the plan.
Block size and memory budget are operational settings and are not in the preamble.
A matching preamble does not prove that the external model evaluated the right inputs.

{pstd}After the external model has produced response.csv, use the original plan
or recreate exactly the same plan in a fresh session:
{phang2}{cmd:mata: fit=sobol_sa_import(p,"response.csv")}
{phang2}{cmd:mata: fit.S1}
{phang2}{cmd:mata: fit.ST}
{phang2}{cmd:mata: fit.se}

{title:Also see}
{pstd}{help mf_sobol}, {help mf_sobol_integrate}, {help sobol}
