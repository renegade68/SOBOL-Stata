# Sobol sequences for Stata

`sobol` generates Sobol sequence coordinates in Stata. The Mata functions also
provide stateful streams, repeated scrambled integration, and first-order (S1)
and total-effect (ST) sensitivity indices for independent uniform inputs.

Version **0.4.2** uses pure Stata/Mata and requires **Stata 16 or newer**.
No Python, R, Fortran, or compiled plugin is needed. This is experimental
software, tested on Stata/MP 16.1 for macOS; other platforms and Stata releases
have not been verified.

## Start using the commands

Use GitHub's **Code → Download ZIP**, extract it, and keep **all 11 files in
`ado/` together**. The five `.mata` files are required implementation files;
`sobol_setup` loads them automatically. They are not separate examples.

In Stata, replace the path below with the extracted folder on your computer:

```stata
adopath ++ "/your/path/SOBOL-Stata/ado"
which sobol
help sobol
```

The path applies to the current session; include it in your analysis do-file
for later sessions. After an update, replace the complete `ado/` folder and
restart Stata. Use `which sobol` and `which sobol_setup` to check for an older
installation taking precedence. This repository uses local-folder installation;
it does not provide `ssc install` or `net install` metadata.

The examples below create their own inputs. Save current work before running
them: `clear` removes the dataset, and `clear all` also removes loaded programs
and Mata objects. The `adopath` setting remains available.

## Generate points

Create a reproducible scrambled 1,024-point net in three dimensions:

```stata
clear
sobol u1-u3, n(1024) scramble seed(20261023)
list u1-u3 in 1/8
```

Each variable is stored as `double`. Omit both `scramble` and `seed()` for the
unscrambled sequence. Its first point is the origin. If observations already
exist, their number is used; an explicit `n()` must agree with that number.
Existing variables are retained, and new variable names must be unused.

For a Mata matrix, after loading the functions:

```stata
sobol_setup
mata: U = sobol_base2(10,3,20261023)
mata: assert(U==st_data(.,("u1","u2","u3")))
```

Here `10` means `2^10 = 1,024` points. This assertion uses the variables from
the preceding example. For large calculations, `sobol_init()` and `sobol_draw()`
generate consecutive blocks without storing the full matrix; see the
[streaming help](ado/mf_sobol.sthlp).

## Integrate a function

This integrates `exp(u1+u2+u3)` over the unit cube using 32 scrambles, each with
4,096 points. A callback returns one row per input point; multiple output
columns represent multiple integrals.

```stata
clear all
sobol_setup
mata:
real matrix sobol_demo_integrand(real matrix u) return(exp(rowsum(u)))
fit = sobol_integrate(&sobol_demo_integrand(),12,3,20261000:+(1..32))
fit.estimate
fit.se
fit.ci
(exp(1)-1)^3
end
```

The final expression is the exact integral, approximately **5.073214**.
`fit.se` and the approximate t interval use variation between scrambles.
They describe numerical integration uncertainty, not uncertainty in estimated
model parameters. They do not guarantee coverage for every integrand.

## Calculate sensitivity indices

For independent uniform inputs on `[0,1]`, this example evaluates
`x1 + 2*x2 + 3*x3`:

```stata
clear all
sobol_setup
mata:
real matrix sobol_demo_model(real matrix X) return(X*(1\2\3))
bounds = (0,1\0,1\0,1)
fit = sobol_sensitivity(&sobol_demo_model(),12,bounds,20261000:+(1..16))
fit.S1
fit.ST
fit.se
end
```

The exact S1 and ST indices are both `(1,4,9)/14`, approximately
**(0.071429, 0.285714, 0.642857)**. Finite-sample estimates vary. The callback
must return one finite scalar per input row, deterministically and independently
of block size. The method uses `R*N*(D+2)` model evaluations; inspect
`sobol_sa_plan()` before an expensive run.

## Choose an interface

| Task | Interface | Help |
| --- | --- | --- |
| Create Stata variables | `sobol` | [syntax and options](ado/sobol.sthlp) |
| Matrices, blocks, indexed access, recorded scrambles | `sobol_base2()`, `sobol_init()`, `sobol_draw()`, `sobol_seek()` | [Mata generator](ado/mf_sobol.sthlp) |
| Repeated integration | `sobol_integrate()` | [integration](ado/mf_sobol_integrate.sthlp) |
| S1/ST, separate sampling/analysis, external model CSV exchange | `sobol_sensitivity()`, `sobol_sa_*()` | [sensitivity](ado/mf_sobol_sensitivity.sthlp) |

Use `help mf_sobol`, `help mf_sobol_integrate`, or `help mf_sobol_sensitivity`
in Stata for formatted help. Examples and the external-model file contract are
included there; no separate example script or dataset is needed.

## Numerical conventions and limits

- Generation supports **1–1,024 dimensions**; sensitivity supports **1–512
  independent uniform inputs** and deterministic scalar output. Correlated-input
  indices, S2, and stochastic-model noise decomposition are not implemented.
- Coordinates use 32-bit numerators on `[0,1)`. Prefer complete prefixes of
  length `2^m`; arbitrary skipping or thinning can lose their balance properties.
  A single scramble's points are dependent observations.
- Zero can occur, including with scrambling. Define an endpoint strategy before
  applying `log()` or an inverse distribution function. No automatic clipping or
  point deletion is applied. Keep coordinates as `double`.
- The explicit seed reproduces this implementation's points without changing
  Stata's RNG. The same numeric seed in another library need not give the same
  scramble. Matching direction numbers, ordering, LMS and shift is required.
- Sensitivity estimates are not clipped to `[0,1]`. Finite-N bias and approximate
  interval coverage require judgment. Sampled zero output variance is an error.
  Memory estimates are workspace budgets, not hard limits on process memory.
- Old sensitivity plans/CSV files retain their package version; process them
  with that version rather than relabeling them.
- The package is not a drop-in replacement for the historical Fortran LPTAU
  routine discussed on Statalist.

Checks include aligned SciPy sequence references, streaming and RNG behavior,
integration examples, analytic sensitivity models, and external-model exchange.
They do not establish accuracy or speed for every application.

## Sources and related Stata work

| Role in this package | Source |
| --- | --- |
| Direction numbers and generator recurrence | [Joe and Kuo's data and algorithm notes](https://web.maths.unsw.edu.au/~fkuo/sobol/), including [Joe and Kuo (2008)](https://doi.org/10.1137/070709359) |
| Private MT19937 adaptation used for scrambling | [Matsumoto and Nishimura's 2002 reference code](https://math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html); its notice is retained in `ado/_sobol_mt19937.mata` |
| Sequence validation reference | [SciPy 1.13.1 `qmc.Sobol`](https://docs.scipy.org/doc/scipy-1.13.1/reference/generated/scipy.stats.qmc.Sobol.html), with aligned direction numbers, bit depth and scramble configuration |
| Sensitivity estimators and validation reference | [Saltelli et al. (2010)](https://doi.org/10.1016/j.cpc.2009.09.018) and [SALib 1.5.1's Sobol analyzer](https://github.com/SALib/SALib/blob/6fb0a9c144692be5993e50dc99575eab44c192be/src/SALib/analyze/sobol.py); see the [SALib citation guidance](https://salib.readthedocs.io/en/latest/citations.html) |

Related Stata/Mata quasi-Monte Carlo work includes
[Drukker and Gates (2006)](https://www.stata-journal.com/article.html?article=st0103)
on Halton/Hammersley sets and [Kolenikov (2012)](https://doi.org/10.1177/1536867X1201200103)
on scrambled Halton sequences. These are related methods, not Sobol code used by
this package. No claim of being the first Stata Sobol implementation is made.
SciPy and SALib are development references, not runtime dependencies; their
seeding and confidence-interval conventions need not match this package.

## Feedback and license

For a problem report, use [GitHub Issues](https://github.com/renegade68/SOBOL-Stata/issues).
Include the Git commit, Stata version, `which sobol` output, complete command,
expected versus actual results, and a small reproducible example. Synthetic
inputs are sufficient.

Project code is provided under the [MIT License](LICENSE). The Joe–Kuo direction
data have their [own BSD-style notice](JOE_KUO_LICENSE.txt); the MT19937 adaptation
retains the Matsumoto–Nishimura notice in its source file. Preserve these notices.
The help files describe the sequence and sensitivity methods.
