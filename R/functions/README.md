# R conversions for `MCSVAR_HS_SV_A`

This directory contains an R translation of the MATLAB routine implemented in
`functions/MCVAR/MCVAR_HS_SV_A.m` (whose function declaration is
`MCSVAR_HS_SV_A`) and the helper routines that it calls.

From the repository root, load the converted function with:

```r
source("R/functions/load_mcsvar_hs_sv_a.R")
```

The main entry point is:

```r
Y_pred <- MCSVAR_HS_SV_A(data, reps, burnin, N, G, L, h)
```

For convenience, `MCVAR_HS_SV_A` is also provided as an alias for the same
function.

The translation uses base R matrix routines only; no external R packages are
required.
