# Reproduction package: two estimations (gen1, taubar=0.14, chi=0.5)

Snapshot copied from `C:\Users\gawater\Documents\GitHub\Shotgun` on 2026-08-08,
corresponding to the results documented in
`mmt8bOccBin3_covid_bT_baa_epsL_summary.tex`/`.pdf` in this folder.

## Contents

**Document**
- `mmt8bOccBin3_covid_bT_baa_epsL_summary.tex` / `.pdf` — the write-up.
- `prior_posterior_covid_bT_baa.pdf`, `det_estab_map_mmt8bOccBin3_covid_bT_baa.pdf` — the two figures the `.tex` includes.

**Estimation model** (in `mmt8bOccBin3_covid_bT_baa/`)
- `mmt8bOccBin3_covid_bT_baa.mod` — main estimation file. Its final `estimation(...)`
  block is set to the **accepted (gen1)** configuration: `mode_compute=0`,
  `mode_file` pointing at the Stage-1 mode below, `mh_jscale=0.25`,
  `mcmc_jumping_covariance='posterior_jumping_covariance'` (gen1's file, not gen2's).
- `fred_data_debt_baa.mat` — the actual estimation dataset (`FYGFDPUN`-based debt/GDP,
  BAA-based loan rate, log-demeaned throughout).
- `posterior_jumping_covariance.mat` — gen1's custom MCMC jumping covariance
  (built from an earlier hessian-covariance round's own pooled posterior draws;
  see `build_posterior_jumping_covariance.m` and Section 2.3 of the `.tex`).
- `mmt8bOccBin3_covid_bT_baa/Output/mmt8bOccBin3_covid_bT_baa_mode.mat` — the
  Stage-1 (linear-KF) posterior mode that `mode_file` points to. Required:
  `mode_compute=0` has nothing to optimize from without it.
- `run_stage1_epsL.m`, `run_stage2_occbin.m`, `run_production.m` — MATLAB
  batch wrappers (paths already updated to this folder). Toggle the model's
  `occbin_val` flag and edit the literal `mode_compute`/`mh_replic`/`mh_jscale`
  values directly in the `.mod` file for a from-scratch two-stage run
  (Dynare requires literal tokens there, not variables — see the `.mod`'s
  own comments).

**Determinacy / E-stability**
- `mmt8bOccBin3_covid_bT_baa_bkonly.mod` — RELAX-regime-only companion model
  (no `occbin_constraints`, no estimation) used for the Blanchard-Kahn/
  E-stability grid.
- `det_estab_mmt8bOccBin3_covid_bT_baa.m` — runs the grid, loads posterior
  means from the saved results `.mat`, produces the combined determinacy/
  E-stability figure and the `phib_val` HPD sensitivity check.

**Posterior figure**
- `plot_prior_posterior_covid_bT_baa.m` — builds the prior/posterior density
  figure from the raw MCMC draw files (`.../metropolis/*_mh1_blck*.mat`).
  Note: those raw draw files are **not** included here (they're large and
  were also overwritten in the source repo by a later gen2 run); this script
  will need a completed estimation run's `metropolis/` output to regenerate
  the figure. The already-rendered PDF is included above.

**Data provenance** (root level)
- `build_debt_dataset_fygfdpun.m` — builds `bT_obs` from FRED `FYGFDPUN`
  (debt held by the public) and `GDP`, log-demeaned. This is the *current*
  data-build script; it fetches directly from FRED via `curl` (no API key).
- `fetch_fred_data_baa.m` — builds the BAA-based `rl_obs`, requires `BAA.csv`
  (included) and `fred_data_debt.mat` (included) as inputs.
- `build_debt_dataset.m`, `grossdebt.csv`, `fred_data_occbin.mat`,
  `fred_data_debt.mat` — earlier links in the data-provenance chain, kept
  for full traceability even though the current `fred_data_debt_baa.mat`
  (already built) is what estimation actually reads.
- `build_posterior_jumping_covariance.m` — documents how gen1's
  `posterior_jumping_covariance.mat` was built, and the critical Dynare
  semantics (the saved matrix must be the *precision*, i.e.
  `inv(Sigma_posterior)`, not the covariance itself — Dynare inverts
  whatever is loaded under the variable name `jumping_covariance`).

## Reproducing the estimation

The included artifacts (mode file + jumping covariance file) let you
re-run the **production MCMC step directly** without redoing Stage 1/Stage 2:

```matlab
addpath('C:\dynare\7.1\matlab');
cd('C:\Users\gawater\Documents\GitHub\Shotgun 2\mmt8bOccBin3_covid_bT_baa');
run('run_production.m')
```

**Caveat**: Dynare's MCMC sampler is not seeded in these scripts, so a fresh
run will not reproduce bit-identical draws — but starting from the same
mode, data, and jumping covariance, it should reproduce the same posterior
distribution and the same qualitative conclusions (acceptance ratios in the
0.1-0.3 range, determinacy/E-stability results, etc.), matching what's
reported in the `.tex`.

To reproduce fully from scratch (re-running Stage 1 and Stage 2 too, e.g.
after further changing the calibration), see the full jscale-tuning
history and kill/re-tune procedure documented inline in the `.mod` file's
comments and in Section 2 of the `.tex`.
