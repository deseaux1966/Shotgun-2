# rn measurement-error prior / mh_jscale investigation — summary (2026-09-19 to 2026-09-26)

## Starting point
Section 2.4's regime-classification check found the OccBin filter/smoother flag
2017Q1-2019Q4 as ZLB-binding, which conflicts with the data (positive observed
policy rates). Root cause (confirmed via rule decomposition): a strong,
non-inertial Taylor rule (bpi=3.47) with a small cushion above the floor
(steady-state rn=1.0022) pushes the notional rate below 1 whenever inflation
dips modestly below its sample mean. The `rn` measurement-error shock absorbs
the resulting misfit (~1 posterior SD in 2017-2019).

## Attempts tried, in order

1. **Interest-rate smoothing (Version A)**: `rn_notional = rn(-1)^rho * [rule]^(1-rho)`,
   lag applied only in the relaxed regime. Linear-filter Stage 1 strongly favored
   heavy smoothing (rho~0.95-0.97) and fixed the misclassification in the linear
   check. But the OccBin likelihood is infinite at that mode and only finite in a
   narrow (rho, bpi) band; a profile search there topped out at logpost 1260.0
   vs. 1266.2 for the original no-smoothing OccBin mode, and the resulting regime
   path only recovered 21 of 41 true ZLB quarters. Abandoned — files:
   `mmt8bOccBin3_covid_bT_baa_smooth*.mod`.

2. **Tighter prior on stderr eps_rn_me**: fixed-parameter sweeps showed the
   OccBin log-posterior kept rising as this SD fell from 0.0029 toward ~0.0005,
   then a profile re-optimizing (bpi, by, stderr epsr, ali, br) at SD=0.0015
   recovered the exact true ZLB dates (2009Q1-2016Q4, 2020Q2-2022Q2) with only
   a ~20-32 log-point fit cost. SD=0.001 did not work (re-optimized logpost
   1184.7 vs. 1266.6 baseline; likelihood infinite at/near the optimum).

3. **Stage 2 with InvGamma(0.0012, 0.0004) prior on stderr eps_rn_me**,
   mh_jscale=0.05→0.15→0.18 (diagnostics), then **switched to Uniform(0.0002,0.002)**
   after the inverse-gamma's lower tail looked like it might be doing the work.
   Diagnostics under the uniform prior (800-3000 draws) at jscale=0.08, 0.12, 0.15
   all looked reasonable or fixable.

4. **Full 50,000-draw x 2-chain production runs** (the actual test that mattered):
   - **jscale=0.18**: 3k-draw diagnostic looked clean (28.6%/24.1% acceptance,
     interior eps_rn_me=0.0006). Full production (chunked, ~9h) collapsed:
     acceptance fell to 7.4%/8.8%, eps_rn_me pinned at its prior floor (0.0002,
     HPD collapsed to a point) by ~10,000 draws in. Geweke 5/8, 2/8 at 15% taper.
   - **jscale=0.15**: pinned even faster — already at the 0.0002 boundary by
     3,000 draws (never reached production).
   - **jscale=0.12**: looked clean at both 800 and 3,000 draws (interior
     eps_rn_me~0.0006, mild acceptance decline 49%->42-45%). Full 50,000-draw
     production (2026-09-25/26, chunked, ~7.5h) **also pinned**: final
     acceptance 14.2%/16.6% (chain 1 below the 15% floor), eps_rn_me collapsed
     to [0.0002,0.0002] again, Geweke chain 1 only 3/8 passed at 15% taper
     (eta, by, br, rhoz, phib_val all failed).

## Conclusion (as of 2026-09-26)
Every jscale tried (0.12, 0.15, 0.18) eventually pins stderr eps_rn_me at the
Uniform(0.0002,0.002) prior's lower bound given enough draws — only the draw
count at which it happens differs (3k for 0.15; ~10k for 0.18; somewhere
between 3k and 50k for 0.12). Short diagnostics are NOT predictive of whether
a jscale will survive to 50,000 draws. This looks like a structural property
of the OccBin likelihood surface near this boundary, not a proposal-tuning
problem fixable by picking a better jscale.

**Not yet done**: a direct likelihood sweep at SD < 0.0002 using the ACTUAL
posterior-mean parameters from the jscale=0.12 production run (an earlier
sweep at SD=0.0001-0.000001 used rounded/approximate parameter values and
showed catastrophic likelihood collapse below 0.0002 — but that check should
be redone with this run's real posterior mean to settle whether the pinning
reflects a genuine preference for SD<0.0002 or a sampler artifact at the
boundary regardless of true value).

## Recommendation on the table (unresolved)
Two options were being weighed when this investigation paused for a computer
restart:
1. Redo the near-zero likelihood sweep with the jscale=0.12 run's actual
   posterior means, to determine once and for all whether SD<0.0002 is
   genuinely favored.
2. Abandon the tightened-prior approach entirely; report the original
   accepted estimation (log data density 1158.40, Table 5 of the summary
   .tex) as the main result, and document the 2017-2019 ZLB misclassification
   plus this whole investigation as an unresolved limitation in Section 2.4.

## Key files
- `mmt8bOccBin3_covid_bT_baa_rnu012_prod.mod` — the jscale=0.12 production
  model (Uniform(0.0002,0.002) prior). Output/chains in
  `mmt8bOccBin3_covid_bT_baa_rnu012_prod/` (Output/ and metropolis/
  subfolders) — NOT deleted, kept on disk.
  `mmt8bOccBin3_covid_bT_baa_rnu012_prod_results.mat` in Output/ has the
  final oo_/M_ posterior structures.
- `mmt8bOccBin3_covid_bT_baa_rnu.mod`, `..._j08.mod`, `..._j12.mod`,
  `..._j12b.mod`, `..._j15.mod` — the standalone diagnostic runs referenced
  above (3000-draw comparisons at each jscale).
- `s2_rn/mmt8bOccBin3_covid_bT_baa_rn_start_mode.mat` — the shared start point
  (OccBin posterior means with bpi, by, stderr epsr, ali, br re-optimized at
  stderr eps_rn_me=0.0015) used as `mode_file` for all of the above.
- `profile_rnsd.m`, `profile_rnsd05.m` — the fixed-parameter SD sweep/profile
  scripts from step 2 above.
- Section 2.4 of `mmt8bOccBin3_covid_bT_baa_epsL_summary.tex` documents the
  ZLB misclassification finding and the smoothing-attempt writeup (step 1
  above); it does NOT yet reflect the Uniform-prior/jscale saga (steps 2-4).
