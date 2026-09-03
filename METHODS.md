# A Blind Transit-Detection Pipeline with an Error-Correcting-Code Contamination Screen

**Status:** exploratory research tool, validated on n=2 real targets plus synthetic
stress tests, including a full pipeline re-validation (Section 12) that found and
fixed two additional real bugs beyond the original development history.
Not yet benchmarked against established pipelines (BLS/TLS/wotan).
Appropriate for a Zenodo software release; the contamination-screening result is
narrow and honest enough for an RNAAS research note; not yet ready for a full
peer-reviewed ApJ/AJ submission (see Limitations).

## 1. Summary

This project builds, from first principles, a box-least-squares-style transit
detector for TESS photometry, with a bootstrap-calibrated significance
threshold and a novel contamination-screening stage based on Reed-Solomon-style
algebraic error location (the Berlekamp-Welch algorithm), informed by TESS's
own SPOC quality flags. Every claim below was verified computationally against
real TESS light curves or controlled synthetic tests built from real noise —
including several points where an initial implementation was wrong, caught by
testing against a case with a known right answer, and fixed. Those failures are
reported here alongside the successes, because they are part of what makes the
final numbers trustworthy.

## 2. Data

- **WASP-18 b** (TIC 100100827): TESS Sectors 2-3, SPOC 2-minute cadence,
  31,154 quality-zero rows. Public MAST data. Known parameters: P=0.94145 d,
  depth ~0.9-1.0%.
- **Pi Mensae c** (TIC 261136679): TESS Sectors 1, 4, 8, SPOC 2-minute cadence,
  45,841 quality-zero rows. Public MAST data. Known parameters: P=6.2679 d,
  depth ~300 ppm (Huang et al. 2018).

Both packages were checksum-verified (SHA-256) on receipt.

## 3. Core detection method

- **Detrending:** per-sector rolling median, 0.5-day window.
- **Noise estimate:** robust (MAD-based) scatter of the detrended residual, not
  the packaged flux-uncertainty column — which was found to be corrupted in
  both datasets (WASP-18: mis-scaled by ~5 orders of magnitude, consistent with
  a divide-by-median-squared bug; Pi Mensae: identically zero for every row).
- **Search statistic:** a box-least-squares approximation — phase-fold at each
  trial period, bin into 100 phase bins, scan box widths of 2-15 bins, take the
  maximum-SNR box.
- **Significance:** a **bootstrap null**, not a fixed SNR cutoff — circularly
  shift the residual, rerun the identical search, repeat ~20-30 times, and use
  the 99th percentile of the resulting SNR distribution as the detection
  threshold. This was necessary: naive SNR thresholds were shown empirically to
  produce false positives once the search space grew large enough (Section 6).

## 4. Result 1 — WASP-18 b (positive control)

Blind search (no ephemeris used) recovered:

| Quantity | Recovered | Published |
|---|---|---|
| Period | 0.94150-0.94160 d | 0.94145 d |
| Depth | 9,880-9,979 ppm | ~9,000-10,000 ppm |
| SNR | 758-765 | — |

Robustness checks: independent detection in both sectors separately
(Sector 2: SNR 571; Sector 3: SNR 499); odd/even transit depths consistent
(9,913 vs 9,850 ppm — rules out a background eclipsing binary at half the
period); true period beats both P/2 and 2P aliases; small (~460 ppm) secondary
eclipse signal, consistent with the star's known phase curve rather than a red
flag.

## 5. Sensitivity floor (injection-recovery)

The real WASP-18 transit and its ~225 ppm phase-curve modulation were both
explicitly modeled and removed (fitted 1x/2x-orbital-frequency sinusoid, plus
deletion — not patching — of a 3-hour window around each transit; patching
with a constant offset was tried first and left a detectable ~3,900 ppm
remnant, which is why deletion was used instead). Synthetic transits of
varying depth were then injected into this real, cleaned noise:

| Depth | Recovery |
|---|---|
| ≥500 ppm | 100% |
| 350 ppm | 50% (detection threshold) |
| ≤250 ppm | 0% |

Noise-only false-alarm threshold (bootstrap, cleaned data): SNR 9.2.

## 6. Result 2 — Pi Mensae c (harder target, contaminated)

Chosen because its published depth (~300 ppm) sits close to the measured
sensitivity floor above. The blind search's first result was wrong: **P=10.34 d,
SNR 105.3** — nowhere near the known 6.27 d period. Per-sector breakdown
immediately showed why: depths of -11, 1109, and 484 ppm across the three
sectors — wildly inconsistent, the signature of a single non-periodic event
rather than a real signal. A diagnostic plot located a single ~8,000 ppm,
20-40-minute event in Sector 4 (t=1421.288 BTJD). Deleting a 3-hour window
around it and rerunning:

| Quantity | Recovered | Published |
|---|---|---|
| Period | 6.26842-6.26594 d | 6.2679 d |
| Depth | 175-315 ppm | ~300 ppm |
| SNR | 38.5-46.6 | — |

Per-sector depths after cleaning: 190, 249, 212 ppm — consistent, the expected
signature of a real periodic signal.

Depth estimates from the box-search statistic are consistently biased low
relative to the literature value (a known limitation of a box fit against a
curved, limb-darkened transit profile — see Limitations).

## 7. Contamination screening: three designs, tested and compared

### 7a. Complementary-pair screen (32x32 most-perfect magic square)

A 32x32 most-perfect magic square was constructed (Ollerenshaw-Brée
reversible-square method) and verified computationally to satisfy every
defining property (row/column/diagonal magic sums, all 2x2 block sums, both
complementary-pair diagonal directions). Its exact algebraic guarantee — the
complement of value v is always 1025-v — was used to split 1024-point blocks
of the light curve into two matched halves and flag blocks with large
first-half/second-half asymmetry.

- Noise-only threshold (bootstrap): 3,083 ppm.
- **False-flag rate on a real periodic transit (WASP-18): 29/30 blocks (97%).**
  This makes it unsafe as an automatic exclusion filter — it cannot distinguish
  a real recurring signal from a single glitch.

### 7b. Berlekamp-Welch error-locator

The classical algorithm was verified exactly on a noiseless synthetic case
(planted errors located to machine precision, recovered polynomial exact)
before being adapted for real, noisy data: fit a local low-degree polynomial
plus one assumed error per ~3-hour window, and score each window by the
residual (in ppm) between the data and the recovered polynomial at the
worst-fitting point.

An initial scoring statistic (a ratio of algebraic values) was numerically
unstable — it scored an obvious, deliberately-injected 8,000 ppm anomaly as
"not significant" while occasionally producing enormous scores from pure
noise. This was caught by testing against synthetic cases with a known right
answer and fixed by switching to a physically interpretable ppm-residual
score.

- Noise-only threshold: 2,431 ppm.
- **False-flag rate on WASP-18: 80/342 windows (23%)** — worse than ideal, but
  a 4x improvement over 7a, and unlike 7a it reports an exact flagged point,
  not just a block.
- Sensitivity floor: ~2,000-3,000 ppm for a single ~1-hour event.

### 7c. Hybrid: SPOC quality flags as a prior, Berlekamp-Welch for exact location

Investigation of the raw FITS `QUALITY` column revealed that TESS's own SPOC
pipeline had already flagged cadences near (though not exactly at) the known
Pi Mensae anomaly with bit 512 (`ImpulsiveOutlier`) — a flag that exists in
the official TESS Data Products documentation but is deliberately excluded
from SPOC's own recommended default exclusion bitmask.

Using `ImpulsiveOutlier` cadences to generate a small number of candidate
regions, then running Berlekamp-Welch only inside those regions:

- Pi Mensae: **26 candidate regions** (vs. 503 blocks for a blind scan — a 19x
  reduction in search cost). The single highest-scoring candidate (4,343 ppm,
  ~7x the runner-up) located the peak at **t=1421.2880**, matching the true
  event to four decimal places.
- WASP-18: **0 candidate regions, 0 false positives** — but this is partly a
  favorable setup rather than proof of discrimination, since `ImpulsiveOutlier`
  never fires in WASP-18's data at all (see Section 8).

## 8. Stress test: does the screen survive a flag that coincides with a real transit?

A controlled synthetic scenario was built from real noise: a genuine periodic
transit (P=2.17 d, depth 3,000 ppm, 12 events) plus a standalone 6,000 ppm
contaminant, plus a SPOC-style flag placed deliberately on **one of the real
transit events** (simulating an unlucky coincidence) and on the standalone
contaminant.

- **v1 (naive hybrid):** excised both — including 76 of 1,169 real transit
  cadences (6.5% of all transit data), costing 15.4 SNR on the eventual
  detection.
- **v2/v3 (periodicity-aware, independent alias search):** correctly spared
  the real transit (found its true period, 2.16-2.17 d) but **also incorrectly
  spared the genuine contaminant** — its own periodicity search was fooled by
  an alias of the real transit's period (≈P/4) appearing to "recur" at the
  contaminant's location, purely by chance overlap, unrelated to the
  contaminant itself. Requiring individually-significant, well-separated
  cycles (v3) did not fix this — the individual cycles really were
  significant, just aliasing someone else's signal.
- **v4 (structural fix):** rather than inventing a second, independent
  periodicity test, each candidate is checked against the period the main BLS
  search has *already found and bootstrap-validated*. This correctly spared
  the real transit **and** correctly excised the contaminant, with **zero**
  real transit data lost and a small net SNR **gain** (+2.6) from removing the
  genuine contaminant.

## 9. Limitations (honest accounting)

- **Real-world validation sample is n=2** targets, both already-published,
  known planets. This is calibration, not a validated general-purpose survey
  tool.
- **The screening pipeline's real-anomaly validation is n=1** (the Pi Mensae
  Sector 4 event). Its behavior on other contamination types — cosmic rays,
  momentum-dump ramps, scattered light, stellar flares — is untested.
- **No baseline comparison** to established tools (`astropy.timeseries.BoxLeastSquares`,
  TLS, `wotan`-style detrending, iterative sigma-clipping). The claims here are
  relative to this project's own earlier iterations, not to the field's
  standard tools.
- **Depth precision is coarse.** The box-search statistic systematically
  underestimates depth against a real, curved, limb-darkened transit profile.
  A production tool would need a proper transit model (e.g. `batman`,
  `exoplanet`) for parameter estimation, not just detection.
- **The v4 structural fix requires an already-validated periodic signal to
  reference.** For single-transit or few-transit targets — the hardest and
  most scientifically interesting case — there is no prior validated period to
  check candidates against, and this problem is unsolved here.
- **The Berlekamp-Welch screen's standalone false-flag rate (23%) is still
  non-trivial** without the SPOC-flag prior; the hybrid design's strong
  results depend on the prior being available and reasonably complete.

## 10. Reproducibility

All code referenced here — the BLS search, bootstrap calibration, magic-square
construction and verification, Berlekamp-Welch implementation and its
synthetic verification, the hybrid screen, and the overlap stress test — was
written, run, and its output captured over the course of this project. Scripts
are included in this repository/release; each analysis script reads directly
from the packaged TESS CSV/FITS data and reproduces the numbers quoted above.

## 11. Suggested next steps

1. Benchmark against `wotan`/TLS on the same two targets to establish whether
   this pipeline adds anything existing tools don't already do better.
2. Test the contamination screen against real (not synthetic) examples of
   other TESS quality-flag categories (`Desat`, `Straylight`, `CollateralCosmic`).
3. Address the single-transit chicken-and-egg problem for the v4 structural
   fix — perhaps via a population-level prior (typical hot-Jupiter/sub-Neptune
   duty cycles) rather than a per-target validated period.
4. Replace the box-search depth/duration statistic with a proper transit
   model for any result intended for publication-quality parameter estimates.

## 12. Update: end-to-end re-validation exposed two more real bugs

Section 8's structural fix (v4) was validated on a *synthetic* overlap
scenario. Running the same fix on the *real* Pi Mensae dataset end-to-end
(the natural next test) failed: the blind search recovered P=8.889 d, not
the correct 6.267 d. This section documents why, because the failure and
its fix are as informative as the earlier successes.

### 12a. Root cause: the quality-flag screen has a blind spot

Isolating each previously-known contaminant one at a time showed a **second
real event, near t=1535 BTJD**, was necessary to exclude to recover the
correct period -- and it was invisible to the v4 screen entirely. Raw-flux
inspection confirmed it: a genuine ~4-6 hour, ~2,900-4,000 ppm dip-then-
overshoot, the same general shape as the Sector 4 anomaly (Section 6) and
almost certainly the same class of event. But SPOC's `ImpulsiveOutlier`
(bit 512) never fires there -- only `quality=0`, plus one unrelated
`ManualExclude`/`CoarsePoint` flag at the window's edge.

`ImpulsiveOutlier` is tuned for sharp, few-cadence, cosmic-ray-like spikes,
not slower multi-hour ramps. A screen that *only* looks where that specific
flag points will systematically miss this category of contamination --
this is a structural limitation of the design in Sections 7c/8, not a bug
in the Berlekamp-Welch localization itself (which, given the chance, can
find such events fine -- it never got the chance here).

### 12b. Fix: merge two candidate sources instead of relying on one

Rather than trusting the quality-flag prior alone, candidates are now drawn
from **both** sources and merged (deduplicating events found by both):

- **Source A (quality-flag prior):** as in Section 7c -- high precision,
  category-limited.
- **Source B (blind self-calibrated scan):** every fixed-size block is
  scored directly by Berlekamp-Welch, with a bootstrap-calibrated
  per-dataset threshold. Lower precision alone (this is the same design
  characterized in Section 7b), but no blind spot tied to a specific SPOC
  flag category.

Re-running end-to-end with the merged candidate list and the same
structural periodicity check (Section 8) recovered **P=6.266-6.270 d**
(literature: 6.2679 d) correctly, while WASP-18 b remained correctly
detected and essentially unchanged (its light curve has zero
`ImpulsiveOutlier` flags at all, so source A never contributed there).

### 12c. Two more bugs found while validating the fix

**Numerical instability at sector boundaries.** The blind scan (source B)
produced physically nonsensical scores -- 2,745,173 ppm and 8,992,114 ppm --
at exactly the last cadence of Sectors 1 and 4. A block spanning the
multi-day inter-sector gap gives Berlekamp-Welch a wildly uneven x-spacing,
and the polynomial division can blow up. Fixed with a stability guard:
detect windows with either (a) an internal gap more than 3x the local
median cadence, or (b) an ill-conditioned design matrix (`cond > 1e8`), and
fall back to a simple robust (MAD-based) deviation statistic instead of
trusting the algebraic fit -- with a physical sanity cap (no real dip
exceeds 50%) applied either way, belt-and-suspenders.

**A null-calibration mismatch that the guard exposed.** After adding the
guard, the blind scan's candidate count *dropped* unexpectedly, losing two
genuine contaminants along with the two numerical artifacts. Cause: the
null-threshold calibration sampled *scattered random indices* across the
entire baseline to build "noise-only" windows -- which have enormous
internal gaps by construction, so every null trial tripped the new gap
guard, while real contiguous blocks mostly didn't. That mismatch was
latent even before the guard existed; the guard just made its effect
visible. Fixed by calibrating the null from genuinely contiguous blocks
(circularly shift the residual, take a real contiguous slice), matching
the bootstrap technique already validated and used elsewhere in this
pipeline (Section 3). With both fixes in place, the blind scan correctly
recovered the t=1535 contaminant via the `robust_fallback` path, and the
final result matched the literature period.

### 12d. Final, re-validated results

| Target | Recovered | Published | Notes |
|---|---|---|---|
| WASP-18 b | P=0.9415-0.9416 d, depth 9,900-10,000 ppm, SNR 729-765 | P=0.94145 d | Unchanged by any fix; zero `ImpulsiveOutlier` flags present |
| Pi Mensae c | P=6.266-6.270 d, depth 147-315 ppm, SNR 37.5-45.6 | P=6.2679 d | Required the merged-source fix; wrong (8.889 d) with the quality-flag-only screen |

### 12e. Lesson

A screen that trades recall for precision by trusting a single, specific
external label (here, one SPOC quality bit) inherits that label's blind
spots. The fix wasn't "make the periodicity check smarter" (Section 8
already solved the false-positive side) -- it was "don't throw away the
broader, noisier detector that has no such blind spot; use both, and let
the validated-signal check adjudicate between them."

## 13. Final packaged deliverable

The validated pipeline (contamination screen with merged candidate
sources, stability guard, and structural periodicity check; blind BLS
search; bootstrap significance) is packaged as `src/pipeline.py` and
`src/fits_reader.py`, with a command-line interface and a Python API. See
`README.md` for build and usage instructions. `examples/` contains the
exact script used to produce the Section 12d numbers, for reproducibility.
