#!/usr/bin/env python3
"""
Blind transit-detection pipeline with a hybrid contamination screen.

Combines:
  - A box-least-squares-style blind period search with bootstrap-calibrated
    significance (no fixed SNR cutoff -- the threshold is derived from the
    data itself).
  - A contamination screen that merges two candidate sources: TESS's own
    SPOC quality flags (high precision, but only catches sharp/impulsive
    events) and a blind self-calibrated Berlekamp-Welch scan (lower
    precision, but catches contamination outside SPOC's specific flag
    categories -- validated necessary in testing, see METHODS.md).
  - A structural periodicity check that decides whether to spare or excise
    each flagged candidate by referencing the pipeline's OWN bootstrap-
    validated best period, rather than an independent alias-prone search.

See METHODS.md for the full validation history, including two real bugs
found and fixed during testing (an unstable scoring statistic, and a null-
calibration mismatch that hid a real contaminant). Every default parameter
below was empirically validated, not guessed -- do not change them without
re-running the validation tests in tests/.

Usage:
    python pipeline.py --csv data.csv --fits sector1.fits sector2.fits \
        --period-min 0.5 --period-max 15.0 --label "My Target"

Or import as a module:
    from pipeline import run_pipeline
    result = run_pipeline(csv_path, fits_paths, period_range=(0.5, 15.0))
"""
import argparse
import sys
import time
import numpy as np
import pandas as pd

from fits_reader import read_lightcurve_fits, IMPULSIVE_OUTLIER_BIT

# ---- validated constants -- see METHODS.md Section 7-9 for how these were chosen ----
K = 4                          # local polynomial degree for Berlekamp-Welch
E_ORDER = 1                    # assumed errors per BW window
DETREND_WINDOW_DAYS = 0.5      # per-sector rolling-median detrend window
BLIND_BLOCK = 91               # ~3h at 2-min cadence, blind-scan block size
CLUSTER_GAP_DAYS = 1.0 / 24    # merge quality-flagged cadences within 1h into one region
CANDIDATE_MARGIN_DAYS = 1.5 / 24
EXCISE_MARGIN_DAYS = 3.0 / 24  # validated: must cover the real event's recovery tail
DEDUP_DAYS = 0.3               # candidates from both sources within this = same event
DEFAULT_N_PERIODS = 5000
DEFAULT_N_BOOT = 20


# =====================================================================
# core numerics
# =====================================================================
def berlekamp_welch(x, y, k, e):
    """Classical Berlekamp-Welch error-locating decoder, adapted for
    overdetermined real-valued windows. Returns (Q, E) coefficient arrays
    (low-to-high degree); P = Q/E is the recovered smooth polynomial."""
    n = len(x)
    n_q, n_e = k + e, e
    A = np.zeros((n, n_q + n_e))
    b = np.zeros(n)
    for i in range(n):
        A[i, :n_q] = [x[i] ** j for j in range(n_q)]
        A[i, n_q:] = [-y[i] * x[i] ** j for j in range(n_e)]
        b[i] = y[i] * x[i] ** e
    coeffs, *_ = np.linalg.lstsq(A, b, rcond=None)
    Q = coeffs[:n_q]
    E = np.concatenate([coeffs[n_q:], [1.0]])
    return Q, E


def window_score(t_win, y_win, gap_factor=3.0, cond_thresh=1e8):
    """Local anomaly score for one window, with a numerical stability
    guard. Validated finding: a window spanning a large internal gap (e.g.
    a sector boundary) or an ill-conditioned fit can make the Berlekamp-
    Welch polynomial division blow up to physically nonsensical values
    (>1,000,000 ppm was observed at exact sector-boundary blocks in
    testing). When that's detected, falls back to a simple robust
    (MAD-based) deviation statistic instead of trusting the fit.
    Returns (score, index_of_flagged_point, method_used)."""
    dt = np.diff(np.sort(t_win))
    median_dt = np.median(dt) if len(dt) else 0
    has_big_gap = median_dt > 0 and dt.max() > gap_factor * median_dt

    x0 = t_win - t_win.mean()
    A_check = np.vander(x0, K + E_ORDER, increasing=True)
    cond = np.linalg.cond(A_check) if len(x0) >= K + E_ORDER else np.inf
    ill_conditioned = not np.isfinite(cond) or cond > cond_thresh

    if has_big_gap or ill_conditioned:
        med = np.median(y_win)
        j = np.argmax(np.abs(y_win - med))
        score = min(np.abs(y_win[j] - med), 0.5)  # 0.5 = physical sanity cap (50%)
        return score, j, "robust_fallback"

    Q, E = berlekamp_welch(x0, y_win, k=K, e=E_ORDER)
    Qp = np.poly1d(Q[::-1])
    Ep = np.poly1d(E[::-1])
    Pp, _ = np.polydiv(Qp, Ep)
    pred = np.array([Pp(xi) for xi in x0])
    resid_pts = y_win - pred
    j = np.argmax(np.abs(resid_pts))
    score = min(np.abs(resid_pts[j]), 0.5)
    return score, j, "berlekamp_welch"


def rolling_median_detrend(t, f, window_days):
    n = len(t)
    out = np.empty(n)
    half = window_days / 2
    idx_lo = idx_hi = 0
    for i in range(n):
        while idx_lo < n and t[idx_lo] < t[i] - half:
            idx_lo += 1
        if idx_hi < i:
            idx_hi = i
        while idx_hi < n - 1 and t[idx_hi + 1] <= t[i] + half:
            idx_hi += 1
        out[i] = np.median(f[idx_lo:idx_hi + 1])
    return out


def bls_search(t, resid, sigma, period_grid, n_bins=100, widths=(2, 3, 4, 5, 6, 8, 10, 12, 15)):
    """Box-least-squares-style blind period search. Returns the highest-SNR
    (period, box width, box location) combination across the grid."""
    best = {"snr": -np.inf}
    total_sum = resid.sum()
    total_cnt = len(resid)
    for P in period_grid:
        ph = ((t - t[0]) / P) % 1.0
        bin_idx = np.minimum((ph * n_bins).astype(int), n_bins - 1)
        bin_sum = np.bincount(bin_idx, weights=resid, minlength=n_bins)
        bin_cnt = np.bincount(bin_idx, minlength=n_bins)
        for w in widths:
            ext_sum = np.concatenate([bin_sum, bin_sum[:w]])
            ext_cnt = np.concatenate([bin_cnt, bin_cnt[:w]])
            win_sum = np.convolve(ext_sum, np.ones(w), mode="valid")[:n_bins]
            win_cnt = np.convolve(ext_cnt, np.ones(w), mode="valid")[:n_bins]
            win_cnt_safe = np.where(win_cnt == 0, 1, win_cnt)
            out_cnt = total_cnt - win_cnt
            out_cnt_safe = np.where(out_cnt == 0, 1, out_cnt)
            depth = (total_sum - win_sum) / out_cnt_safe - win_sum / win_cnt_safe
            snr = depth * np.sqrt(win_cnt_safe) / sigma
            snr = np.where(win_cnt > 15, snr, -np.inf)
            j = np.argmax(snr)
            if snr[j] > best["snr"]:
                best = {"snr": snr[j], "period": P, "depth": depth[j], "n_in": win_cnt[j],
                        "width_bins": w, "phase_start_bin": j, "n_bins": n_bins}
    return best


def bootstrap_threshold(t, resid, sigma, period_grid, n_boot=DEFAULT_N_BOOT, seed=0, pctl=99):
    """Circular-shift bootstrap null for BLS significance -- the honest
    alternative to a fixed SNR cutoff. See METHODS.md Section 3."""
    rng = np.random.default_rng(seed)
    snrs = []
    for _ in range(n_boot):
        shift = rng.integers(1000, len(resid) - 1000)
        b = bls_search(t, np.roll(resid, shift), sigma, period_grid)
        snrs.append(b["snr"])
    return float(np.percentile(snrs, pctl))


# =====================================================================
# contamination screen: two candidate sources, merged
# =====================================================================
def quality_flag_candidates(t, resid, raw_t, quality, flag_bit=IMPULSIVE_OUTLIER_BIT):
    """Candidate source A: cluster SPOC quality-flagged cadences, localize
    the exact peak within each cluster via Berlekamp-Welch. High precision,
    but only sees whatever category `flag_bit` covers (default:
    ImpulsiveOutlier, i.e. sharp/impulsive events -- NOT slower multi-hour
    ramps; see METHODS.md for the real case this missed)."""
    flagged_idx = np.where((quality.astype(int) & flag_bit) != 0)[0]
    if len(flagged_idx) == 0:
        return []
    flagged_t = np.sort(raw_t[flagged_idx])
    clusters = [[flagged_t[0]]]
    for tt in flagged_t[1:]:
        if tt - clusters[-1][-1] <= CLUSTER_GAP_DAYS:
            clusters[-1].append(tt)
        else:
            clusters.append([tt])
    out = []
    for c in clusters:
        center = np.mean(c)
        m = np.abs(t - center) <= CANDIDATE_MARGIN_DAYS
        if m.sum() < K + 2 * E_ORDER + 5:
            continue
        tw, yw = t[m], resid[m]
        score, j, method = window_score(tw, yw)
        out.append({"source": "quality_flag", "region_t": center, "located_t": tw[j],
                     "score": score, "method": method})
    return out


def blind_scan_candidates(t, resid, sigma, n_null=150, seed=0):
    """Candidate source B: scan every block, self-calibrate a threshold from
    genuinely contiguous (rolled, not scattered-index) null blocks, flag
    anything that clears it. Lower precision than source A, but has no
    blind spot tied to a specific SPOC flag category -- validated necessary
    in testing (see METHODS.md, the t~1535 case)."""
    rng = np.random.default_rng(seed)
    n = len(resid)
    null_scores = np.empty(n_null)
    for i in range(n_null):
        shift = rng.integers(0, n)
        resid_shifted = np.roll(resid, shift)
        start = rng.integers(0, n - BLIND_BLOCK)
        s_, _, _ = window_score(t[start:start + BLIND_BLOCK], resid_shifted[start:start + BLIND_BLOCK])
        null_scores[i] = s_
    thresh = np.percentile(null_scores, 99)

    n_blocks = n // BLIND_BLOCK
    out = []
    for b in range(n_blocks):
        lo, hi = b * BLIND_BLOCK, (b + 1) * BLIND_BLOCK
        score, j, method = window_score(t[lo:hi], resid[lo:hi])
        if score >= thresh:
            out.append({"source": "blind_scan", "region_t": t[lo:hi][j], "located_t": t[lo:hi][j],
                         "score": score, "method": method})
    return out, thresh


def merge_candidates(cands_a, cands_b):
    merged = list(cands_a)
    for cb in cands_b:
        match = next((ca for ca in merged if abs(cb["located_t"] - ca["located_t"]) < DEDUP_DAYS), None)
        if match is None:
            merged.append(cb)
        else:
            match["source"] = match["source"] + "+blind_scan"
    return merged


def apply_periodicity_check(candidates, t, resid, sigma, period_grid, seed=0, n_boot=15):
    """Structural fix (validated over two failed alternatives -- see
    METHODS.md Section 8): decide SPARE vs EXCISE for each candidate by
    checking it against the pipeline's OWN bootstrap-validated best period,
    not an independent alias-prone periodicity search."""
    if not candidates:
        return candidates, None
    ref_best = bls_search(t, resid, sigma, period_grid)
    ref_thresh = bootstrap_threshold(t, resid, sigma, period_grid, n_boot=n_boot, seed=seed)
    ref_valid = ref_best["snr"] >= ref_thresh

    if ref_valid:
        ref_period = ref_best["period"]
        ref_dur = ref_period * ref_best["width_bins"] / ref_best["n_bins"]
        ref_phase_c = (ref_best["phase_start_bin"] + ref_best["width_bins"] / 2) / ref_best["n_bins"]
        ref_epoch = t[0] + ref_phase_c * ref_period
        for c in candidates:
            ph = ((c["located_t"] - ref_epoch) / ref_period + 0.5) % 1.0 - 0.5
            half_w = 1.5 * (ref_dur / ref_period) / 2
            c["verdict"] = "SPARE" if abs(ph) < half_w else "EXCISE"
    else:
        for c in candidates:
            c["verdict"] = "EXCISE"  # no validated signal to protect -- conservative default

    return candidates, {"ref_valid": ref_valid, "ref_snr": ref_best["snr"], "ref_thresh": ref_thresh,
                         "ref_period": ref_best.get("period")}


# =====================================================================
# top-level pipeline
# =====================================================================
def run_pipeline(csv_path, fits_paths, period_range, n_periods=DEFAULT_N_PERIODS,
                  n_boot=DEFAULT_N_BOOT, seed=0, verbose=True, label=None):
    """Run the full pipeline: detrend -> merged contamination screen ->
    structural periodicity check -> excise -> blind BLS search ->
    bootstrap significance. Returns a dict with the final result."""
    def log(msg):
        if verbose:
            print(msg)

    if label:
        log(f"\n{'='*72}\n{label}\n{'='*72}")

    df = pd.read_csv(csv_path)
    t = df.time_btjd.values
    f = df.pdcsap_flux_normalized.values
    sector = df.sector.values

    resid = np.empty(len(df))
    for s in sorted(df.sector.unique()):
        m = sector == s
        resid[m] = f[m] - rolling_median_detrend(t[m], f[m], DETREND_WINDOW_DAYS)
    sigma = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    log(f"rows={len(df)}  sigma={sigma*1e6:.1f}ppm")

    all_t, all_q = [], []
    for p in fits_paths:
        tr, qr, _ = read_lightcurve_fits(p)
        good = ~np.isnan(tr)
        all_t.append(tr[good])
        all_q.append(qr[good])
    raw_t = np.concatenate(all_t)
    raw_q = np.concatenate(all_q)
    order = np.argsort(raw_t)
    raw_t, raw_q = raw_t[order], raw_q[order]

    PERIOD_GRID = np.exp(np.linspace(np.log(period_range[0]), np.log(period_range[1]), n_periods))

    t0 = time.time()
    cand_a = quality_flag_candidates(t, resid, raw_t, raw_q)
    cand_b, blind_thresh = blind_scan_candidates(t, resid, sigma, seed=seed)
    candidates = merge_candidates(cand_a, cand_b)
    log(f"screen: {time.time()-t0:.0f}s | quality_flag={len(cand_a)} blind_scan={len(cand_b)} "
        f"(thresh {blind_thresh*1e6:.0f}ppm) -> {len(candidates)} merged candidates")

    candidates, ref_info = apply_periodicity_check(candidates, t, resid, sigma, PERIOD_GRID,
                                                     seed=seed, n_boot=max(6, n_boot // 3))
    if ref_info:
        log(f"reference: P={ref_info['ref_period']}  SNR={ref_info['ref_snr']:.1f} "
            f"vs {ref_info['ref_thresh']:.1f} -> "
            f"{'validated' if ref_info['ref_valid'] else 'not validated (conservative excise-all)'}")
    for c in candidates:
        log(f"  [{c['source']:>20}] t={c['region_t']:.4f} -> {c['located_t']:.4f}  "
            f"score={c['score']*1e6:.0f}ppm ({c['method']})  {c['verdict']}")

    keep = np.ones(len(t), dtype=bool)
    for c in candidates:
        if c["verdict"] == "EXCISE":
            keep &= np.abs(t - c["located_t"]) >= EXCISE_MARGIN_DAYS
    log(f"excised {(~keep).sum()} cadences -> {keep.sum()} remain")

    t_clean, resid_clean = t[keep], resid[keep]
    sigma_clean = 1.4826 * np.median(np.abs(resid_clean - np.median(resid_clean)))
    best = bls_search(t_clean, resid_clean, sigma_clean, PERIOD_GRID)
    duration_h = best["period"] * best["width_bins"] / best["n_bins"] * 24
    thresh = bootstrap_threshold(t_clean, resid_clean, sigma_clean, PERIOD_GRID,
                                  n_boot=n_boot, seed=seed + 100)
    verdict = "DETECTED" if best["snr"] >= thresh else "NOT SIGNIFICANT"

    log(f"FINAL: P={best['period']:.5f}d depth={best['depth']*1e6:.0f}ppm duration={duration_h:.2f}h "
        f"SNR={best['snr']:.1f} n_in={best['n_in']:.0f}")
    log(f"bootstrap threshold: {thresh:.1f} -> {verdict}")

    return {"period": best["period"], "depth": best["depth"], "duration_hours": duration_h,
            "snr": best["snr"], "n_in": best["n_in"], "threshold": thresh, "verdict": verdict,
            "candidates": candidates, "n_excised": int((~keep).sum())}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="cleaned light-curve CSV (columns: time_btjd, "
                                                    "pdcsap_flux_normalized, quality, sector)")
    ap.add_argument("--fits", required=True, nargs="+", help="one or more original SPOC FITS files "
                                                                "(for the QUALITY column)")
    ap.add_argument("--period-min", type=float, default=0.5, help="days")
    ap.add_argument("--period-max", type=float, default=15.0, help="days")
    ap.add_argument("--n-periods", type=int, default=DEFAULT_N_PERIODS)
    ap.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--label", default=None)
    args = ap.parse_args()

    result = run_pipeline(args.csv, args.fits, (args.period_min, args.period_max),
                           n_periods=args.n_periods, n_boot=args.n_boot, seed=args.seed,
                           label=args.label or args.csv)
    sys.exit(0 if result["verdict"] == "DETECTED" else 1)


if __name__ == "__main__":
    main()
