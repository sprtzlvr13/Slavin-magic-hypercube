import numpy as np
import pandas as pd
import time
import sys
sys.path.insert(0, '/home/claude')
from fits_reader import read_lightcurve_fits, IMPULSIVE_OUTLIER_BIT

K = 4
E_ORDER = 1
BLIND_BLOCK = 91          # ~3h, same scale as the quality-flag candidate window
CLUSTER_GAP_DAYS = 1.0/24
CANDIDATE_MARGIN_DAYS = 1.5/24
EXCISE_MARGIN_DAYS = 3.0/24     # validated -- must cover the real event's recovery tail, not just the peak
DEDUP_DAYS = 0.3                 # candidates from both sources within this are the same event

def berlekamp_welch(x, y, k, e):
    n = len(x)
    n_q, n_e = k + e, e
    A = np.zeros((n, n_q + n_e)); b = np.zeros(n)
    for i in range(n):
        A[i, :n_q] = [x[i]**j for j in range(n_q)]
        A[i, n_q:] = [-y[i] * x[i]**j for j in range(n_e)]
        b[i] = y[i] * x[i]**e
    coeffs, *_ = np.linalg.lstsq(A, b, rcond=None)
    Q = coeffs[:n_q]
    E = np.concatenate([coeffs[n_q:], [1.0]])
    return Q, E

def window_score(t_win, y_win, gap_factor=3.0, cond_thresh=1e8):
    """BW-based local anomaly score, with a stability guard: if the window
    spans an unusually large internal gap (e.g. a sector boundary) or the
    fit is numerically ill-conditioned, the polynomial division can blow up
    to physically nonsensical values (confirmed: >1,000,000 ppm at exactly
    the last cadence of Sectors 1 and 4, where fixed-size blocks span the
    inter-sector gap). In that case, fall back to a simple robust statistic
    (max deviation from the local median, MAD-scaled) instead of trusting
    the ill-conditioned algebraic fit."""
    dt = np.diff(np.sort(t_win))
    median_dt = np.median(dt) if len(dt) else 0
    has_big_gap = median_dt > 0 and dt.max() > gap_factor * median_dt

    x0 = t_win - t_win.mean()
    A_check = np.vander(x0, K + E_ORDER, increasing=True)  # same scale as the BW design matrix
    cond = np.linalg.cond(A_check) if len(x0) >= K + E_ORDER else np.inf
    ill_conditioned = not np.isfinite(cond) or cond > cond_thresh

    if has_big_gap or ill_conditioned:
        med = np.median(y_win)
        mad = 1.4826 * np.median(np.abs(y_win - med))
        j = np.argmax(np.abs(y_win - med))
        score = np.abs(y_win[j] - med)
        # physically sane cap: nothing real dims a star by more than ~50%
        score = min(score, 0.5)
        return score, j, "robust_fallback"

    Q, E = berlekamp_welch(x0, y_win, k=K, e=E_ORDER)
    Qp = np.poly1d(Q[::-1]); Ep = np.poly1d(E[::-1])
    Pp, _ = np.polydiv(Qp, Ep)
    pred = np.array([Pp(xi) for xi in x0])
    resid_pts = y_win - pred
    j = np.argmax(np.abs(resid_pts))
    score = np.abs(resid_pts[j])
    score = min(score, 0.5)  # same physical sanity cap, belt-and-suspenders
    return score, j, "berlekamp_welch"

def rolling_median_detrend(t, f, window_days):
    n = len(t); out = np.empty(n); half = window_days / 2
    idx_lo = idx_hi = 0
    for i in range(n):
        while idx_lo < n and t[idx_lo] < t[i] - half: idx_lo += 1
        if idx_hi < i: idx_hi = i
        while idx_hi < n - 1 and t[idx_hi + 1] <= t[i] + half: idx_hi += 1
        out[i] = np.median(f[idx_lo:idx_hi + 1])
    return out

def bls_search(t, resid, sigma, period_grid, n_bins=100, widths=(2,3,4,5,6,8,10,12,15)):
    best = {"snr": -np.inf}
    total_sum = resid.sum(); total_cnt = len(resid)
    for P in period_grid:
        ph = ((t - t[0]) / P) % 1.0
        bin_idx = np.minimum((ph * n_bins).astype(int), n_bins - 1)
        bin_sum = np.bincount(bin_idx, weights=resid, minlength=n_bins)
        bin_cnt = np.bincount(bin_idx, minlength=n_bins)
        for w in widths:
            ext_sum = np.concatenate([bin_sum, bin_sum[:w]])
            ext_cnt = np.concatenate([bin_cnt, bin_cnt[:w]])
            win_sum = np.convolve(ext_sum, np.ones(w), mode='valid')[:n_bins]
            win_cnt = np.convolve(ext_cnt, np.ones(w), mode='valid')[:n_bins]
            win_cnt_safe = np.where(win_cnt == 0, 1, win_cnt)
            depth = (total_sum-win_sum)/np.where((total_cnt-win_cnt)==0,1,total_cnt-win_cnt) - win_sum/win_cnt_safe
            snr = depth * np.sqrt(win_cnt_safe) / sigma
            snr = np.where(win_cnt > 15, snr, -np.inf)
            j = np.argmax(snr)
            if snr[j] > best["snr"]:
                best = {"snr": snr[j], "period": P, "depth": depth[j], "n_in": win_cnt[j],
                        "width_bins": w, "phase_start_bin": j, "n_bins": n_bins}
    return best

# =====================================================================
# candidate source A: quality-flag prior (high precision, category-limited)
# =====================================================================
def quality_flag_candidates(t, resid, raw_t, quality, flag_bit=IMPULSIVE_OUTLIER_BIT):
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
        if m.sum() < K + 2*E_ORDER + 5:
            continue
        tw, yw = t[m], resid[m]
        score, j, method = window_score(tw, yw)
        out.append({"source": "quality_flag", "region_t": center, "located_t": tw[j],
                     "score": score, "method": method})
    return out

# =====================================================================
# candidate source B: blind self-calibrated scan (catches what SPOC's
# specific flag categories don't -- the actual gap found in testing)
# =====================================================================
def blind_scan_candidates(t, resid, sigma, n_null=80, seed=0):
    """Null calibration must use genuinely CONTIGUOUS blocks -- scattered
    random indices across the full baseline would carry huge internal gaps
    and trip the stability guard on every trial, mismatching the real scan
    (which is contiguous and rarely gap-affected). Roll the residual by a
    random offset and take a real contiguous slice, so the null preserves
    the actual local time-spacing structure."""
    rng = np.random.default_rng(seed)
    null_scores = np.empty(n_null)
    n = len(resid)
    for i in range(n_null):
        shift = rng.integers(0, n)
        resid_shifted = np.roll(resid, shift)
        start = rng.integers(0, n - BLIND_BLOCK)
        s_, _, _ = window_score(t[start:start+BLIND_BLOCK], resid_shifted[start:start+BLIND_BLOCK])
        null_scores[i] = s_
    thresh = np.percentile(null_scores, 99)

    n_blocks = len(resid) // BLIND_BLOCK
    out = []
    for b in range(n_blocks):
        lo, hi = b*BLIND_BLOCK, (b+1)*BLIND_BLOCK
        score, j, method = window_score(t[lo:hi], resid[lo:hi])
        if score >= thresh:
            out.append({"source": "blind_scan", "region_t": t[lo:hi][j],
                        "located_t": t[lo:hi][j], "score": score, "method": method})
    return out, thresh

def merge_candidates(cands_a, cands_b):
    merged = list(cands_a)
    for cb in cands_b:
        if not any(abs(cb["located_t"] - ca["located_t"]) < DEDUP_DAYS for ca in merged):
            merged.append(cb)
        else:
            # same event found by both sources -- keep the higher-confidence score, tag both
            for ca in merged:
                if abs(cb["located_t"] - ca["located_t"]) < DEDUP_DAYS:
                    ca["source"] = ca["source"] + "+blind_scan"
    return merged

def run_corrected_pipeline(csv_path, fits_paths, label, period_range, n_periods=2000, seed=0):
    print(f"\n{'='*72}\n{label}\n{'='*72}")
    df = pd.read_csv(csv_path)
    t = df.time_btjd.values
    f = df.pdcsap_flux_normalized.values
    sector = df.sector.values
    resid = np.empty(len(df))
    for s in sorted(df.sector.unique()):
        m = sector == s
        resid[m] = f[m] - rolling_median_detrend(t[m], f[m], 0.5)
    sigma = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    print(f"rows={len(df)}  sigma={sigma*1e6:.1f}ppm")

    all_t, all_q = [], []
    for p in fits_paths:
        tr, qr, _ = read_lightcurve_fits(p)
        good = ~np.isnan(tr)
        all_t.append(tr[good]); all_q.append(qr[good])
    raw_t = np.concatenate(all_t); raw_q = np.concatenate(all_q)
    order = np.argsort(raw_t); raw_t, raw_q = raw_t[order], raw_q[order]

    t0 = time.time()
    cand_a = quality_flag_candidates(t, resid, raw_t, raw_q)
    cand_b, blind_thresh = blind_scan_candidates(t, resid, sigma, seed=seed)
    candidates = merge_candidates(cand_a, cand_b)
    print(f"screen time: {time.time()-t0:.0f}s | quality-flag: {len(cand_a)}  "
          f"blind-scan: {len(cand_b)} (thresh {blind_thresh*1e6:.0f}ppm)  "
          f"-> merged: {len(candidates)} unique candidates")

    PERIOD_GRID = np.exp(np.linspace(np.log(period_range[0]), np.log(period_range[1]), n_periods))

    # reference signal + bootstrap validation (v4 structural periodicity logic)
    if candidates:
        ref_best = bls_search(t, resid, sigma, PERIOD_GRID)
        rng = np.random.default_rng(seed+50)
        boot_snrs = [bls_search(t, np.roll(resid, rng.integers(1000, len(resid)-1000)), sigma, PERIOD_GRID)["snr"]
                     for _ in range(6)]
        ref_thresh = np.percentile(boot_snrs, 99)
        ref_valid = ref_best["snr"] >= ref_thresh
        print(f"reference signal: P={ref_best['period']:.5f}d SNR={ref_best['snr']:.1f} vs {ref_thresh:.1f} "
              f"-> {'validated' if ref_valid else 'not validated'}")

        for c in candidates:
            if ref_valid:
                ref_period = ref_best['period']
                ref_dur = ref_period * ref_best['width_bins'] / ref_best['n_bins']
                ref_phase_c = (ref_best['phase_start_bin'] + ref_best['width_bins']/2) / ref_best['n_bins']
                ref_epoch = t[0] + ref_phase_c * ref_period
                ph = ((c["located_t"] - ref_epoch) / ref_period + 0.5) % 1.0 - 0.5
                half_w = 1.5 * (ref_dur / ref_period) / 2
                c["verdict"] = "SPARE" if abs(ph) < half_w else "EXCISE"
            else:
                c["verdict"] = "EXCISE"

        for c in candidates:
            print(f"  [{c['source']:>20}] t={c['region_t']:.4f} located={c['located_t']:.4f} "
                  f"score={c['score']*1e6:.0f}ppm ({c.get('method','?')}) -> {c['verdict']}")

    keep = np.ones(len(t), dtype=bool)
    for c in candidates:
        if c.get("verdict") == "EXCISE":
            keep &= np.abs(t - c["located_t"]) >= EXCISE_MARGIN_DAYS
    print(f"excised {(~keep).sum()} cadences -> {keep.sum()} remain")

    t_clean, resid_clean = t[keep], resid[keep]
    sigma_clean = 1.4826 * np.median(np.abs(resid_clean - np.median(resid_clean)))
    best = bls_search(t_clean, resid_clean, sigma_clean, PERIOD_GRID)
    duration_h = best['period'] * best['width_bins'] / best['n_bins'] * 24
    print(f"FINAL: P={best['period']:.5f}d depth={best['depth']*1e6:.0f}ppm duration={duration_h:.2f}h "
          f"SNR={best['snr']:.1f} n_in={best['n_in']:.0f}")

    rng2 = np.random.default_rng(seed+100)
    boot2 = [bls_search(t_clean, np.roll(resid_clean, rng2.integers(1000, len(resid_clean)-1000)),
                         sigma_clean, PERIOD_GRID)["snr"] for _ in range(8)]
    thresh2 = np.percentile(boot2, 99)
    verdict = "DETECTED" if best["snr"] >= thresh2 else "NOT SIGNIFICANT"
    print(f"bootstrap threshold: {thresh2:.1f} -> {verdict}")
    return best, verdict

res_wasp = run_corrected_pipeline(
    '/home/claude/wasp18_extract/wasp18_tess_claude_package/cleaned_csv/wasp18_tess_sectors02_03_combined_clean.csv',
    ['/home/claude/wasp18_extract/wasp18_tess_claude_package/original_fits/tess_sector02_wasp18_lc.fits',
     '/home/claude/wasp18_extract/wasp18_tess_claude_package/original_fits/tess_sector03_wasp18_lc.fits'],
    "WASP-18 (expect: unchanged, P~0.9416d depth~9880-9979ppm)",
    period_range=(0.3, 3.0), seed=1)

res_pimen = run_corrected_pipeline(
    '/home/claude/pimen_extract/pi_men_tess_claude_package/cleaned_csv/pi_men_tess_sectors01_04_08_combined_clean.csv',
    ['/home/claude/pimen_extract/pi_men_tess_claude_package/original_fits/tess_sector01_pi_men_lc.fits',
     '/home/claude/pimen_extract/pi_men_tess_claude_package/original_fits/tess_sector04_pi_men_lc.fits',
     '/home/claude/pimen_extract/pi_men_tess_claude_package/original_fits/tess_sector08_pi_men_lc.fits'],
    "Pi Mensae (expect: FIXED, P~6.267d depth~175-315ppm, catching BOTH t~1421 and t~1535)",
    period_range=(0.5, 15.0), seed=2)

print("\n\n" + "="*72)
print("FINAL COMPARISON")
print("="*72)
print(f"WASP-18:   P={res_wasp[0]['period']:.5f}d depth={res_wasp[0]['depth']*1e6:.0f}ppm "
      f"SNR={res_wasp[0]['snr']:.1f} -> {res_wasp[1]}")
print(f"Pi Mensae: P={res_pimen[0]['period']:.5f}d depth={res_pimen[0]['depth']*1e6:.0f}ppm "
      f"SNR={res_pimen[0]['snr']:.1f} -> {res_pimen[1]}   [target: P~6.267d]")
