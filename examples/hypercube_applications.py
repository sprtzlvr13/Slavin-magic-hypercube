import numpy as np
import sys
sys.path.insert(0, '/home/claude')
from siamese_hypercube import siamese_hypercube_cascade_v2

n, D = 31, 4
print(f"building n={n} D={D} verified magic hypercube ({n**D:,} cells)...")
grid = siamese_hypercube_cascade_v2(n, D)  # values 1..923521
N = n ** D

# rank -> coordinate lookup (needed to place ranked tasks into the grid)
rank_to_coord = np.empty((N + 1, D), dtype=np.int32)
it = np.nditer(grid, flags=['multi_index'])
for val in it:
    rank_to_coord[int(val)] = it.multi_index

# =====================================================================
# USE CASE 1: provably fair sharding under a REAL skewed workload
# =====================================================================
print(f"\n{'='*70}\nUSE CASE 1: fair sharding\n{'='*70}")
rng = np.random.default_rng(0)
# realistic skewed task weights -- most tasks cheap, some very expensive
# (Pareto-like, common for real job/task-cost distributions)
weights = rng.pareto(a=1.5, size=N) + 0.1
ranks = np.argsort(np.argsort(weights)) + 1  # rank 1..N by weight, ties broken arbitrarily

N_WORKERS = n  # use one axis (31 values) as "worker id"

# --- method A: hypercube-line assignment ---
coords = rank_to_coord[ranks]  # coords[i] = hypercube coordinate of task i
worker_id_hc = coords[:, 0]    # group by axis 0's line membership
totals_hc = np.array([weights[worker_id_hc == w].sum() for w in range(N_WORKERS)])

# --- method B: naive contiguous rank chunking (a common real-world naive approach) ---
order = np.argsort(ranks)
chunk_size = N // N_WORKERS
worker_id_naive = np.minimum(np.arange(N)[np.argsort(order)] // chunk_size, N_WORKERS - 1)
totals_naive = np.array([weights[worker_id_naive == w].sum() for w in range(N_WORKERS)])

# --- method C: pure random assignment ---
worker_id_rand = rng.integers(0, N_WORKERS, size=N)
totals_rand = np.array([weights[worker_id_rand == w].sum() for w in range(N_WORKERS)])

for name, totals in [("hypercube-line", totals_hc), ("naive contiguous-chunk", totals_naive),
                      ("random", totals_rand)]:
    print(f"{name:>24}: mean={totals.mean():10.1f}  std={totals.std():10.1f}  "
          f"(std/mean={100*totals.std()/totals.mean():5.2f}%)  "
          f"min={totals.min():10.1f}  max={totals.max():10.1f}")

# =====================================================================
# USE CASE 2: corruption detection via line-sum checksums
# =====================================================================
print(f"\n{'='*70}\nUSE CASE 2: corruption detection\n{'='*70}")
M = n * (N + 1) // 2  # magic constant

def line_sums_axis0(g):
    return g.sum(axis=0)  # sum along axis 0 for every fixed (i1,...,iD-1) -- gives n^(D-1) line sums

expected = line_sums_axis0(grid)
n_lines = expected.size
print(f"checking {n_lines:,} lines along axis 0, each must sum to {M:,}")
assert np.all(expected == M), "sanity check failed"

def corrupt_and_detect(grid, n_corruptions, n_trials, rng, simultaneous_pairs=False):
    detected = 0
    for _ in range(n_trials):
        g = grid.copy()
        flat_idx = rng.choice(g.size, size=n_corruptions, replace=False)
        if simultaneous_pairs and n_corruptions == 2:
            # adversarial case: corrupt two cells in the SAME line with opposite
            # offsets that exactly cancel -- the known blind spot
            idx0 = np.unravel_index(flat_idx[0], g.shape)
            other_axes = idx0[1:]
            a1, a2 = rng.choice(n, size=2, replace=False)
            i1 = (a1,) + other_axes
            i2 = (a2,) + other_axes
            delta = int(rng.integers(1, 1000))
            g[i1] += delta
            g[i2] -= delta
        else:
            for fi in flat_idx:
                idx = np.unravel_index(fi, g.shape)
                g[idx] += int(rng.integers(-1000, 1000))
        sums = line_sums_axis0(g)
        if not np.all(sums == M):
            detected += 1
    return detected / n_trials

rng2 = np.random.default_rng(1)
print("\nrandom (non-adversarial) corruption detection rate:")
for n_corr in (1, 2, 3, 5):
    rate = corrupt_and_detect(grid, n_corr, 500, rng2)
    print(f"  {n_corr} random cell(s) corrupted: {rate*100:.1f}% detected")

print("\nADVERSARIAL blind spot -- two corruptions in the SAME line, exactly canceling:")
rate_adv = corrupt_and_detect(grid, 2, 500, rng2, simultaneous_pairs=True)
print(f"  detection rate: {rate_adv*100:.1f}% (this is the known, real blind spot)")
