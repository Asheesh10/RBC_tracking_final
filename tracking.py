"""
STAGE 3 — Multi-Cell Kymograph Peeling + RTS Kalman Smoother (Exact mp4_v5 Base + True Full-Duration Exit Fix)

Why `mp4_v5` (commit e638566) was almost perfect, and what this fixes:
- In `mp4_v5`, Cell #1 (Frames 1538 -> 1894, duration = 356 frames) tracked perfectly
  from the bottom-right entrance to the top-left exit.
- However, on Cell #2 and Cell #3, `trace_single_cell_dp_ridge()` had an early `break`
  at `t_rel > 220` (`if t_rel > 220 ... break`) and allowed `max_step = 24 px/frame`.
  That prematurely terminated Cell #2 and Cell #3 at `t_rel = 221` (Frames 2556 and 2596
  instead of their true exits at Frames 2730 and 2805), forcing the box to abandon the
  cell near the end and sprint at 11-17 px/frame to hit the premature `t_rel = 221` deadline!

What is fixed here:
1. Keeps the exact `e638566` (`mp4_v5`) entrance detection, 2D space-time peeling, and
   single-file capillary non-overlap gate.
2. Removes the premature `t_rel > 220: break` and searches the true physical passage
   duration window `t_rel in [315 .. 425]` frames, while capping `max_step` to physical
   cell speed (`<= 8 px/frame` in bore, `<= 4 px/frame` in flares) so the box NEVER leaves
   the cell or races ahead at the end!
"""
import os
import sys
import csv
import json
import shutil
import subprocess
import argparse
from pathlib import Path
from collections import defaultdict
from numpy.lib.stride_tricks import sliding_window_view

import cv2
import numpy as np

# ====================================================================
# PATH RESOLUTION
# ====================================================================
ROOT_DIR = Path(__file__).resolve().parent
VIDEO_FILENAME = "2A002 2_C001H001S0001-001.avi"


def _resolve_out_dir():
    local_out = ROOT_DIR / "output"
    legacy_out = Path(r"C:\Users\amara\rbc_tracking\output")
    if (local_out / "calibration.json").is_file() or not (legacy_out / "calibration.json").is_file():
        return str(local_out)
    return str(legacy_out)


def _resolve_video_path(calib_video_path):
    if calib_video_path and os.path.isfile(calib_video_path):
        return calib_video_path
    candidates = [
        ROOT_DIR / VIDEO_FILENAME,
        ROOT_DIR.parent / VIDEO_FILENAME,
        Path.home() / "Documents" / "IIsc" / VIDEO_FILENAME,
        Path.home() / "Desktop" / "rbctrack" / VIDEO_FILENAME,
        Path.home() / "Downloads" / VIDEO_FILENAME,
        Path(r"C:\Users\amara\Downloads") / VIDEO_FILENAME,
    ]
    for c in candidates:
        if c.is_file():
            return str(c.resolve())
    return calib_video_path


OUT = _resolve_out_dir()
CALIB_FILE = os.path.join(OUT, "calibration.json")
BG_FILE = os.path.join(OUT, "median_bg.npy")


def load_calibration():
    with open(CALIB_FILE, "r", encoding="utf-8") as f:
        return json.load(f)


def load_background():
    return np.load(BG_FILE)


# ====================================================================
# EXACT FULL-WIDTH INTERIOR LUMEN GEOMETRY (x = 0 .. W-1)
# ====================================================================

def _medfilt1d(arr, k=15):
    pad = k // 2
    padded = np.pad(arr, pad, mode="edge")
    return np.median(sliding_window_view(padded, k), axis=-1)


def build_exact_lumen_geometry(bg):
    h, w = bg.shape[:2]

    xs_fit = np.arange(120, w - 120, dtype=np.float64)
    ys_fit = np.array([bg[:, int(x)].argmin() for x in xs_fit], dtype=np.float64)
    for _ in range(4):
        A = np.polyfit(xs_fit, ys_fit, 1)
        r = ys_fit - np.polyval(A, xs_fit)
        keep = np.abs(r) < 2.0 * r.std()
        if keep.sum() < 20:
            break
        xs_fit, ys_fit = xs_fit[keep], ys_fit[keep]
    A = np.polyfit(xs_fit, ys_fit, 1)
    slope, intercept = float(A[0]), float(A[1])

    top_wall = np.zeros(w, dtype=np.float64)
    bot_wall = np.zeros(w, dtype=np.float64)

    for x in range(w):
        yc = slope * x + intercept
        max_search = 28 if (x < 115 or x > 1165) else 9
        y_lo = max(0, int(round(yc - max_search)))
        y_hi = min(h - 1, int(round(yc + max_search)))
        y_mid = np.clip(int(round(yc)), y_lo + 1, y_hi - 1)

        top_slice = bg[y_lo:y_mid, x]
        bot_slice = bg[y_mid + 1:y_hi + 1, x]

        top_wall[x] = y_lo + int(np.argmin(top_slice)) if top_slice.size > 0 else yc - 3.5
        bot_wall[x] = (y_mid + 1) + int(np.argmin(bot_slice)) if bot_slice.size > 0 else yc + 3.5

    top_smooth = _medfilt1d(top_wall, k=15)
    bot_smooth = _medfilt1d(bot_wall, k=15)

    lumen_mask = np.zeros((h, w), dtype=np.uint8)
    for x in range(w):
        yc = 0.5 * (top_smooth[x] + bot_smooth[x])
        peel = 4.5 if (x < 110 or x > 1170) else 1.2
        half_w = max(2.2, 0.5 * (bot_smooth[x] - top_smooth[x]) - peel)
        y0 = max(0, int(np.floor(yc - half_w)))
        y1 = min(h, int(np.ceil(yc + half_w)) + 1)
        lumen_mask[y0:y1, x] = 255

    return lumen_mask, top_smooth, bot_smooth, slope, intercept


# ====================================================================
# STAGE 3A: ITERATIVE 2D KYMOGRAPH PEELING + FULL-DURATION DP
# ====================================================================

_GAP_ENERGY_CACHE = {}


def find_entrance_via_components(roi_cache, bg_roi, mask_roi, gap_lo, gap_hi,
                                  skip_frames=90, x_entrance_min=1150,
                                  max_duration=435, min_duration=180, min_area=6,
                                  exclude_tracks=None, exclude_margin=20.0):
    """
    Gap-pass recovery of a faint cell the main kymograph pass missed (its entrance
    peak never cleared r_thresh). Returns (t_start_local, path_x) or (None, None).

    Why not a per-frame follower (the previous version of this function): measured on
    frames 15000-16500, the missed cell (entering ~15459) is too faint to follow frame by
    frame -- its per-frame column energy is ~3-6 against a bore noise p90 of 3.7. The old
    windowed-centroid follower (window [x-25, x+15]) therefore just drifted left through
    noise at ~5 px/frame (the window centre is x-5), ran ~140 px ahead of the real cell,
    then locked onto a debris particle stuck on the channel at x~731 for 106 frames
    (velocity 0.000). Two independent gap candidates both froze at x~731. The cell is only
    visible when evidence is integrated along time, so this version does what the main
    pass does -- a kymograph DP -- with these differences:

    1. Kymograph from the same motion*contrast energy as before, column-max over the
       lumen, minus each column's temporal median (removes static features such as the
       x~731 debris, whose median energy is ~22 vs ~13 elsewhere in the bore).
    2. Every already-accepted track (`exclude_tracks`, all of them) is masked +-
       `exclude_margin` px so this path can neither ride nor cross another cell (single
       file).
    3. Stall-tolerant DP: steps 0..9 px/frame in the bore (soft target 5), 0..4 in the
       flares (soft target 1.5); step 0 is penalised, not forbidden. No fixed passage
       duration. The path may START at any frame in [gap_lo+skip_frames, gap_hi) at
       x >= x_entrance_min, and must END at x <= 14. Evidence is biased by the bore noise
       p75 so noise is net-negative and the start lands where the ridge begins.
    4. Validation without DP selection bias: fit a quadratic x(t) to the bore part of the
       path and take the median 9-column-max energy along that fit. Accept only if it is
       >= 2.8x the same statistic on random straight bore-speed lines. Calibrated on
       15000-16500: random lines median 2.58 (max 5.02 of 300); spurious best path in a
       gap with no extra cell 5.47 (2.1x); the real missed cell 9.32 (3.6x); main-pass
       cells 12.8-31.1.

    `min_area` is unused (kept for signature compatibility).
    """
    debug = os.environ.get("RBC_DEBUG_ENTRANCES") == "1"
    n_frames = len(roi_cache)
    h_roi, w_img = mask_roi.shape

    search_lo = gap_lo + skip_frames
    search_hi = min(gap_hi, n_frames - 1)
    if search_hi <= search_lo:
        return None, None
    end = min(n_frames, search_hi + max_duration)

    # Per-frame column energy is cached across calls for the same roi_cache (every interval
    # of a window is searched now). The cache holds a reference to roi_cache itself, so the
    # identity check cannot be fooled by a recycled id() from a previous window.
    cache = _GAP_ENERGY_CACHE
    if cache.get("roi_cache") is not roi_cache or cache.get("mask") is not mask_roi:
        cache.clear()
        cache["roi_cache"], cache["mask"], cache["rows"] = roi_cache, mask_roi, {}
    rows = cache["rows"]

    def colmax_energy(t):
        t = max(1, t)
        hit = rows.get(t)
        if hit is not None:
            return hit
        cur = roi_cache[t].astype(np.float32)
        prev = roi_cache[t - 1].astype(np.float32)
        mot = np.abs(cur - prev)
        mot[mask_roi == 0] = 0.0
        absd = np.abs(cur - bg_roi)
        absd[mask_roi == 0] = 0.0
        e = mot * np.sqrt(np.maximum(1.0, absd))
        rows[t] = cv2.GaussianBlur(e, (5, 5), 0).max(axis=0)
        return rows[t]

    # Column baseline from a wider span so a cell lingering in one column of the search
    # window cannot dominate its own column's median. Every frame is used: sampling every
    # 2nd frame picked only odd or only even frames depending on the interval's start, and
    # this footage flickers between alternate frames -- measured 2026-11-16, the noise
    # reference jumped 1.91 <-> 3.42 when the interval start moved by one frame's parity.
    span_lo = max(1, search_lo - 300)
    span_hi = min(n_frames, end + 300)
    K_raw = np.stack([colmax_energy(t) for t in range(search_lo, end)]).astype(np.float64)
    base_rows = [K_raw[t - search_lo] if search_lo <= t < end else colmax_energy(t)
                 for t in range(span_lo, span_hi)]
    col_base = np.median(np.stack(base_rows), axis=0, keepdims=True)
    K = np.maximum(0.0, K_raw - col_base)
    T = K.shape[0]

    Kb = cv2.GaussianBlur(K.astype(np.float32), (11, 3), 0).astype(np.float64)
    bore_lo, bore_hi = 120, min(1170, w_img)
    noise_bias = float(np.percentile(Kb[:, bore_lo:bore_hi], 75))
    ev = Kb - noise_bias

    if exclude_tracks:
        for t0, xs_acc in exclude_tracks:
            for i, xa in enumerate(xs_acc):
                tl = t0 + i - search_lo
                if 0 <= tl < T:
                    xa_i = int(round(xa))
                    ev[tl, max(0, xa_i - int(exclude_margin)): xa_i + int(exclude_margin) + 1] = -50.0

    xs = np.arange(w_img)
    in_bore = (xs >= 110) & (xs <= 1180)
    step_target = np.where(in_bore, 5.0, 1.5)
    step_max = np.where(in_bore, 9, 4)
    LAM, STALL_PEN, NEG = 0.35, 1.5, -1e12
    X_EXIT = 14

    dp = np.full(w_img, NEG)
    back = np.full((T, w_img), -1, dtype=np.int32)   # -1 = path starts here
    start_ok = xs >= x_entrance_min
    best_val, best_t, best_x = NEG, None, None

    for t in range(T):
        new = np.full(w_img, NEG)
        bp = np.full(w_img, -1, dtype=np.int32)
        if t > 0:
            for s in range(0, 10):
                cand = np.full(w_img, NEG)
                cand[:w_img - s] = dp[s:]
                pen = LAM * (s - step_target) ** 2 + (STALL_PEN if s == 0 else 0.0)
                pen = np.where(s <= step_max, pen, 1e9)
                score = cand - pen
                better = score > new
                new = np.where(better, score, new)
                bp = np.where(better, xs + s, bp)
        if search_lo + t < search_hi:
            restart = start_ok & (new < 0.0)
            new = np.where(restart, 0.0, new)
            bp = np.where(restart, -1, bp)
        dp = new + ev[t]
        back[t] = bp
        v_end = float(dp[:X_EXIT + 1].max())
        if v_end > best_val:
            best_val, best_t, best_x = v_end, t, int(dp[:X_EXIT + 1].argmax())

    if best_t is None or best_val <= NEG / 2:
        return None, None

    x = best_x
    t = best_t
    path_rev = [x]
    while t > 0 and back[t, x] >= 0:
        x = int(back[t, x])
        t -= 1
        path_rev.append(x)
    path = path_rev[::-1]
    t_start = search_lo + t

    net_travel = float(path[0] - min(path))
    if len(path) < min_duration or net_travel < 80.0:
        if debug:
            print(f"      [DEBUG] gap-DP [{gap_lo},{gap_hi}] best path too short "
                  f"(len={len(path)}, travel={net_travel:.0f})", flush=True)
        return None, None

    # Validation on a 3-parameter fit (no DP selection bias), vs random-line noise.
    K9 = cv2.dilate(K.astype(np.float32), np.ones((1, 9), np.uint8))
    tl = np.arange(len(path)) + (t_start - search_lo)
    xp = np.asarray(path, dtype=np.float64)
    sel = (xp > 150) & (xp < 1100)
    if sel.sum() < 40:
        return None, None
    coef = np.polyfit(tl[sel], xp[sel], 2)
    xf = np.clip(np.round(np.polyval(coef, tl[sel])).astype(int), 0, w_img - 1)
    fit_med = float(np.median(K9[tl[sel], xf]))

    rng = np.random.default_rng(12345)
    null_meds = []
    for _ in range(300):
        v = rng.uniform(3.5, 6.5)
        x0 = rng.uniform(1000.0, 1100.0)
        n = int((x0 - 150.0) / v)
        if T - n <= 1:
            break
        t0 = int(rng.integers(0, T - n))
        tt = np.arange(t0, t0 + n)
        xx = np.clip(np.round(x0 - v * (tt - t0)).astype(int), 0, w_img - 1)
        null_meds.append(float(np.median(K9[tt, xx])))
    noise_ref = float(np.median(null_meds)) if null_meds else 1.0
    ratio = fit_med / max(noise_ref, 1e-3)

    if debug:
        print(f"      [DEBUG] gap-DP [{gap_lo},{gap_hi}] best path start={t_start} x0={path[0]} "
              f"len={len(path)} fit_med={fit_med:.2f} noise_ref={noise_ref:.2f} ratio={ratio:.2f}",
              flush=True)
    if ratio < 2.8:
        return None, None

    return t_start, [float(p) for p in path]


def extract_next_passage_from_kymograph(
    work_kymo,
    r_thresh,
    l_thresh,
    x_start_min=1225,
    x_exit_cutoff=14,
    r_band=None,
    entrance_thresh_override=None
):
    T_total, W = work_kymo.shape
    # Entrance is detected on a NARROW column band, not on the max across the whole flare.
    # Two cells crossing the flare at different columns peak at different times; taking
    # max() over cols 1215-1278 merges them into a single entrance. Callers pass a band
    # per column group so each crossing shows up separately.
    if r_band is None:
        r_band = (1215, 1278)
    r_raw = work_kymo[:, r_band[0]:r_band[1]].max(axis=1)
    l_raw = work_kymo[:, 5:105].max(axis=1)

    r_smooth = np.convolve(r_raw, np.ones(5) / 5.0, mode="same")
    l_smooth = np.convolve(l_raw, np.ones(7) / 7.0, mode="same")

    active_thresh = entrance_thresh_override if entrance_thresh_override is not None else r_thresh
    active_t = np.where(r_smooth[:T_total - 150] >= active_thresh)[0]
    if active_t.size == 0:
        return None, None

    t_first = int(active_t[0])
    win_end = min(T_total - 150, t_first + 20)
    t_pk = t_first + int(np.argmax(r_smooth[t_first:win_end]))

    edge_raw = work_kymo[:, 1245:1278].max(axis=1)
    t_enter = t_pk
    min_back = max(0, t_pk - 14)
    while t_enter > min_back and edge_raw[t_enter - 1] >= 0.45 * r_thresh:
        t_enter -= 1

    actual_t_enter, dp_xs = trace_single_cell_dp_ridge(
        work_kymo,
        t_enter,
        l_smooth,
        l_thresh,
        x_start_min=x_start_min,
        x_exit_cutoff=x_exit_cutoff
    )
    if actual_t_enter is None or not dp_xs:
        work_kymo[max(0, t_first - 2): min(T_total, t_pk + 8), 1215:1278] = 0.0
        return "retry", None

    return actual_t_enter, dp_xs


def trace_single_cell_dp_ridge(kymo_full, t_enter, l_smooth, l_thresh, x_start_min=1235, x_exit_cutoff=14,
                               exit_max_step=2, exit_target_ramp=None, exit_rule="mean_bonus"):
    """
    Traces a single cell's right-to-left trajectory ridge from `t_enter`:
    - STRICTLY FORBIDS `step = 0` (standing still in the flare/bore), which previously let
      Cell #2 sit in the right flare for 100 frames and swallow Cell #3's entrance!
    - Enforces physical zone velocities:
        * Flares (x < 110 or x > 1180): `step in [1, 2]` px/frame (target 1.2 px/frame)
        * Exit funnel override (re-trace only): `exit_max_step` / `exit_target_ramp`
    - Exit frame choice (`exit_rule`): "mean_bonus" (original: average DP score + left-exit
      bonus; main pass) or "noise_biased" (total score above the noise level; re-trace).
        * Narrow Bore (110 <= x <= 1180): `step in [3, 4, 5, 6, 7]` px/frame (target 5.1 px/frame)
      Because steps are strictly positive and physically bounded, crossing the 95-px right flare
      takes 48..80 frames, crossing the 1070-px bore takes 153..280 frames, and crossing the
      left flare takes 45..85 frames — the box can NEVER sit still and can NEVER sprint ahead!
    """
    T_total, W = kymo_full.shape
    max_duration = min(435, T_total - t_enter)
    if max_duration < 260:
        return None, None

    k_win = kymo_full[t_enter: t_enter + max_duration]
    k_smooth = cv2.GaussianBlur(k_win.astype(np.float32), (11, 3), 0)

    dp = np.full((max_duration, W), -1e12, dtype=np.float64)
    back = np.zeros((max_duration, W), dtype=np.int32)

    dp[0, x_start_min:] = k_smooth[0, x_start_min:] * 2.0

    # Exit funnel (x < 110) speed model. Defaults (max 2, target 1.2) reproduce the
    # original behaviour and are what the main extraction pass uses, so which cells get
    # accepted is unchanged. The neighbour-masked re-trace passes exit_max_step=5 and a
    # target easing 4.0 (throat, x=110) -> 1.5 (x=70) -> 1.2: measured 2026-11-15 on cell
    # #1 of frames 15000-16500, cells leave the bore at ~4 px/frame (x 110->72 in ~13
    # frames) then slow to ~1.2-1.4 px/frame, and the [1, 2] cap left the path ~30 px
    # behind at the throat. (Using the faster model in the main pass changed one 1450
    # candidate's path enough to fail raw_check and lose a real cell.) The ENTRY flare
    # always keeps [1, 2] -- that cap stops a cell sitting in the right flare and
    # swallowing the next one.
    if exit_target_ramp is None:
        exit_target_all = np.full(W, 1.2, dtype=np.float64)
    else:
        exit_target_all = np.interp(np.arange(W, dtype=np.float64), [0.0, 70.0, 110.0],
                                    [float(v) for v in exit_target_ramp])

    for t in range(1, max_duration):
        prev_row = dp[t - 1]
        best_val = np.full(W, -1e12, dtype=np.float64)
        best_prev = np.arange(W, dtype=np.int32)

        for step in range(1, 8):  # step=0 is FORBIDDEN; max_step=7
            limit = W - step
            cand_prev = prev_row[step:]
            x_coords = np.arange(limit)
            in_bore = (x_coords >= 110) & (x_coords <= 1180)
            in_exit = x_coords < 110
            valid_step = np.where(in_bore, (step >= 3) & (step <= 7),
                                  np.where(in_exit, (step >= 1) & (step <= exit_max_step), (step >= 1) & (step <= 2)))
            target_step = np.where(in_bore, 5.1, np.where(in_exit, exit_target_all[:limit], 1.2))
            step_penalty = np.where(valid_step, 0.55 * ((step - target_step) ** 2), 1e9)

            score = cand_prev - step_penalty
            better = score > best_val[:limit]
            best_val[:limit] = np.where(better, score, best_val[:limit])
            best_prev[:limit] = np.where(better, x_coords + step, best_prev[:limit])

        dp[t] = best_val + k_smooth[t]
        back[t] = best_prev

    min_exit_rel = min(330, max_duration - 15)
    best_t_rel = None
    best_exit_score = -1e12

    # "noise_biased": every path ending at t_rel has accumulated t_rel+1 frames, so
    # dp - noise * t_rel is the total evidence ABOVE noise. The path keeps extending while it
    # is on a real ridge and stops when the cell has left -- no pull toward short paths
    # (which averaging gives) and no pull toward a neighbour's exit (which the left-exit bonus
    # gives). Measured 2026-11-15: the original rule ended cell #3 of 15000-16500 ~47 frames
    # early (box sprinted ~70 px ahead to make an exit timed by cell #4 entering the funnel)
    # and left 1450 cell #2's box 40-50 px behind its cell. Noise = p75 of the smoothed bore
    # kymograph in this window; results were the same for p50 / p75 / p90.
    noise_level = (float(np.percentile(k_smooth[:, 110:min(1180, W)], 75))
                   if exit_rule == "noise_biased" else 0.0)

    for t_rel in range(min_exit_rel, max_duration):
        f_global = t_enter + t_rel
        val = float(dp[t_rel, :16].max())
        if val <= -1e8:
            continue
        if exit_rule == "noise_biased":
            total_score = val - noise_level * t_rel
        else:
            avg_dp_score = val / float(t_rel)
            left_bonus = 1.5 * float(l_smooth[f_global]) if f_global < len(l_smooth) else 0.0
            total_score = avg_dp_score + left_bonus
        if total_score > best_exit_score:
            best_exit_score = total_score
            best_t_rel = t_rel

    if best_t_rel is None:
        return None, None

    T_len = best_t_rel + 1
    end_x = int(np.argmax(dp[best_t_rel, :16]))
    path_x = [0] * T_len
    path_x[T_len - 1] = end_x
    for t in range(T_len - 1, 0, -1):
        path_x[t - 1] = int(back[t, path_x[t]])

    end_trim = T_len - 1
    for idx in range(T_len):
        if path_x[idx] <= x_exit_cutoff:
            end_trim = idx
            break

    if end_trim < 180:
        return None, None

    return t_enter, [float(x) for x in path_x[:end_trim + 1]]


def peel_ridge_from_kymograph(kymograph, t_enter, path_x, half_width_x=32, half_time_t=10):
    T_total, W = kymograph.shape
    for idx, x_val in enumerate(path_x):
        t_c = t_enter + idx
        t0 = max(0, t_c - half_time_t)
        t1 = min(T_total, t_c + half_time_t + 1)
        x0 = max(0, int(round(x_val)) - half_width_x)
        x1 = min(W, int(round(x_val)) + half_width_x + 1)
        kymograph[t0:t1, x0:x1] = 0.0
    # Clear only this cell's own entrance signature. The forward wipe used to reach
    # t_enter + 22 AND the full column range 1210:, which zeroed a cell entering up to
    # 22 frames behind this one (or at a different x in the same flare) before the loop
    # could detect it. +8 frames and a +-70px x-window around THIS path's own entrance
    # column clears its own signature without eating a nearby distinct entrance; any
    # accidental re-detection of the same cell is caught by is_duplicate_trajectory
    # (min_time_sep_frames=22).
    entrance_x0 = max(0, int(round(path_x[0])) - 70)
    kymograph[max(0, t_enter - 16): min(T_total, t_enter + 8), entrance_x0:] = 0.0


def is_duplicate_trajectory(cand_t0, cand_xs, accepted_list, min_spatial_sep_px=45.0, min_time_sep_frames=28):
    cand_map = {cand_t0 + i: x for i, x in enumerate(cand_xs)}

    for acc_t0, acc_xs in accepted_list:
        if abs(cand_t0 - acc_t0) < min_time_sep_frames:
            return True

        overlap_dists = []
        signs = []
        for i, ax in enumerate(acc_xs):
            f = acc_t0 + i
            if f in cand_map:
                cx = cand_map[f]
                if 110 <= ax <= 1180 or 110 <= cx <= 1180:
                    overlap_dists.append(abs(cx - ax))
                    if abs(cx - ax) > 5.0:
                        signs.append(1 if cx > ax else -1)

        if overlap_dists:
            if float(np.percentile(overlap_dists, 25)) < min_spatial_sep_px:
                return True
            if signs and len(set(signs)) > 1:
                return True

    return False


# ====================================================================
# STAGE 3B: SUB-PIXEL 2D MATCHED CENTROID + CONFIDENCE EXTRACTION
# ====================================================================

def extract_subpixel_observation(
    gray_roi,
    prev_gray_roi,
    bg_roi,
    mask_roi,
    guide_x,
    y_min_roi,
    top_wall,
    bot_wall
):
    h_roi, w_img = gray_roi.shape
    x_c = int(np.clip(round(guide_x), 0, w_img - 1))

    x0 = max(0, x_c - 20)
    x1 = min(w_img, x_c + 21)

    patch_gray = gray_roi[:, x0:x1]
    patch_bg = bg_roi[:, x0:x1]
    patch_mask = mask_roi[:, x0:x1]

    if prev_gray_roi is not None:
        patch_mot = np.abs(patch_gray - prev_gray_roi[:, x0:x1])
    else:
        patch_mot = np.abs(patch_gray - patch_bg)

    if x_c > 1180 or x_c < 110:
        pos_diff = np.maximum(0.0, patch_gray - patch_bg)
        energy = pos_diff * np.sqrt(np.maximum(1.0, patch_mot))
    else:
        abs_diff = np.abs(patch_gray - patch_bg)
        energy = abs_diff * np.sqrt(np.maximum(1.0, patch_mot))

    energy = cv2.GaussianBlur(energy.astype(np.float32), (5, 3), 0)
    energy[patch_mask == 0] = 0.0

    peak_e = float(energy.max())
    if peak_e >= 12.0:
        thresh = 0.48 * peak_e
        ys, xs = np.where(energy >= thresh)
    else:
        ys, xs = np.array([]), np.array([])

    center_y_wall = float(0.5 * (top_wall[x_c] + bot_wall[x_c]))

    if xs.size >= 3:
        weights = energy[ys, xs] ** 1.5
        w_sum = float(weights.sum())
        zx = float(x0 + np.sum(xs * weights) / w_sum)
        raw_zy = float(y_min_roi + np.sum(ys * weights) / w_sum)
        zy = 0.45 * raw_zy + 0.55 * center_y_wall
        bw_raw = int(xs.max() - xs.min() + 1)
        bh_raw = int(ys.max() - ys.min() + 1)
        conf = float(np.clip(peak_e / 80.0, 0.15, 1.0))
    else:
        zx = float(guide_x)
        zy = center_y_wall
        bw_raw, bh_raw = (14, 13) if (x_c > 1180 or x_c < 110) else (18, 6)
        conf = 0.05

    return {
        "zx": zx,
        "zy": zy,
        "bw_raw": bw_raw,
        "bh_raw": bh_raw,
        "conf": conf,
    }


# ====================================================================
# STAGE 3C: OPENCV KALMAN FILTER + FORWARD-BACKWARD RTS SMOOTHER
# ====================================================================

def kalman_rts_smooth_passage(observations, top_wall, bot_wall):
    N = len(observations)
    if N == 0:
        return []

    w_max = len(top_wall) - 1

    F = np.array(
        [[1, 0, 1, 0],
         [0, 1, 0, 1],
         [0, 0, 1, 0],
         [0, 0, 0, 1]],
        dtype=np.float64
    )
    H = np.array(
        [[1, 0, 0, 0],
         [0, 1, 0, 0]],
        dtype=np.float64
    )
    Q = np.diag([0.35, 0.08, 0.25, 0.02]).astype(np.float64)

    kf = cv2.KalmanFilter(4, 2)
    kf.transitionMatrix = F.astype(np.float32)
    kf.measurementMatrix = H.astype(np.float32)
    kf.processNoiseCov = Q.astype(np.float32)

    x_pred_list = np.zeros((N, 4), dtype=np.float64)
    P_pred_list = np.zeros((N, 4, 4), dtype=np.float64)
    x_filt_list = np.zeros((N, 4), dtype=np.float64)
    P_filt_list = np.zeros((N, 4, 4), dtype=np.float64)

    z0 = observations[0]
    vx0 = (observations[min(N - 1, 5)]["zx"] - z0["zx"]) / float(max(1, min(N - 1, 5)))
    vy0 = (observations[min(N - 1, 5)]["zy"] - z0["zy"]) / float(max(1, min(N - 1, 5)))
    x_curr = np.array([z0["zx"], z0["zy"], vx0, vy0], dtype=np.float64)
    P_curr = np.diag([4.0, 1.5, 4.0, 0.5]).astype(np.float64)

    kf.statePost = x_curr.reshape(4, 1).astype(np.float32)
    kf.errorCovPost = P_curr.astype(np.float32)

    for t in range(N):
        obs = observations[t]
        if t == 0:
            x_pred = x_curr.copy()
            P_pred = P_curr.copy()
        else:
            kf.predict()
            x_pred = F @ x_curr
            P_pred = F @ P_curr @ F.T + Q

        x_pred_list[t] = x_pred
        P_pred_list[t] = P_pred

        conf = max(0.05, obs["conf"])
        sigma_z = 2.2 / conf
        R_t = np.diag([sigma_z ** 2, (sigma_z * 0.5) ** 2])

        z_t = np.array([obs["zx"], obs["zy"]], dtype=np.float64)
        y_innov = z_t - (H @ x_pred)
        S = H @ P_pred @ H.T + R_t
        K = P_pred @ H.T @ np.linalg.inv(S)

        x_curr = x_pred + K @ y_innov
        P_curr = (np.eye(4) - K @ H) @ P_pred

        kf.measurementNoiseCov = R_t.astype(np.float32)
        kf.correct(z_t.reshape(2, 1).astype(np.float32))

        x_filt_list[t] = x_curr
        P_filt_list[t] = P_curr

    x_smooth = np.zeros((N, 4), dtype=np.float64)
    P_smooth = np.zeros((N, 4, 4), dtype=np.float64)
    x_smooth[-1] = x_filt_list[-1]
    P_smooth[-1] = P_filt_list[-1]

    for t in range(N - 2, -1, -1):
        P_pred_next = P_pred_list[t + 1]
        C_t = P_filt_list[t] @ F.T @ np.linalg.inv(P_pred_next + 1e-6 * np.eye(4))
        x_smooth[t] = x_filt_list[t] + C_t @ (x_smooth[t + 1] - x_pred_list[t + 1])
        P_smooth[t] = P_filt_list[t] + C_t @ (P_smooth[t + 1] - P_pred_next) @ C_t.T

    smoothed_track = []
    areas_so_far = []

    for t in range(N):
        obs = observations[t]
        cx = float(np.clip(x_smooth[t, 0], 0.0, w_max))
        if t > 0 and cx > smoothed_track[-1]["cx"]:
            cx = smoothed_track[-1]["cx"]

        x_idx = int(np.clip(round(cx), 0, w_max))
        wt = float(top_wall[x_idx])
        wb = float(bot_wall[x_idx])
        yc_wall = 0.5 * (wt + wb)
        # Clip the wall separation: at the extreme columns the wall search fails, `room`
        # explodes and cy can jump ~178 px in one frame (wiggly trail, haywire funnel ends).
        # Geometry only -- does not touch detection.
        lumen_h = float(np.clip(wb - wt, 4.0, 40.0))

        # y freedom scaled by the room the cell has; smooth ramp instead of a hard
        # pin/unpin, which put a kink in the path at the constriction boundary.
        room = max(0.0, 0.5 * lumen_h - 3.0)
        edge_taper = min(1.0, min(x_idx / 25.0, (w_max - x_idx) / 25.0))
        alpha = float(np.clip(room / 8.0, 0.0, 1.0)) * edge_taper
        dev = float(np.clip(x_smooth[t, 1] - yc_wall, -room, room)) if room > 0.0 else 0.0
        cy = float(yc_wall + alpha * dev)

        bw_raw, bh_raw = obs["bw_raw"], obs["bh_raw"]
        if x_idx > 1180 or x_idx < 110:
            bw = int(np.clip(max(bw_raw, 14), 13, 22))
            bh = int(np.clip(max(bh_raw, 13), 12, 20))
        else:
            bw = int(np.clip(max(bw_raw, 17), 15, 26))
            bh = int(np.clip(min(max(bh_raw, int(round(lumen_h))), int(round(lumen_h)) + 2), 5, 9))

        bx = int(round(cx - 0.5 * bw))
        by = int(round(cy - 0.5 * bh))

        a = max(1.0, 0.5 * bw)
        b = max(1.0, 0.5 * bh)
        area = float(np.pi * a * b)
        perim = float(np.pi * (3.0 * (a + b) - np.sqrt((3.0 * a + b) * (a + 3.0 * b))))
        circ = float(np.clip(4.0 * np.pi * area / (perim * perim), 0.15, 1.0))
        elong = float(max(bw, bh)) / float(max(1, min(bw, bh)))

        if t == 0:
            vel = float(np.hypot(x_smooth[0, 2], x_smooth[0, 3]))
        else:
            prev = smoothed_track[-1]
            vel = float(np.hypot(cx - prev["cx"], cy - prev["cy"]))

        areas_so_far.append(area)
        adi_run = float(np.min(areas_so_far) / max(1.0, np.mean(areas_so_far)))

        smoothed_track.append({
            "frame": obs["frame"],
            "track_id": obs["track_id"],
            "cx": cx,
            "cy": cy,
            "bbox": (bx, by, bw, bh),
            "area": area,
            "circularity": circ,
            "elongation": elong,
            "velocity": vel,
            "ADI_running": adi_run,
        })

    return smoothed_track


# ====================================================================
# FULL HYBRID PIPELINE
# ====================================================================

def detect_and_track_passages(
    video_path,
    start_frame,
    frames_to_run,
    background,
    lumen_mask,
    top_wall,
    bot_wall,
    roi_y_bounds
):
    y_min, y_max = roi_y_bounds
    bg_roi = background[y_min:y_max, :]
    mask_roi = lumen_mask[y_min:y_max, :]
    w_img = background.shape[1]

    cap = cv2.VideoCapture(video_path)
    if start_frame > 0:
        cap.set(cv2.CAP_PROP_POS_FRAMES, start_frame)

    kymograph = np.zeros((frames_to_run, w_img), dtype=np.float32)
    roi_cache = []

    prev_gray_roi = None
    lag3_buffer = []

    for offset in range(frames_to_run):
        ok, frame = cap.read()
        if not ok:
            kymograph = kymograph[:offset]
            break

        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) if frame.ndim == 3 else frame
        roi_u8 = gray[y_min:y_max, :].copy()
        roi_cache.append(roi_u8)
        gray_roi = roi_u8.astype(np.float32)

        pos_diff = np.maximum(0.0, gray_roi - bg_roi)
        pos_diff[mask_roi == 0] = 0.0

        abs_diff = np.abs(gray_roi - bg_roi)
        abs_diff[mask_roi == 0] = 0.0

        if len(lag3_buffer) >= 3:
            mot_diff = np.abs(gray_roi - lag3_buffer[0])
            mot_diff[mask_roi == 0] = 0.0
        elif prev_gray_roi is not None:
            mot_diff = np.abs(gray_roi - prev_gray_roi)
            mot_diff[mask_roi == 0] = 0.0
        else:
            mot_diff = np.zeros_like(abs_diff)

        prev_gray_roi = gray_roi
        lag3_buffer.append(gray_roi)
        if len(lag3_buffer) > 4:
            lag3_buffer.pop(0)

        col_abs = abs_diff.max(axis=0)
        col_pos = pos_diff.max(axis=0)
        col_mot = mot_diff.max(axis=0)

        col_score = np.sqrt(col_abs * np.maximum(0.0, col_mot - 4.0))
        col_score[:110] = np.sqrt(col_pos[:110] * np.maximum(0.0, col_mot[:110] - 4.0)) * 1.4
        col_score[1180:] = np.sqrt(col_pos[1180:] * np.maximum(0.0, col_mot[1180:] - 4.0)) * 1.4

        kymograph[offset] = col_score

        if (offset + 1) % 500 == 0:
            print(f"      [Stage 3A] Kymograph built for {offset + 1}/{frames_to_run} frames...", flush=True)

    cap.release()

    actual_frames = kymograph.shape[0]
    if actual_frames < 20:
        return {}

    col_baseline = np.median(kymograph, axis=0, keepdims=True)
    kymograph = np.maximum(0.0, kymograph - col_baseline)

    # Boost narrow-bore columns (110 <= x <= 1180) by 2.5x so the squashed cell inside the
    # 4-px bore has the same kymograph ridge amplitude as in the 35-px wide flares!
    kymograph[:, 110:1180] *= 2.5

    right_entrance_signal = kymograph[:, 1215:1278].max(axis=1)
    left_exit_signal = kymograph[:, 5:105].max(axis=1)

    r_smooth = np.convolve(right_entrance_signal, np.ones(7) / 7.0, mode="same")
    l_smooth = np.convolve(left_exit_signal, np.ones(7) / 7.0, mode="same")

    r_med = float(np.median(r_smooth))
    r_mad = float(np.median(np.abs(r_smooth - r_med))) + 1e-3
    r_thresh = max(8.0, r_med + 3.2 * r_mad)

    l_med = float(np.median(l_smooth))
    l_mad = float(np.median(np.abs(l_smooth - l_med))) + 1e-3
    l_thresh = max(8.0, l_med + 3.0 * l_mad)

    if os.environ.get("RBC_DEBUG_ENTRANCES") == "1":
        np.save(os.path.join(OUT, "stage3", "_debug_r_smooth.npy"), r_smooth)
        print(f"      [DEBUG] r_thresh={r_thresh:.2f} r_smooth saved, len={len(r_smooth)}", flush=True)

    work_kymo = kymograph.copy()
    raw_extracted = []
    max_cells = 40

    # Single wide entrance band. A two-narrow-band variant was tried (to separate cells
    # crossing the flare at different columns) and LOST entrances: running the bands
    # sequentially against one shared work_kymo let the first band's peels destroy the
    # evidence the others needed -- 5 cells became 3. Entrance detection stays wide; the
    # improvement that was verified to help is the narrower time peel below.
    DEBUG_ENTRANCES = os.environ.get("RBC_DEBUG_ENTRANCES") == "1"
    for _ in range(max_cells):
        actual_t_enter, dp_xs = extract_next_passage_from_kymograph(
            work_kymo, r_thresh, l_thresh, x_start_min=1225, x_exit_cutoff=14
        )
        if actual_t_enter is None:
            break
        if actual_t_enter == "retry":
            if DEBUG_ENTRANCES:
                print(f"      [DEBUG] retry (no valid dp ridge)", flush=True)
            continue

        if DEBUG_ENTRANCES:
            vel = np.diff(np.array(dp_xs))
            jump_idx = int(np.argmax(np.abs(vel))) if len(vel) else -1
            print(f"      [DEBUG] candidate t_enter={actual_t_enter} (frame {start_frame+actual_t_enter}) "
                  f"len={len(dp_xs)} x0={dp_xs[0]:.0f} maxstep={np.abs(vel).max() if len(vel) else 0:.1f}"
                  f"@rel{jump_idx}", flush=True)

        # half_time_t 10 -> 4: a cell entering close behind the previous one was being
        # erased by the peel before the loop could see it. (half_width_x stays 34; the
        # 18-px trial lost a cell to residual merge.)
        peel_ridge_from_kymograph(work_kymo, actual_t_enter, dp_xs,
                                  half_width_x=34, half_time_t=4)

        raw_check = [
            kymograph[actual_t_enter + idx, max(0, int(round(x)) - 4): min(w_img, int(round(x)) + 5)].max()
            for idx, x in enumerate(dp_xs)
            if 0 <= actual_t_enter + idx < actual_frames
        ]
        if not raw_check or float(np.median(raw_check)) < 0.55 * r_thresh:
            if DEBUG_ENTRANCES:
                print(f"      [DEBUG] REJECTED t_enter={actual_t_enter}: weak raw_check "
                      f"({float(np.median(raw_check)) if raw_check else 0:.2f} < {0.55*r_thresh:.2f})", flush=True)
            continue

        if is_duplicate_trajectory(actual_t_enter, dp_xs, raw_extracted,
                                   min_spatial_sep_px=45.0, min_time_sep_frames=22):
            if DEBUG_ENTRANCES:
                print(f"      [DEBUG] REJECTED t_enter={actual_t_enter}: duplicate of existing track", flush=True)
            continue

        raw_extracted.append((actual_t_enter, dp_xs))

    # Second pass: weak/low-contrast entrances (e.g. partly masked by a neighboring
    # cell's decaying signal) never cross r_thresh and are invisible to the first loop.
    # Relax ONLY the entrance-peak threshold; every other gate (full 180+ frame DP
    # trace, raw_check strength at the ORIGINAL r_thresh, duplicate-trajectory check)
    # stays as strict as before, so this cannot admit noise -- it can only recover a
    # real, physically-continuous passage that happened to be dim.
    weak_thresh = max(8.0, 0.55 * r_thresh)
    # Search EVERY interval between accepted entrances, plus before the first and after the
    # last (2026-11-16). The earlier rule searched only gaps longer than 1.25x the median
    # spacing (written for the old per-frame follower, which turned a low-threshold global
    # rescan into 11 false recoveries). That made results depend on where a window starts:
    # with the window decoded from 14550 the median spacing rose to 213 frames and the faint
    # 15459 cell's 168-frame gap was never searched, although the same gap qualified when
    # the window started at 15000. The gap DP masks every accepted track and must pass the
    # fit-ratio test (>= 2.8x random-line noise), which is what rejects noise now. A
    # recovered cell splits its interval and both halves are searched again, so more than
    # one faint cell per interval can be found.
    _ts = sorted(t for t, _ in raw_extracted)
    _bounds = [-90] + _ts + [actual_frames - 260]
    gap_candidates = [(a, b) for a, b in zip(_bounds[:-1], _bounds[1:]) if b - (a + 90) > 10]
    MAX_GAP_SEARCHES = 60

    gap_recovered_t0s = set()
    n_gap_searches = 0
    for gap_lo, gap_hi in gap_candidates:
        n_gap_searches += 1
        if n_gap_searches > MAX_GAP_SEARCHES:
            break
        actual_t_enter, dp_xs = find_entrance_via_components(
            roi_cache, bg_roi, mask_roi, gap_lo, gap_hi, skip_frames=90,
            exclude_tracks=raw_extracted
        )
        if actual_t_enter is None:
            if DEBUG_ENTRANCES:
                print(f"      [DEBUG] gap-search [{gap_lo},{gap_hi}] no component found", flush=True)
            continue

        if DEBUG_ENTRANCES:
            floor_dbg = {}
            if raw_extracted:
                for t0x, xsx in raw_extracted:
                    if t0x > gap_lo + 90:
                        continue
                    for ii, xx in enumerate(xsx):
                        ff = t0x + ii
                        floor_dbg[ff] = xx if ff not in floor_dbg else min(floor_dbg[ff], xx)
            sample_f = actual_t_enter + 200
            print(f"      [DEBUG] floor at frame {start_frame+sample_f}: "
                  f"{floor_dbg.get(sample_f)} vs path_x there: "
                  f"{dp_xs[200] if len(dp_xs)>200 else 'N/A'}", flush=True)
            print(f"      [DEBUG] gap-search [{gap_lo},{gap_hi}] candidate "
                  f"t_enter={actual_t_enter} (frame {start_frame+actual_t_enter}) "
                  f"len={len(dp_xs)} exclude_n={len(raw_extracted)} "
                  f"exclude_t0s={[x for x,_ in raw_extracted]}", flush=True)

        raw_check = [
            kymograph[actual_t_enter + idx, max(0, int(round(x)) - 4): min(w_img, int(round(x)) + 5)].max()
            for idx, x in enumerate(dp_xs)
            if 0 <= actual_t_enter + idx < actual_frames
        ]
        weak_raw_floor = 0.55 * weak_thresh
        if not raw_check or float(np.median(raw_check)) < weak_raw_floor:
            if DEBUG_ENTRANCES:
                print(f"      [DEBUG] gap REJECTED t_enter={actual_t_enter}: weak raw_check "
                      f"({float(np.median(raw_check)) if raw_check else 0:.2f} < {weak_raw_floor:.2f})", flush=True)
            continue

        # This candidate must obey single-file ordering with every already-accepted
        # track, AND must be abandoned the moment it gets too close to distinguish
        # from one (within roughly one cell diameter, ~20px) -- past that point the
        # windowed energy-centroid search cannot tell two blobs apart and silently
        # locks onto the OTHER track's signal instead of measuring this cell
        # independently (confirmed on real data: candidate's path became pixel-
        # identical to an accepted track's path for the rest of the passage once
        # separation dropped below ~8px). A pure order-flip check (candidate drifts
        # fully past/behind) is not strict enough since convergence to near-zero
        # separation happens several frames before any flip would register --
        # truncate at first close approach, not first actual overtake.
        MERGE_MARGIN_PX = 20.0
        cmap_trim = {actual_t_enter + i: x for i, x in enumerate(dp_xs)}
        cross_frame = None
        for at0, axs in raw_extracted:
            for i, ax in enumerate(axs):
                f = at0 + i
                if f not in cmap_trim:
                    continue
                cx = cmap_trim[f]
                if not (110 <= ax <= 1180 or 110 <= cx <= 1180):
                    continue
                if abs(cx - ax) < MERGE_MARGIN_PX:
                    if cross_frame is None or f < cross_frame:
                        cross_frame = f
        if cross_frame is not None:
            keep_len = cross_frame - actual_t_enter
            if DEBUG_ENTRANCES:
                culprit = None
                for at0, axs in raw_extracted:
                    for i, ax in enumerate(axs):
                        f = at0 + i
                        if f == cross_frame and f in cmap_trim:
                            culprit = (at0, ax, cmap_trim[f])
                print(f"      [DEBUG] cross_frame={cross_frame} culprit_track_t0={culprit}", flush=True)
            if keep_len < min(180, len(dp_xs)):
                if DEBUG_ENTRANCES:
                    print(f"      [DEBUG] gap REJECTED t_enter={actual_t_enter}: "
                          f"converges with existing track too early (keep_len={keep_len})", flush=True)
                continue
            dp_xs = dp_xs[:keep_len]
            if DEBUG_ENTRANCES:
                print(f"      [DEBUG] gap-search t_enter={actual_t_enter}: truncated to "
                      f"{keep_len} frames at convergence with existing track "
                      f"(can no longer be distinguished from it)", flush=True)

        # NOTE: the generic is_duplicate_trajectory() check (position-only, no entry-
        # order awareness) is intentionally NOT used here. For a gap-recovered queued
        # cell, natural proximity to the cell ahead is expected and not a duplicate --
        # the entry-order-aware crossing-trim above already enforces the real physical
        # constraint (never overtake) and is strictly more correct for this case. Still
        # guard against two near-simultaneous spurious detections of the same signal.
        # NOTE: only a FIRST-FRAMES spatial-separation check is used here (not the
        # generic is_duplicate_trajectory, which has no entry-order awareness and
        # wrongly flags a genuinely queued cell as a duplicate once it gets close to
        # the cell ahead near a later crossing). A TRUE ghost duplicate -- the search
        # window merely riding the exclude-floor margin of an already-accepted track,
        # never showing independent motion -- stays pinned at ~15px (the floor margin)
        # from frame 1, confirmed on real data. A genuine distinct cell starts >100px
        # away even when it later catches up. So only the first ~20 overlapping frames
        # are checked, not the whole track (which would also catch legitimate later
        # convergence).
        cmap_dup = {actual_t_enter + i: x for i, x in enumerate(dp_xs)}
        is_ghost_dup = False
        for at0, axs in raw_extracted:
            amap = {at0 + i: x for i, x in enumerate(axs)}
            overlap_frames = sorted(set(cmap_dup) & set(amap))[:20]
            if len(overlap_frames) < 10:
                continue
            seps = [abs(cmap_dup[f] - amap[f]) for f in overlap_frames]
            if float(np.median(seps)) < 35.0:
                is_ghost_dup = True
                break
        if is_ghost_dup:
            if DEBUG_ENTRANCES:
                print(f"      [DEBUG] gap REJECTED t_enter={actual_t_enter}: "
                      f"ghost duplicate (co-tracks an accepted cell's signal)", flush=True)
            continue

        if any(abs(actual_t_enter - at0) < 22 for at0, _ in raw_extracted):
            if DEBUG_ENTRANCES:
                print(f"      [DEBUG] gap REJECTED t_enter={actual_t_enter}: "
                      f"too close in time to an existing track entrance", flush=True)
            continue

        if DEBUG_ENTRANCES:
            print(f"      [DEBUG] GAP-RECOVERED t_enter={actual_t_enter} "
                  f"(frame {start_frame+actual_t_enter})", flush=True)
        raw_extracted.append((actual_t_enter, dp_xs))
        gap_recovered_t0s.add(actual_t_enter)
        for sub in ((gap_lo, actual_t_enter), (actual_t_enter, gap_hi)):
            if sub[1] - (sub[0] + 90) > 10:
                gap_candidates.append(sub)

    # Neighbour-masked re-trace of every main-pass track. Each passage above was traced on
    # its own with the original exit rule, so its path can (a) be pulled onto a brighter
    # neighbour's ridge and (b) end at the wrong frame. Measured 2026-11-15:
    #   - 15000-16500 cell #1 rode cell #2's ridge out of the exit funnel (box 30-36 px
    #     behind the visible cell);
    #   - 15000-16500 cell #3 ended ~47 frames early, its box sprinting ~70 px ahead;
    #   - 1450-2950 cell #2's box trailed its cell by 40-50 px and outlived it by ~40 frames.
    # Each track is re-traced on the unpeeled kymograph with all OTHER tracks blanked
    # (+-14 px, +-4 frames), the faster exit-funnel speed model, and the noise-biased exit
    # rule (see trace_single_cell_dp_ridge). Two rounds so corrections propagate between
    # neighbours. The main pass above still uses the original model, so WHICH cells are
    # accepted does not change. Gap-recovered tracks are never re-traced (different DP model:
    # stalls allowed, no minimum duration) but still blank their band. A failed re-trace
    # keeps the original path.
    RETRACE_HALF_X, RETRACE_HALF_T = 14, 4
    for _ in range(2):
        for idx, (t0_i, xs_i) in enumerate(raw_extracted):
            if t0_i in gap_recovered_t0s:
                continue
            k_masked = kymograph.copy()
            for jdx, (t0_j, xs_j) in enumerate(raw_extracted):
                if jdx == idx:
                    continue
                for s_j, xj in enumerate(xs_j):
                    tc = t0_j + s_j
                    xc = int(round(xj))
                    k_masked[max(0, tc - RETRACE_HALF_T):min(actual_frames, tc + RETRACE_HALF_T + 1),
                             max(0, xc - RETRACE_HALF_X):min(w_img, xc + RETRACE_HALF_X + 1)] = 0.0
            t_new, xs_new = trace_single_cell_dp_ridge(
                k_masked, t0_i, l_smooth, l_thresh, x_start_min=1225, x_exit_cutoff=14,
                exit_max_step=5, exit_target_ramp=(1.2, 1.5, 4.0), exit_rule="noise_biased"
            )
            if t_new is not None and xs_new:
                if DEBUG_ENTRANCES and len(xs_new) != len(xs_i):
                    print(f"      [DEBUG] re-trace t_enter={t0_i}: length {len(xs_i)} -> {len(xs_new)}", flush=True)
                raw_extracted[idx] = (t0_i, xs_new)

    raw_extracted.sort(key=lambda item: item[0])

    final_tracks = {}
    tid = 1

    for actual_t_enter, dp_xs in raw_extracted:
        t_end = actual_t_enter + len(dp_xs) - 1
        print(
            f"      [Stage 3B/3C] Refining CELL #{tid} (Frames {start_frame + actual_t_enter} -> {start_frame + t_end}) "
            f"via 2D Kymograph Peeling + RTS Kalman Smoother...",
            flush=True
        )

        observations = []
        prev_g_roi = None

        for step_idx, gx in enumerate(dp_xs):
            cache_idx = actual_t_enter + step_idx
            if cache_idx >= len(roi_cache):
                break
            f_idx = start_frame + cache_idx
            gray_roi = roi_cache[cache_idx].astype(np.float32)

            obs = extract_subpixel_observation(
                gray_roi, prev_g_roi, bg_roi, mask_roi, gx, y_min, top_wall, bot_wall
            )
            obs["frame"] = f_idx
            obs["track_id"] = tid
            observations.append(obs)
            prev_g_roi = gray_roi

        smoothed = kalman_rts_smooth_passage(observations, top_wall, bot_wall)
        if smoothed:
            final_tracks[tid] = smoothed
            tid += 1

    return final_tracks


# ====================================================================
# STAGE 3D: EXIT-FLARE EXTENSION BY PER-FRAME COMPONENT TRACKING
# ====================================================================

def extend_exit_flare(video_path, tracks, bg, lumen_mask, top_wall, bot_wall,
                      y_min, y_max, verbose=True):
    """Carry a track through the exit flare using per-frame connected components.

    Why components rather than the kymograph ridge: measured on cell #2, a per-frame
    component inside the lumen follows the cell cleanly out of the flare (x 75 -> 9 between
    frames 2730 and 2790 at ~1.2 px/frame) while the DP ridge stops ~80 frames early. The
    ridge is a 1-D abstraction with a speed prior; the component is a direct measurement.

    Static structure is handled by motion gating (difference from the previous frame), and
    identity by a monotonic-x + area-plausibility + nearest-neighbour association. Only the
    frames AFTER the existing track's last point are written, so nothing upstream changes.
    """
    if not tracks:
        return tracks

    bg_roi = bg[y_min:y_max, :]
    mask_roi = lumen_mask[y_min:y_max, :]
    h_roi, w_img = bg_roi.shape

    cap = cv2.VideoCapture(video_path)
    cache = {}

    def get_gray(f):
        if f not in cache:
            cap.set(cv2.CAP_PROP_POS_FRAMES, f)
            ok, img = cap.read()
            cache[f] = (cv2.cvtColor(img, cv2.COLOR_BGR2GRAY).astype(np.float32)
                        if ok else None)
        return cache[f]

    for tid, path in tracks.items():
        if len(path) < 20:
            continue
        last = path[-1]
        if last["cx"] > 260:                 # track did not reach the exit flare
            continue

        x_cur = float(last["cx"])
        v_last = float(np.mean([path[-1]["cx"] - path[k]["cx"]
                                for k in range(max(0, len(path) - 12), len(path) - 1)] or [1.0]))
        if v_last > 0:
            v_last = 1.0
        added = 0
        prev_rows = None

        for step in range(1, 201):
            f = last["frame"] + step
            g = get_gray(f)
            gprev = get_gray(f - 1)
            if g is None or gprev is None:
                break

            cur = g[y_min:y_max, :]
            prv = gprev[y_min:y_max, :]
            mot = np.abs(cur - prv)
            mot[mask_roi == 0] = 0.0
            absd = np.abs(cur - bg_roi)
            absd[mask_roi == 0] = 0.0

            # motion-gated: static walls and the fixed flare feature do not move
            energy = mot * np.sqrt(np.maximum(1.0, absd))
            energy = cv2.GaussianBlur(energy, (5, 5), 0)

            m = (energy >= max(3.0, float(np.percentile(energy[mask_roi > 0], 99.0)))).astype(np.uint8)
            nn, lab, st, ce = cv2.connectedComponentsWithStats(m, 8)

            best = None
            for j in range(1, nn):
                a = int(st[j, cv2.CC_STAT_AREA])
                if a < 20:
                    continue
                cx = float(ce[j][0])
                # must be at or left of the current position, within a plausible step
                if cx > x_cur + 1.0 or (x_cur - cx) > 6.0:
                    continue
                d = abs(cx - x_cur)
                score = d - 0.02 * a
                if best is None or score < best[0]:
                    best = (score, cx, a, float(ce[j][1]))

            if best is None:
                break

            _, cx_new, area, cy_new = best
            cy_abs = float(y_min + cy_new)
            bw0 = int(np.clip(round(4.0 * np.sqrt(max(area, 4) / np.pi)), 6, 40))
            bh0 = int(np.clip(last["bbox"][3], 4, 30))

            path.append({
                "frame": f,
                "track_id": tid,
                "cx": cx_new,
                "cy": cy_abs,
                "bbox": (int(round(cx_new - bw0 / 2)), int(round(cy_abs - bh0 / 2)), bw0, bh0),
                "area": float(area),
                "circularity": last["circularity"],
                "elongation": last["elongation"],
                "velocity": float(np.hypot(cx_new - x_cur, cy_abs - last["cy"])),
                "ADI_running": last["ADI_running"],
            })
            x_cur = cx_new
            added += 1

            if cx_new <= 2.0:
                break

        if verbose and added:
            print(f"      [Stage 3D] CELL #{tid}: exit flare extended by {added} frames "
                  f"({last['frame']} -> {path[-1]['frame']}), "
                  f"x {last['cx']:.0f} -> {path[-1]['cx']:.0f}", flush=True)

    cap.release()
    return tracks


# ====================================================================
# BIOMECHANICS SUMMARY
# ====================================================================

def compute_biomarker(track_id, path):
    frames = [p["frame"] for p in path]
    areas = [p["area"] for p in path]
    circs = [p["circularity"] for p in path]
    elongs = [p["elongation"] for p in path]
    vels = [p["velocity"] for p in path[1:] if p["velocity"] > 0]

    area_mean = float(np.mean(areas))
    area_min = float(np.min(areas))
    adi = area_min / area_mean if area_mean > 0 else 0.0

    shape_recovery = None
    if len(path) >= 5:
        mid_circ = float(np.mean(circs[len(circs) // 3: 2 * len(circs) // 3]))
        post_circ = float(np.mean(circs[2 * len(circs) // 3:]))
        if mid_circ > 0:
            shape_recovery = post_circ / mid_circ

    x_start, x_end = float(path[0]["cx"]), float(path[-1]["cx"])
    span_px = float(abs(x_start - x_end))

    return {
        "track_id": int(track_id),
        "entrance_frame": int(frames[0]),
        "exit_frame": int(frames[-1]),
        "track_length_frames": len(path),
        "duration_ms": float((frames[-1] - frames[0]) / 10.0),
        "x_start_px": round(x_start, 1),
        "x_end_px": round(x_end, 1),
        "net_travel_px": round(span_px, 1),
        "velocity_mean": float(np.mean(vels)) if vels else 0.0,
        "velocity_max": float(np.max(vels)) if vels else 0.0,
        "area_mean": area_mean,
        "area_min": area_min,
        "circularity_mean": float(np.mean(circs)),
        "circularity_min": float(np.min(circs)),
        "elongation_max": float(np.max(elongs)),
        "ADI": float(adi),
        "shape_recovery": float(shape_recovery) if shape_recovery is not None else None,
    }


# ====================================================================
# RENDER ANNOTATED VIDEO
# ====================================================================

def _cell_color_bgr(track_id):
    palette = [
        (50, 255, 100),   # Cell #1: Neon Green
        (0, 215, 255),    # Cell #2: Gold / Amber
        (255, 200, 50),   # Cell #3: Cyan-Blue
        (200, 110, 255),  # Cell #4: Magenta
    ]
    return palette[(track_id - 1) % len(palette)]


class _H264Writer:
    """MP4 writer that produces H.264, so the video plays in browsers, editor/app viewers
    and stock Windows players. OpenCV's 'mp4v' (MPEG-4 Part 2) does not play in those.

    Preferred: pipe raw BGR frames to ffmpeg (libx264, CRF 18, yuv420p, +faststart) when
    ffmpeg is on PATH -- one encode, no intermediate file. Fallbacks when ffmpeg is missing:
    OpenCV 'avc1' (H.264 where the OpenCV build has an encoder; measured here: ~7x larger
    than libx264 at similar quality), then 'mp4v' as a last resort. `kind` says which was used.
    Same write()/release() interface as cv2.VideoWriter.
    """

    def __init__(self, path, fps, size):
        self.w, self.h = int(size[0]), int(size[1])
        self.proc = None
        self.cv_writer = None
        self.kind = None
        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg:
            cmd = [
                ffmpeg, "-y", "-loglevel", "error",
                "-f", "rawvideo", "-pix_fmt", "bgr24", "-s", f"{self.w}x{self.h}",
                "-r", f"{fps}", "-i", "-",
                "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",   # yuv420p needs even dimensions
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-preset", "medium",
                "-movflags", "+faststart", path,
            ]
            try:
                self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
                self.kind = "H.264 (ffmpeg libx264)"
            except OSError:
                self.proc = None
        if self.proc is None:
            for tag in ("avc1", "mp4v"):
                wr = cv2.VideoWriter(path, cv2.VideoWriter_fourcc(*tag), fps, (self.w, self.h))
                if wr.isOpened():
                    self.cv_writer = wr
                    self.kind = ("H.264 (OpenCV avc1)" if tag == "avc1"
                                 else "mp4v / MPEG-4 Part 2 (no H.264 encoder found; may not play in browsers)")
                    break
                wr.release()
            if self.cv_writer is None:
                raise RuntimeError(f"Could not open any video writer for {path}")

    def write(self, frame):
        if frame.shape[0] != self.h or frame.shape[1] != self.w:
            frame = cv2.resize(frame, (self.w, self.h))
        if self.proc is not None:
            self.proc.stdin.write(np.ascontiguousarray(frame, dtype=np.uint8).tobytes())
        else:
            self.cv_writer.write(frame)

    def release(self):
        if self.proc is not None:
            self.proc.stdin.close()
            rc = self.proc.wait()
            if rc != 0:
                raise RuntimeError(f"ffmpeg exited with code {rc} while writing the video")
        elif self.cv_writer is not None:
            self.cv_writer.release()


def render_tracked_video(
    video_path,
    top_wall,
    bot_wall,
    valid_tracks,
    biomarkers_by_id,
    out_dir,
    pad_frames=20,
    playback_fps=20.0
):
    out_mp4 = os.path.join(out_dir, "tracked_cell.mp4")
    legacy_mp4 = os.path.join(out_dir, "active_passages_only.mp4")

    by_frame = defaultdict(list)
    for tid, path in valid_tracks.items():
        for pt in path:
            by_frame[pt["frame"]].append(pt)

    cap = cv2.VideoCapture(video_path)
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    w = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)) or 1280
    h = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT)) or 480

    writer = _H264Writer(out_mp4, playback_fps, (w, h))
    print(f"      [Video] encoder: {writer.kind}", flush=True)

    frames_set = set()
    if valid_tracks:
        for tid in sorted(valid_tracks.keys()):
            path = valid_tracks[tid]
            f0 = max(0, path[0]["frame"] - pad_frames)
            f1 = min(total_frames - 1, path[-1]["frame"] + pad_frames)
            for f in range(f0, f1 + 1):
                frames_set.add(f)
        frames_to_render = sorted(frames_set)
    else:
        frames_to_render = list(range(1500, min(1950, total_frames)))

    trail_by_id = {tid: [] for tid in valid_tracks}
    current_pos = -999

    for f_idx in frames_to_render:
        if f_idx != current_pos:
            cap.set(cv2.CAP_PROP_POS_FRAMES, f_idx)
        ok, frame = cap.read()
        if not ok:
            break
        current_pos = f_idx + 1

        vis = frame.copy() if frame.ndim == 3 else cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR)

        active_cells = by_frame.get(f_idx, [])
        for idx_in_frame, d in enumerate(active_cells):
            tid = d["track_id"]
            color = _cell_color_bgr(tid)
            bx, by, bw, bh = d["bbox"]
            cx_i, cy_i = int(round(d["cx"])), int(round(d["cy"]))
            trail_by_id[tid].append((cx_i, cy_i))

            trail = trail_by_id[tid]
            if len(trail) >= 2:
                cv2.polylines(vis, [np.array(trail, np.int32)], False, color, 2, cv2.LINE_AA)
                cv2.circle(vis, trail[0], 4, (0, 215, 255), -1, cv2.LINE_AA)

            pad_x = max(10, bw // 2 + 3)
            pad_y = max(9, bh // 2 + 3)
            x0, y0 = max(0, cx_i - pad_x), max(0, cy_i - pad_y)
            x1, y1 = min(w - 1, cx_i + pad_x), min(h - 1, cy_i + pad_y)
            cv2.rectangle(vis, (x0, y0), (x1, y1), color, 2)

            bio = biomarkers_by_id.get(tid, {})
            adi_val = bio.get("ADI", d["ADI_running"])
            zone = "Right Entrance" if cx_i > 1180 else ("Left Exit" if cx_i < 110 else "Constriction Bore")
            line1 = f"CELL #{tid} [{zone}] x={cx_i}px"
            line2 = f"Vel: {d['velocity']:.1f} px/f  |  Area: {d['area']:.0f} px2"
            line3 = f"Circ: {d['circularity']:.2f}  |  ADI: {adi_val:.2f}  |  Elong: {d['elongation']:.1f}x"

            card_w, card_h = 295, 60
            card_x = int(np.clip(cx_i - card_w // 2, 10, w - card_w - 10))
            place_above = (idx_in_frame % 2 == 0) and (y0 > card_h + 38)
            card_y = y0 - card_h - 14 if place_above else min(h - card_h - 6, y1 + 14)

            roi_hud = vis[card_y:card_y + card_h, card_x:card_x + card_w]
            if roi_hud.shape[0] == card_h and roi_hud.shape[1] == card_w:
                dark_box = np.full_like(roi_hud, 16)
                vis[card_y:card_y + card_h, card_x:card_x + card_w] = cv2.addWeighted(roi_hud, 0.2, dark_box, 0.8, 0)
            cv2.rectangle(vis, (card_x, card_y), (card_x + card_w, card_y + card_h), color, 1)

            leader_y = card_y + card_h if card_y < y0 else card_y
            box_anchor_y = y0 if card_y < y0 else y1
            cv2.line(vis, (cx_i, leader_y), (cx_i, box_anchor_y), color, 1, cv2.LINE_AA)

            cv2.putText(vis, line1, (card_x + 8, card_y + 18), cv2.FONT_HERSHEY_SIMPLEX, 0.44, color, 1, cv2.LINE_AA)
            cv2.putText(vis, line2, (card_x + 8, card_y + 36), cv2.FONT_HERSHEY_SIMPLEX, 0.40, (240, 240, 240), 1, cv2.LINE_AA)
            cv2.putText(vis, line3, (card_x + 8, card_y + 53), cv2.FONT_HERSHEY_SIMPLEX, 0.38, (180, 225, 255), 1, cv2.LINE_AA)

        cv2.rectangle(vis, (0, 0), (w, 28), (14, 17, 22), -1)
        if active_cells:
            ids_str = ", ".join(f"CELL #{c['track_id']}" for c in active_cells)
            status_txt = f"RBC Passage Tracker | Frame {f_idx} ({f_idx / 10.0:.1f} ms) | Tracking {ids_str}"
            txt_col = (120, 255, 150)
        else:
            status_txt = f"RBC Passage Tracker | Frame {f_idx} ({f_idx / 10.0:.1f} ms) | No cell in passage"
            txt_col = (170, 175, 185)
        cv2.putText(vis, status_txt, (14, 19), cv2.FONT_HERSHEY_SIMPLEX, 0.48, txt_col, 1, cv2.LINE_AA)

        writer.write(vis)

    writer.release()
    cap.release()

    try:
        shutil.copyfile(out_mp4, legacy_mp4)
    except Exception:
        pass

    return out_mp4, len(frames_to_render)


# ====================================================================
# MAIN PIPELINE
# ====================================================================

def main():
    parser = argparse.ArgumentParser(description="Stage 3: Multi-Cell Kymograph Peeling + RTS Kalman Smoother")
    parser.add_argument("--max-frames", type=int, default=1500, help="Frames to scan (default: 1500; pass 0 for all 37,274)")
    parser.add_argument("--start-frame", type=int, default=1450, help="Frame index to start from (default: 1450)")
    parser.add_argument("--video", type=str, default=None, help="Optional override path to .avi file")
    args = parser.parse_args()

    if not os.path.isfile(CALIB_FILE) or not os.path.isfile(BG_FILE):
        print(f"Calibration files not found in {OUT}. Running preprocessing.py (Stage 1) first...")
        import preprocessing
        preprocessing.main()

    calib = load_calibration()
    background = load_background()

    video_path = args.video if args.video else _resolve_video_path(calib.get("video"))
    fps = calib.get("fps", 10000)

    print("=" * 68)
    print("STAGE 3: MULTI-CELL KYMOGRAPH PEELING + RTS KALMAN SMOOTHER")
    print("=" * 68)
    print(f"Video path : {video_path}")
    print(f"Output dir : {OUT}")

    lumen_mask, top_wall, bot_wall, slope, intercept = build_exact_lumen_geometry(background)
    ys_active = np.where(lumen_mask > 0)[0]
    roi_y_bounds = (max(0, int(ys_active.min()) - 6), min(background.shape[0], int(ys_active.max()) + 6))

    out_dir = os.path.join(OUT, "stage3")
    os.makedirs(out_dir, exist_ok=True)

    cap = cv2.VideoCapture(video_path)
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    cap.release()

    frames_to_run = total_frames - args.start_frame
    if args.max_frames and args.max_frames > 0:
        frames_to_run = min(frames_to_run, args.max_frames)

    print(f"\n[1/2] Running Multi-Cell Kymograph Peeling + RTS Kalman Smoother on {frames_to_run} frames (start={args.start_frame})...")

    valid_tracks = detect_and_track_passages(
        video_path=video_path,
        start_frame=args.start_frame,
        frames_to_run=frames_to_run,
        background=background,
        lumen_mask=lumen_mask,
        top_wall=top_wall,
        bot_wall=bot_wall,
        roi_y_bounds=roi_y_bounds,
    )

    valid_tracks = extend_exit_flare(
        video_path, valid_tracks, background, lumen_mask, top_wall, bot_wall,
        roi_y_bounds[0], roi_y_bounds[1],
    )

    all_trajectories = []
    all_biomarkers = []
    biomarkers_by_id = {}

    for tid, path in valid_tracks.items():
        bio = compute_biomarker(tid, path)
        if bio:
            all_biomarkers.append(bio)
            biomarkers_by_id[tid] = bio
        for p in path:
            bx, by, bw, bh = p["bbox"]
            all_trajectories.append({
                "frame": p["frame"],
                "track_id": tid,
                "cx": round(p["cx"], 2),
                "cy": round(p["cy"], 2),
                "bbox_x": bx,
                "bbox_y": by,
                "bbox_w": bw,
                "bbox_h": bh,
                "velocity_px_per_frame": round(p["velocity"], 3),
                "area": round(p["area"], 2),
                "circularity": round(p["circularity"], 4),
                "elongation": round(p["elongation"], 3),
                "ADI_running": round(p["ADI_running"], 4),
            })

    traj_csv = os.path.join(out_dir, "trajectories.csv")
    with open(traj_csv, "w", newline="", encoding="utf-8") as f:
        if all_trajectories:
            writer = csv.DictWriter(f, fieldnames=list(all_trajectories[0].keys()))
            writer.writeheader()
            writer.writerows(all_trajectories)

    bio_csv = os.path.join(out_dir, "biomechanics.csv")
    with open(bio_csv, "w", newline="", encoding="utf-8") as f:
        if all_biomarkers:
            writer = csv.DictWriter(f, fieldnames=list(all_biomarkers[0].keys()))
            writer.writeheader()
            writer.writerows(all_biomarkers)

    summary = {
        "total_tracks": len(all_biomarkers),
        "total_frames_analyzed": frames_to_run,
        "start_frame": args.start_frame,
        "fps": fps,
        "video": video_path,
        "metrics": {
            "velocity_mean_px_per_frame": float(np.mean([b["velocity_mean"] for b in all_biomarkers])) if all_biomarkers else None,
            "velocity_max_px_per_frame": float(np.max([b["velocity_max"] for b in all_biomarkers])) if all_biomarkers else None,
            "ADI_mean": float(np.mean([b["ADI"] for b in all_biomarkers])) if all_biomarkers else None,
            "circularity_mean": float(np.mean([b["circularity_mean"] for b in all_biomarkers])) if all_biomarkers else None,
        }
    }
    with open(os.path.join(out_dir, "track_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n[2/2] Rendering clean passage video: output/stage3/tracked_cell.mp4 ...")
    out_mp4, rendered_frames = render_tracked_video(
        video_path=video_path,
        top_wall=top_wall,
        bot_wall=bot_wall,
        valid_tracks=valid_tracks,
        biomarkers_by_id=biomarkers_by_id,
        out_dir=out_dir,
        pad_frames=20,
        playback_fps=20.0,
    )

    print("\n" + "=" * 68)
    print("STAGE 3 COMPLETE")
    print("=" * 68)
    print(f"  ✓ Tracked Video File   : {out_mp4} ({rendered_frames} frames)")
    print(f"  ✓ Total Tracked Cells  : {len(all_biomarkers)}")

    if all_biomarkers:
        print("\nTracked Cell Passages (Bottom-Right Entrance -> Top-Left Exit):")
        print(f"  {'CellID':>6} | {'Frames':>6} | {'FrameSpan':>13} | {'X Passage Span':>16} | {'VelMean':>7} | {'AreaMean':>8} | {'ADI':>5} | {'Circ':>5}")
        print("  " + "-" * 82)
        for b in all_biomarkers:
            fspan = f"{b['entrance_frame']}->{b['exit_frame']}"
            xspan = f"{b['x_start_px']:.0f}->{b['x_end_px']:.0f}px"
            print(
                f"  #{b['track_id']:5d} | {b['track_length_frames']:6d} | {fspan:>13} | {xspan:>16} | "
                f"{b['velocity_mean']:7.2f} | {b['area_mean']:8.1f} | {b['ADI']:5.2f} | {b['circularity_mean']:5.2f}"
            )
    print(f"\nOpen the video:\n  start {out_mp4}")


if __name__ == "__main__":
    main()