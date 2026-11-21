# What did not go right

Newest at the bottom. Log misses with the same care as successes. A working overlay is not a reason to delete earlier failures.

Template:

```markdown
## YYYY-MM-DD — short name

- **Tried:**
- **What went wrong:** (what we saw, not a vibe)
- **Why (we think):**
- **What we changed:**
- **Still wrong:**
- **Evidence:**
```

---

## 2026-09-16 — Annotated screenshot as if it were a raw frame

- **Tried:** Use `image (3).png` to invent detection.
- **What went wrong:** Red circle sits on the right flare. It is not in the experiment image. Gray conversion would turn it into a dark ring that can look like a wall.
- **Why (we think):** The circle is a human annotation; screenshot coordinates also may not match the AVI.
- **What we changed:** Seeked an unannotated frame from the AVI instead (24642).
- **Still wrong:** `image (3).png` is still in the folder; do not feed it to the masker by accident.
- **Evidence:** `image (3).png` vs `lab/evidence/20260916_tube_mask_probe/probe_still.png`

## 2026-09-16 — First lumen fill only kept the left flare

- **Tried:** Adaptive-threshold walls, rotate tube horizontal, fill the **largest** gap between wall runs in each column; keep only the largest connected lumen component.
- **What went wrong:** `lumen_mask.png` was a white triangle in the left bulb only. Overlay was a solid red sausage; green almost absent. Lumen fraction ~0.37%.
- **Why (we think):**
  1. A 5×5 morphological close **fused** the two close walls in the narrow bore into one white band, so there was no gap to fill.
  2. “Largest gap” is the wrong rule: bath-side gaps can beat the bore.
  3. Keeping only the biggest blob dropped the right flare if it was disconnected.
  4. Overlay painted solid walls on top of any lumen, so even a thin bore would have been hidden.
- **What we changed:** Dropped the heavy close; thinner walls (two parallel lines). Fill between wall hits near the **centerline**, peel 2 px, punch walls out. Keep all sizable lumen pieces. Overlay uses wall **edges** plus green fill.
- **Still wrong:** See dirt leak and short right tip below.
- **Evidence:** first overlay/lumen pair was overwritten by later runs on the same still. The failure mode is described here; the saved files are the **after** state.

## 2026-09-16 — Thin bore deleted by cleanup / wrong gap pick

- **Tried:** After two walls were visible, still fill “largest valid gap,” then 5×5 morphological open on the lumen.
- **What went wrong:** Flares filled; most of the narrow channel stayed empty except a broken line on the right half. Mid-tube dirt became extra blobs.
- **Why (we think):** Open 5×5 is larger than a ~few-pixel bore, so it erases the channel. Largest-gap still prefers the wrong interval when debris adds extra runs. Max gap 0.45×height is wide enough to include bath.
- **What we changed:** Envelope between min/max wall rows within 0.18×height of the median wall row. Horizontal close (3×31) in rotated space. No 5×5 open.
- **Still wrong:** Dirt on the glass still inflates the envelope locally. Far-right flare tip underfilled (wall geometry is a Y / open tip).
- **Evidence:** `lab/evidence/20260916_tube_mask_probe/overlay.png` (current), `lumen_mask.png`

## 2026-09-16 — Dirt on the tube counted as lumen/walls

- **Tried:** Same wall detector; no extra “ignore compact blobs on the wall” rule.
- **What went wrong:** Mid-channel specks (visible on the still) produce extra red rings and small green islands off the true bore.
- **Why (we think):** Adaptive threshold does not know dirt from glass. Compact dark objects attached to the wall survive the “long component” filter if they connect to the tube, or they locally widen the min–max envelope.
- **What we changed:** Not fixed yet.
- **Still wrong:** Those leaks remain on the current overlay.
- **Evidence:** mid-tube blobs on `overlay.png` / `lumen_mask.png`

## 2026-11-07 — Per-label ROI contour optimization was slower

- **Tried:** Compute each watershed label's bounding box once, then call `findContours` on a cropped label mask.
- **What went wrong:** The cropped approach took about 0.58 s/frame vs about 0.45 s/frame for the original detector on three seeked frames. An initial centroid-offset mistake also changed detections; it was caught by comparison and not retained.
- **Why (we think):** The per-pixel coordinate grouping overhead outweighed avoiding repeated full-frame scans for this data.
- **What we changed:** Replaced it with a single contour extraction over the watershed foreground mask, mapping contours back to labels to preserve the original detection order.
- **Still wrong:** The new approach is only checked on three frames; full-pipeline runtime and broader output equivalence remain unverified.
- **Evidence:** Exact detection output comparison against the pre-change `tracking.py` on frames 0, 18637, and 37273; those checks are summarized in `lab/log.md`.

## 2026-11-07 — Stage 3 saves no tracks or trajectory rows

- **Tried:** Inspected Stage 3 output and traced `CellTracker` initialization and result collection.
- **What went wrong:** `trajectories.csv` and `biomechanics.csv` are zero-byte files; summary reports zero tracks and zero analyzed frames.
- **Why:** `add_new_tracks` created tracks with `state=None`, while `associate` excludes those tracks from matching. `update` only initialized Kalman state for matched tracks, so new tracks could never reach three hits. The Kalman filter also used invalid property names and vector shapes for the installed OpenCV API. Separately, `main` initialized `frame_count = 0` but never incremented it.
- **What we changed:** Initialize each track's Kalman filter at creation using OpenCV's `statePost`/`errorCovPost` and column-vector measurement/state shapes; update matched filters through that state; increment the successful decoded-frame count. Also pass the actual track ID into `compute_biomarker` instead of incorrectly using the first frame number as its ID.
- **Still wrong:** Existing empty Stage 3 output files are not regenerated automatically; rerun Stage 3 to replace them. The full video pipeline has not been rerun.
- **Evidence:** Before the fix, local `output/stage3/trajectories.csv` and `biomechanics.csv` were 0 bytes and `track_summary.json` reported zero. After the fix, a five-detection synthetic sequence maintained one track for five hits and produced its biomarker with ID 0. On five consecutive real frames, detector counts were 35–39 and the tracker emitted 37 tracks by the third frame.

## 2026-11-15 — Frames 15000–16500: 5 of 6 cells tracked

- **Tried:** Current `tracking.py` (kymograph DP + peel + gap search + RTS smoother) on `--start-frame 15000 --max-frames 1500`.
- **What went wrong:** User counted 6 cell passages in the window; pipeline output 5 tracks. Which cell was missed is not yet identified (the run's `output/stage3` files were later overwritten by a 18000 run).
- **Why (we think) — code-review hypotheses, NOT yet confirmed by a debug run:**
  1. Window edges: DP must start at x≥1225 and needs ≥260 (really ~330–435) frames to reach x<16, so a cell already in the bore at frame 15000, or entering after ~16070, can never become a track.
  2. Stalls / queues: kymograph signal is motion-gated (`sqrt(abs * (mot-4))`), so a stalled cell has ~zero signal; DP forbids step=0 and <3 px/f in the bore, so it cannot follow a stall.
  3. Back-to-back cells rejected as duplicates: `is_duplicate_trajectory` (25th pct separation <45 px, or any order flip, or Δt_enter <22) and the gap-pass `MERGE_MARGIN_PX=20` + `keep_len<180` reject a real follower that closes up on a stalled leader. Earlier snapshots today (`C:\Users\amara\rbc_tracking\verify2_683.png`, `finalcheck_850.png`) show this pattern: cell #5 at vel 0.0 in the bore with #6 closing in, then both boxes on the same cell at the exit.
  4. Gap pass only searches gaps >1.25× median spacing and >120 frames, and skips the first 90 frames after the leader — a queued cell right behind another is never searched.
- **Debug run (same day, `RBC_DEBUG_ENTRANCES=1`, existing code unchanged) — actual cause:**
  - Pass 1 (kymograph) found 5 entrances: 15006, 15034, 15315, 15369, 15518. The cell entering at **15459** never crossed `r_thresh=23.8` (it lingers ~20 frames at x≈1215 near the entrance).
  - The gap pass recovered it (gap 369→518). It moved normally through the bore (~4.5–6.7 px/f, x 1248→734 by frame 15578), then the windowed energy-centroid follower **froze at x=731.89 for 106 frames (15578–15683, velocity 0.000)** while the real cell went on.
  - When cell #6 reached x≈750 at frame 15684, `MERGE_MARGIN_PX=20` cut #5 there. The run says "6 tracks", but #5 is a stub ending mid-bore (15459→15683, x 1253→732). Nothing covers that cell's second half (x 732→0, ~15580→~15780), so only 5 full passages exist.
  - Likely a fixed feature at x≈732: the other gap candidate (entering 15124, rejected) was also sitting at x=731.36 at its 200th frame, and the raw frames show a speck on the channel near x≈730. The follower has no motion model, so a stationary energy source in its window beats a faint moving cell.
  - Hypotheses 1 (window edges) and 4 were **not** the cause here. The "vel 0.0 at x=732" in `verify2_683.png` was this same trap, not a real stall.
- **Further finding while fixing:** the follower never sat on the real cell at all. Per-frame energy of this cell is ~3–6 vs bore noise p90 3.7, so a per-frame follower cannot hold it. The old window `[x-25, x+15]` is centred at x-5, so in pure noise it drifts left ~5 px/f (cell-like speed). Track #5 ran ~140 px ahead of the real ridge, then stuck on a debris particle on the channel at x≈731 (seen in the raw + bg-diff crop; median energy there ~22 vs ~13 elsewhere). The planned "better window" fix would not have been enough.
- **What we changed (Sai approved):** Rewrote `find_entrance_via_components` as a kymograph DP: column-max motion×contrast energy minus each column's temporal median, accepted tracks masked ±20 px, stall-tolerant steps (bore 0–9, flares 0–4, soft targets 5 / 1.5), free start in the gap, end at x≤14. Validation fits a quadratic to the bore part of the path and needs ≥2.8× the median energy found on random bore-speed lines.
- **Result:** 15000: #5 is now a full passage, 15459→15768 (x 1269→18), on the faint ridge and through the debris without stopping (ratio 4.44). The noise path in gap [34,315] is rejected (2.16). Other 5 tracks unchanged. 1450: the old code produced a ghost 4th track (1633→2067, stuck at x=730.43 for 300+ frames). New code rejects it (1.74); the 3 real tracks are identical. 18000: identical (gap pass never runs with 2 cells).
- **Still wrong:** (a) The gap pass needs ≥3 strong cells and only searches between them. It cannot run in a 2-cell window or before the first / after the last strong entrance. (b) Window edges: the 18000 window has a bright cell already in the bore at frame 18000 (exits ~18140) that cannot be tracked by design. (c) Exit-flare speed prior (1.5 px/f) may make the last ~50 frames of a track lag the cell. (d) #1/#2 in 15000 end at the same spot (x≈31, frame 15360) — not checked. (e) The 2.8× threshold is calibrated on one real faint cell and two noise paths only.
- **Evidence:** `lab/evidence/20261115_debug_15000_16500.log`; per-frame tracks in `C:\Users\amara\rbc_tracking\output\stage3\trajectories.csv` (this run; the 18000 run is backed up in `stage3_backup_18000_run`). Earlier renders: `C:\Users\amara\rbc_tracking\verify_cell5.png`, `verify2_683.png`, `finalcheck_850.png`.

## 2026-11-15 — Exit funnel: cell #1's box stopped while the cell went ahead (15000–16500)

- **Tried:** Gap-fixed `tracking.py` on 15000–16500. Sai saw cell #1's box stop at the end of the funnel with the cell ahead of it.
- **What went wrong:** Box #1 sat at x≈95–105 (funnel throat) from ~15287 to 15321, then rode with box #2 to the exit at 15376. In the background-difference crops the real cell #1 (a "C"-shaped ring) was already at x≈72 by 15297 and x≈20 by 15344, and cell #2 (a "9" ring) was behind it. Box #1 was 30–36 px behind, then on cell #2.
- **Why:** (1) Main-pass DP traces each cell alone. Cell #1's faint path was pulled onto cell #2's much brighter ridge (the two paths were within 14 px for 32 frames). (2) The DP's flare speed cap (1–2 px/f) is too slow at the exit throat, where cells leave at ~4 px/f (x 110→72 in ~13 frames) before slowing to ~1.2. (3) Once the path was >20 px behind, the per-frame measurement (±20 px window) locked onto the throat walls, which flicker in the background difference because the tube jitters frame to frame.
- **Dead ends (kept for the record):**
  - Loosening the exit-funnel speed cap in the main DP: cell #1's path jumped onto cell #2's ridge already in the bore. On 1450 it changed one candidate's path enough to fail `raw_check` and lose a real cell (3 → 2).
  - Re-tracing every track with neighbours masked: on 1450 the right exit for cell #3 (2805) moved to 2729. The left-exit bonus, computed from the unmasked kymograph, pulled it to cell #2's exit.
  - An "accept only if mean path energy does not drop" guard: biased toward shorter paths that skip weak funnel segments. It accepted moving 1450 cell #2's exit 2714→2675, which the kymograph ridge contradicts.
  - Per-frame blob detection in the funnel (local temporal median background, and phase-correlation registration): the wall flicker survived both and dominated detections.
- **What we changed (Sai said "fix that"):** `trace_single_cell_dp_ridge` gets optional `exit_max_step` / `exit_target_ramp`. Defaults keep the original model, so main-pass accept/reject is unchanged. New conflict-triggered re-trace in `detect_and_track_passages`: only tracks sharing a ridge (within ±14 px for ≥5 frames) are re-traced. Each is re-traced on the unpeeled kymograph with other tracks blanked (±14 px, ±4 frames), with the exit bonus recomputed from the masked kymograph, using exit steps 1–5 and a target easing 4.0→1.5→1.2. Two rounds; gap-recovered tracks are never re-traced.
- **Result:** 15000: only #1/#2 flagged. Box #1 at 15297/15310/15330/15346 is 77/56/~43/31 vs real ~72/55/35/22, staying on the "C" cell and ahead of #2; it ends 15353 instead of 15376. Box #2 is 4–9 px behind the "9". #3–#6 are identical. 1450 and 18000 trajectories are identical to before (no shared ridges).
- **Still wrong:** Both boxes trail the cell centre by ~7–12 px in the exit funnel (the box still overlaps the cell). The shared-ridge trigger is tested on one conflict only. Wall flicker at the throat is unaddressed.
- **Evidence:** `lab/evidence/20261115_gap_dp_fix/cell1_exit_before_vs_now_bgdiff.png`, `cell1_exit_funnel_after_fix.png`, `run_*_exit_fix.log`.

## 2026-11-15 — Cell #3 (15000) box runs ahead and ends ~47 frames early; same bug in 1450 cell #2

- **Seen:** Sai unsure about cell #3. Box #3 sprints up to ~70 px ahead of its kymograph ridge in the last part of the bore (~7 px/f), waits at the funnel throat (x≈95), then runs 15–20 px ahead and ends at 15658. In the bg-diff the cell (a "C" ring) is at x≈72/53/40/22/10 at frames 15642/15654/15666/15690/15700, leaving ~15705. Box #4 slides onto cell #3's "C" in its last ~10 frames.
- **Why:** The main DP picks the exit frame by *average* score + a left-exit bonus. The average rewards cutting weak/slow stretches, and the bonus fires when a neighbour (cell #4) enters the funnel (~15656). With the funnel capped at 1–2 px/f, the path must reach the throat early to make that exit, hence the sprint. The shared-ridge trigger missed it because box #3 ended before box #4 arrived, so their paths never overlap.
- **Also found:** 1450 cell #2 has the same bug. The faint cell is at x≈45/28/12 at 2640/2655/2670 and gone by 2685, while its box trails 40–50 px behind and ends 2714. The header comment in `tracking.py` ("true exit 2730") is contradicted by the frames.
- **Fix (Sai approved "apply it"; was prototyped offline first):** exit frame = argmax of (total DP score − noise level × duration), noise = p75 of smoothed bore kymograph; no mean, no exit bonus. Insensitive to p50/p75/p90. Re-trace every main-pass track with it (neighbours masked, faster exit model); main-pass acceptance unchanged. Offline: 15000 #1 unchanged; #3 ends 15686, 81/54/37/25 vs 72/53/40/32; #6 15870→15880 (cell still visible at x≈13 at 15882); 1450 #2 2714→2672 (matches frames); 1450 #1 1873→1885, #3 2788→2790. 18000 not yet run.
- **Evidence:** `lab/evidence/20261115_gap_dp_fix/cell3_exit_*`, `check_1450_cell2_exit_current_vs_candidate.png`, `check_15000_cell6_exit_current_vs_candidate.png`.
- **Applied result (all three windows rerun):** 15000: #3 ends 15686, box 78/59/46/33 vs real 72/53/40/32. #1 unchanged and right; #6 ends 15880; #5 untouched. 1450: #2 ends 2672 (box 62/47/19 vs real 45/28/12, was 95/77/62); #1/#3 final ends unchanged (1899/2805). 18000: #1 18333→18325 (cell gone ~18330); #2 19009→19028 (faint cell still in funnel at 19018/19028). Accepted cells identical in all windows.
- **Still wrong:** #3 still ends ~20 frames before the cell leaves. Box #4 (15000) runs ~10 px ahead in its last ~10 frames, partly onto cell #3. Boxes generally trail the cell centre by ~5–19 px in the exit funnel. The `tracking.py` header comment's "true exits 2730 / 2805" for 1450 is contradicted for cell #2 (frames show exit ~2675) — not edited.

## 2026-11-16 — Gap pass: window-placement dependence and an alternate-frame flicker bug

- **Tried:** The gap pass now searches every interval between accepted cells, plus before the first and after the last (Sai approved). Verified on 1450/15000/18000 (1500-frame windows) and the runner's test window (core 15000–16500, decoded 14550–16950).
- **Bug found and fixed:** the gap DP's column baseline sampled every 2nd frame, so it used only odd or only even frames depending on the interval start. This footage flickers between alternate frames, so the noise reference jumped 1.91 ↔ 3.42 when the start moved by one frame. Now all frames are used, and the noise reference is stable at ~2.4–2.8.
- **Calibration after the fix:** no-cell paths score 1.13–1.86×. Recovered real cells score 3.52× (15459, 1500-frame window) and 3.16× (15075, runner window). The 2.8× threshold still separates these.
- **Still wrong:** (a) The real cell entering ~15776 scores only 2.05× in the 1500-frame window (rejected) — too close to the noise maximum (1.86×) for any safe threshold. (b) In the runner's window placement, the faint 15459 cell is not recovered (1.74×): the cell ahead of it is placed at entrance 15350 there (15369 in the 1500 window), and its entry-flare path overlaps the faint cell, so the ±20 px mask hides the faint cell's start. Net effect: in busy windows the count can differ by about one faint cell depending on window placement.
- **Ideas not yet tried:** run two offset window grids and merge detections (a cell missed in one placement is often found in the other); blank accepted tracks only in the bore, not in the entry flare where DP paths are approximate.
- **Evidence:** test logs in the session scratchpad; runner outputs in `OneDrive\Desktop\RBC_test_window_15000-16500`.

## 2026-11-18 — Bore-only masking in the gap pass (prototype, rejected)

- **Tried:** To recover the faint 15459 cell in placements where the cell ahead (entrance 15350) is blanked ±20 px over the faint cell's start, a scratchpad prototype stripped entry-funnel points (x > 1180) from the gap pass's exclude tracks, so accepted cells were blanked only in the bore and exit funnel. `tracking.py` was not edited.
- **What went wrong:** worse on every busy window. Runner-style window 14550–16950: faint cell still missed, AND the 15075 recovery lost (3.16× → 1.88×). 15000–16500 window: the faint 15459 cell, previously recovered at 3.52×, lost (1.33×). 1450 and 18000 unchanged.
- **Why:** with the entry funnel unmasked, the gap DP latches onto the bright entrance of the NEXT accepted cell. In every interval the best path started exactly at that cell's entrance (local frames 514, 752, 862, 964), then scored low and was rejected. The ±20 px entry-funnel blanking is what keeps the gap search off neighbouring entrances.
- **What we changed:** nothing (prototype discarded).
- **Still wrong:** the faint 15459 cell is found only when the window starts at 15000. The offset grid does not fix this case (both grids place the cell ahead at 15350).
- **Evidence:** prototype log in the session scratchpad (`proto_bore_mask.log`).

## 2026-11-18 — Offset grid's first whole-video merge added noise; base result has noise tracks too

- **Tried:** Whole-video run of grid B (cores offset 750) merged into grid A by path matching only.
- **What went wrong:** 120 cells = 106 (A) + 14 added by B. Entry-funnel crops at each added track's start: only 3 show a real cell (entering 6628, 7260, 23898); 10 show an empty funnel, and 1 (32183) is the cell entering as the next track. Ridge evidence along the added tracks (energy along the track vs random lines): real 2.75–3.34, noise 1.11–1.92. Ten of the 14 lasted exactly 327–331 frames, the main pass's minimum passage length, which is what a path through noise produces.
- **Why (we think):** grid B's windows have their own entrance thresholds; in quiet stretches the main pass accepts paths through noise, and path matching cannot reject a path that matches nothing.
- **What we changed:** `run_full_video.py` adds a grid-B track only if its ridge evidence is ≥ 2.5 (`RIDGE_MIN_RATIO`); otherwise verdict "rejected". Scores cached in `ridge_scores.json`, shown in `merge_report.csv`.
- **Also found — base result:** scoring all 120 tracks gave 29 below 2.0 and 13 at 2.0–2.5. Entry-funnel crops of all 42: grid-A tracks entering 27728, 35260, 35322, 35760, 7830 (#21), ~7700 (#20), ~20000 (#53) and A+B tracks entering 28232, 28292, 29101, ~26800 (#72) show no cell, so roughly 8–11 of the original 106 are likely noise (clustered at 27700–29100 and 35200–35800). But 5 clearly real bright cells also scored 1.2–1.8 (entering 32880, 21583, 28924, 8938, 26164), probably because their paths leave the ridge further down the channel. So the ridge score must not be used to delete grid-A tracks.
- **Still wrong:** the base-result false positives; low-scoring real cells with off-ridge paths. A better per-track check is needed (candidate: detect the cell as a ring in the entry funnel at the track's start, where cells are clearly visible).
- **Evidence:** `lab/evidence/20261118_offset_grid/` (entry-funnel crops for the additions and for all 42 low-scoring tracks, track scores, kymograph of 30000–33000).
