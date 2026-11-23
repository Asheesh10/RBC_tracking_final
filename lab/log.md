# Log

Newest entries at the bottom.

## 2026-09-16 — Stored handoff and mapped local files

- **Done:** Wrote `RBC_handoff_context.txt` (sections 1–12 plus folder inventory). Did not change the notebook or raw companions.
- **Evidence:** none (inventory only)
- **Open:** Workspace search had missed large/ignored files; directory listing found them.
- **Next:** Use the handoff file instead of re-pasting context.

## 2026-09-16 — Local inventory

- **Done:** Listed `C:\Users\ASUS\Documents\IIsc`.
- **Evidence:** none. Facts: AVI is `2A002 2_C001H001S0001-001.avi` (22,901,777,870 bytes). CIHX 16,603 bytes. MII 1,139 bytes. Notebook `Untitled24.ipynb`. No `RBC Videos Longevity India` folder. No separate `RBC_01_Inspect.ipynb` filename (inspection cells live in `Untitled24.ipynb`).
- **Open:** AVI not decoded. Inspection notebook still points at Colab paths and the name without `-001`.
- **Next:** Either retarget inspection or, as later decided, start with single-frame detection on a still.

## 2026-09-16 — Sample frame and detection game plan

- **Done:** User added `image (3).png`. Confirmed diagonal tube, bath full of cell-like objects, target in the **right-hand flare** (red circle = annotation). Agreed the next problem is algorithmic single-frame detection, broken into: tube ROI → optional contrast → blobs → size/shape/location filters. Wrote `lab/gameplan.md`.
- **Evidence:** `image (3).png` (scouting still, annotated; not a raw AVI frame export)
- **Open:** No detector run yet. Mask not drawn. Narrow-section appearance unknown.
- **Next:** Try the four-layer recipe on this PNG without using the red circle as a feature.

## 2026-09-16 — Lab notebook + logging skill

- **Done:** Added `lab/` (CURRENT, log, decisions, gameplan, evidence/) and project skill `.cursor/skills/rbc-lab-log/SKILL.md` so agents keep this notebook updated.
- **Evidence:** none
- **Open:** Skill must be applied in later sessions; humans can still edit `lab/CURRENT.md` by hand.
- **Next:** First detection experiment on `image (3).png`.

## 2026-09-16 — Permission-before-code + tube-ROI proposal only

- **Done:** Wrote `AGENTS.md`: do not create or modify code until Sai permits; propose and wait. Logged that in `lab/decisions.md`. Wrote the tube-mask + probe-frame proposal into `lab/gameplan.md`. No scripts or notebook edits.
- **Evidence:** none
- **Open:** Awaiting permission to implement any tube-detection code.
- **Next:** Sai reviews the proposal (still → overlay; then ~3 seeked AVI frames if needed).

## 2026-09-16 — Tube mask on one random AVI frame

- **Done:** Seeked frame **24642 / 37274** (seed 20260916) from `2A002 2_C001H001S0001-001.avi`. Built wall + lumen masks. Did not scan the full video. Reused the saved still while tuning fill logic.
- **Evidence:** `lab/evidence/20260916_tube_mask_probe/` (`probe_still.png`, `overlay.png`, `lumen_mask.png`, `walls.png`). Green ribbon follows the bore and both flares; small leaks at mid-tube dirt.
- **Open:** Is this mask good enough to restrict cell search? Right flare tip slightly short.
- **Next:** User review of overlay; then blobs inside the lumen only.

## 2026-09-16 — Failure log added

- **Done:** Created `lab/failures.md` and wrote the tube-mask dead-ends (annotated PNG, fused walls / left-flare-only fill, thin bore erased by open + wrong gap rule, dirt leaks). Skill and `AGENTS.md` now require logging misses. Earlier failed overlays were overwritten; the write-up is the record of those attempts.
- **Evidence:** `lab/failures.md`; current after-state still in `lab/evidence/20260916_tube_mask_probe/`
- **Open:** Dirt leaks and short right tip still unfixed.
- **Next:** Overlay review, then cell blobs inside the lumen if the ribbon is accepted.

## 2026-10-05 — Repo, sync script, harness-agnostic rules

- **Done:** Added `.gitignore` (AVI excluded), `README.md` system map, extended `AGENTS.md`, `scripts/sync.ps1` + `sync.sh`, Copilot/CLAUDE/Cline pointers, Cursor skills `rbc-git-sync` and `rbc-writeup`, `lab/related-work.md`, `writeup/README.md`. Tube-mask code unchanged.
- **Evidence:** none (infra). Confirm `git status` does not list the AVI.
- **Open:** GitHub CLI (`gh`) is not installed on this machine; local commit exists, remote not created.
- **Next:** Install GitHub CLI, then `gh repo create iisc-rbc-tracking --private --source=. --remote=origin --push`. Overlay review for the lumen mask; other person clones and copies the AVI locally.

## 2026-10-05 — Private GitHub remote

- **Done:** Installed GitHub CLI; created private repo and pushed `master`. URL: https://github.com/Sai2125/iisc-rbc-tracking . No AVI in `git ls-files`.
- **Evidence:** none (infra)
- **Open:** Collaborator needs repo access + local AVI copy.
- **Next:** Overlay review for the lumen mask.

## 2026-11-07 — Stage 3 detector performance probe

- **Done:** Optimized watershed contour extraction in `tracking.py` to find cell contours together, then restore label ordering. Kept watershed and measurement rules unchanged.
- **Evidence:** Compared against the pre-change detector on seeked frames 0, 18637, and 37273 of 37274; detection outputs matched exactly (34, 37, and 35 detections). Detector time was 0.105–0.108 s/frame vs 0.436–0.479 s/frame on these samples.
- **Open:** This is a three-frame detector-only comparison; full-pipeline runtime was not measured, and the 3-day estimate is not confirmed for this stage.
- **Next:** If the full run remains slow, profile the actual bottleneck before further changes.
- **Follow-up timing:** A 500-frame sequential sample from frame 18637 took 64.96 s (0.1299 s/frame) including decode, detection, and tracking. Linear extrapolation to 37274 frames is about 80.7 minutes; postprocessing and CSV output are excluded.
- **Stage 3 output diagnosis:** The completed output CSVs were empty because `CellTracker` creates `state=None` tracks, but association ignores tracks without a state, so the Kalman initialization path is unreachable. The summary frame counter is also never incremented. See `lab/failures.md`; a code fix and rerun are needed.
- **Stage 3 tracking fix:** Initialized Kalman state when a detection creates a track, corrected OpenCV matrix property names and column-vector shapes, incremented the successfully decoded frame counter, and passed the true track ID to biomarker calculation. A five-detection synthetic sequence reached five hits; on five consecutive real frames, 37 tracks were emitted by the third frame. Full-video rerun not performed; existing empty output files remain until Stage 3 is rerun.

## 2026-11-15 — Zone 1/3 redesign: new module, detection+Kalman+Hungarian prototype

- **Done:** Designed (in chat, not yet in `lab/`) a revamped algorithm: 3 zones by channel geometry (entry funnel / capillary / exit funnel), each with a different association method. Zone 1 & 3 (funnels): background+motion-fused detection, Otsu threshold (per-frame, self-calibrating), bilateral filter (not blur, to protect shape metrics), watershed split for touching cells, per-track Kalman filter + Hungarian assignment (scipy `linear_sum_assignment`), gated to a max-displacement radius - replacing the old plan of raw K-nearest-neighbour matching (which would reintroduce the exact crossing-paths failure already diagnosed against iCLOTS/Trackpy). Zone 2 (capillary): kymograph fully replaces per-frame segmentation (user's call) for free track identity + transit velocity, plus a derived `ADI_proxy = (capillary_width_px * streak_length_in_kymograph) / resting_area` using the real measured ratio (resting RBC diameter ~18px vs capillary width ~9px on this footage, i.e. cells roughly double the channel width and must deform to fit) - labelled a proxy, not true ADI, until Stage 4 validates it against Zone 1/3 segmented ADI for the same cells.
- **New file (not touching `tracking.py`/`preprocessing.py`):** `zone13_tracking.py` - Zone 1/3 detection + Kalman + Hungarian prototype. Reuses existing `calibration.json` / `median_bg.npy` read-only. Zone boundaries hardcoded from calibration: `X_IN_END=110`, `X_OUT_START=1180`. ROI restricted to flare bands only (+-45px around centerline past the zone boundaries) so detection does not pick up the open bath's unrelated free-floating cells.
- **Evidence:** Ran on frames 1500-1649 (150 frames, no full-video pass). `lab/evidence/zone13_test_overlay.png` - overlay confirms detection stays correctly scoped to the two flare regions only (does not leak into the bath). 62 total tracks created; vision-inspected zoom into the crowded flare-tip corner shows only ~2-3 real distinct cell blobs there but ~15-20 track markers clustered on them - i.e. real over-fragmentation, not real cell count.
- **Open:** Per-frame Otsu on a tiny sparse ROI histogram is suspected unstable (causes detection flicker -> track death/rebirth instead of one continuous track); `GATE_RADIUS_PX=20` may also be too tight for actual per-frame displacement in the funnels (not yet measured). Architecture (fused detect + Kalman + Hungarian + ROI scoping) is validated as sound; thresholds/continuity are not yet trustworthy.
- **Next:** Measure actual per-frame cell displacement in the funnels to set a data-driven gate radius; consider a smoothed/hysteresis threshold (e.g. temporal-median Otsu, or relax `MAX_MISSES`) instead of raw per-frame Otsu, before trusting track counts. Zone 2 (kymograph + ADI_proxy) module not yet started.

## 2026-11-15 — Zone 1/3: root-caused and fixed the over-fragmentation

- **Root cause found:** `abs(gray-bg)` fired on a transient BRIGHT artifact at the extreme frame corner (x~1270-1279, y~380-386; debris/reflection/vignetting, ramps 109->255 over ~20 frames then fades - confirmed by direct pixel inspection, not a cell). Dragged per-frame Otsu threshold from ~15 up to 88-255 while active, causing severe detection flicker (tracks dying/respawning on the same real cell).
- **Fix 1 - darkening-only contrast:** RBCs in this brightfield footage are confirmed darker than background (~0-40 vs ~100-115, Stage-1 calibration). Changed detector signal from `abs(gray-bg)` to `max(0, bg-gray)` (darkening only), eliminating the whole class of brightening-artifact false positives. Threshold max dropped from 255 to 88 across the same 150-frame test.
- **Fix 2 - minimum signal floor:** Separately found that when a flare genuinely has **zero visible cells** (vision-confirmed on a raw frame export), Otsu still finds *some* split in pure noise and reports fake blobs - one such frame (fused.max()=34.9) produced 20 "detections" on an empty flare. Measured noise-floor vs real-cell-present fused.max() across 150 frames: noise stays ~35-55, real cells reach ~70-470. Added `MIN_SIGNAL_FLOOR=60` (skip detection entirely below this) and `MIN_OTSU_FLOOR=25` (clamp the computed threshold) as data-derived constants, not guesses.
- **Result:** 62 fragmented tracks -> 2 clean, continuous tracks over the same 150-frame window (1500-1649). Vision-verified: both tracks sit on real visible cell blobs at the funnel tips (not background); the earlier 15-20-marker cluster on ~2-3 real cells is gone entirely with no leakage into the open bath.
- **Gate radius check:** Measured real per-frame displacement on the 63-hit track: 0.1-2.2 px/frame (slow, near the flare tip, as expected pre-constriction-acceleration). A second, shorter track showed occasional ~6-9px jumps across missed-detection gaps. `GATE_RADIUS_PX=20` comfortably covers both - no change needed.
- **Evidence:** `lab/evidence/zone13_test_overlay.png` (updated), `lab/evidence/zone13_frame16_raw.png` / `_raw_zoom.png` (vision-confirmed empty-flare frame used to find the false-positive bug).
- **Still open:** Only tested on 150 frames / 1 funnel-entry window; have not checked exit-flare (Zone 3) separately, nor a window with multiple simultaneous cells (watershed-split path untested on real data so far). Zone 2 (kymograph + ADI_proxy) module not yet started.


## 2026-11-15 — 15000–16500 miss (5/6): code-review diagnosis

- **Done:** Read `tracking.py` Stage 3A–3C and today's diagnostic renders for this window. Listed the four gates that can drop a real cell (window edges, motion-gated kymograph + no-stall DP, duplicate/merge rejection of back-to-back cells, narrow gap-search rules). Written up in `lab/failures.md`. Reviewed Sai's sketch (KNN in funnels, kymograph + RTS in bore): zoning agrees with the 2026-11-15 redesign; recommend Kalman + Hungarian instead of plain KNN, and noted RTS only smooths, it cannot recover a cell that was never extracted.
- **Evidence:** none new (no code run). `scripts/sync.ps1 pull` failed: local `master` ahead 1 / behind 22 of `origin/master`, with uncommitted `tracking.py`; left for Sai.
- **Open:** Which of the 6 cells was missed; which gate rejected it.
- **Next:** With permission, run existing `tracking.py` with `RBC_DEBUG_ENTRANCES=1` on 15000/1500 (back up `output/stage3` first); then decide between patching gates vs the zone pipeline.

## 2026-11-15 — 15000–16500 debug pass: missed cell found

- **Done:** With Sai's permission, ran unchanged `tracking.py` with `RBC_DEBUG_ENTRANCES=1 --start-frame 15000 --max-frames 1500`. Backed up the 18000 outputs first (`output/stage3_backup_18000_run`).
- **Result:** 6 tracks reported, but #5 (entering 15459, gap-pass recovery) freezes at x=731.89 for 106 frames (15578–15683, velocity 0) and is cut at x=732 when #6 arrives. Only 5 full passages. The second half of the real cell (~15580→~15780) has no track. Likely a static feature near x≈732 that traps the energy-centroid follower; window-edge and gap-rule hypotheses ruled out for this miss. Details in `lab/failures.md`.
- **Evidence:** `lab/evidence/20261115_debug_15000_16500.log`; `output/stage3/trajectories.csv` (local).
- **Open:** Static feature at x≈732 not confirmed by a zoomed crop. #1/#2 converge at x≈31 at frame 15360 (possible double box at the exit).
- **Next:** If approved, patch `find_entrance_via_components` (motion-predicted window, static-energy suppression, lost-track re-acquire) and rerun 15000, 1450 and 18000 windows to check 6/6 with no regressions.

## 2026-11-15 — Gap pass rewritten as kymograph DP: 15000 window 6/6

- **Done (Sai approved):** Replaced the per-frame windowed-centroid follower in `tracking.py::find_entrance_via_components` with a stall-tolerant kymograph DP (static-column baseline removed, accepted tracks masked ±20 px, free start in the gap, end at the exit). Added a fit-based validation (quadratic fit on the bore part, ≥2.8× random-line noise). No other function changed. Pre-patch copy kept in the session scratchpad.
- **Result:** 15000–16500: 6 full passages. The missed cell is 15459→15768, x 1269→18, ratio 4.44; the noise path rejected at 2.16. 1450–2950: 3 cells; the old code's ghost 4th track (stuck on debris at x=730) is gone, and the 3 real tracks are byte-identical. 18000–19500: identical (2 cells; gap pass does not run).
- **Evidence:** `lab/evidence/20261115_gap_dp_fix/` (before/after kymographs, debris crop, rejected-noise strip, 18000 kymograph, run logs). Per-run outputs in `C:\Users\amara\rbc_tracking\output\stage3_*` (local); `output/stage3` holds the patched 15000 run.
- **Open:** Gap pass cannot search 2-cell windows or outside the first/last strong entrance; window-edge cells (e.g. a cell already in the bore at 18000) are untracked by design; threshold calibrated on few examples. Not committed — git sync still blocked (local `master` diverged; `tracking.py` also carries other uncommitted changes).
- **Next:** Sai to confirm 6/6 in `tracked_cell.mp4`; resolve the git divergence; then consider exit-side detection / zone pipeline so missed cells are flagged by an unmatched exit.

## 2026-11-15 — H.264 video output; exit-funnel fix for cell #1 (15000 window)

- **H.264 (Sai approved):** Sai could not play `tracked_cell.mp4` — OpenCV wrote MPEG-4 Part 2 (`mp4v`). Added `_H264Writer` in `tracking.py`: pipes frames to ffmpeg (libx264, CRF 18, yuv420p, +faststart) when ffmpeg is on PATH. Falls back to OpenCV `avc1`, then `mp4v`, and logs which encoder was used. Verified: H.264 High profile, 905 frames, 20 fps; trajectories unchanged. OpenCV `avc1` on this machine was ~7× larger at similar PSNR, so ffmpeg is preferred.
- **Exit funnel (Sai: "fix that"):** Cell #1's box stopped at the funnel throat and then rode cell #2's ridge. Root cause and four dead ends in `lab/failures.md`. Fix: conflict-triggered neighbour-masked re-trace with a faster exit-funnel speed model (main pass unchanged).
- **Result:** 15000: box #1 now follows the real cell through the funnel (~7–12 px behind its centre, previously 30–36 px behind and then on cell #2), ending at 15353 instead of 15376. #3–#6 identical. 1450 and 18000 trajectories identical to before.
- **Evidence:** `lab/evidence/20261115_gap_dp_fix/` (`cell1_exit_*`, `run_*_exit_fix.log`). Outputs: `C:\Users\amara\rbc_tracking\output\stage3_*_exitfix3` (local); `output/stage3` = final 15000 run.
- **Open:** Small residual trailing in the exit funnel. Nothing committed — git sync still blocked.

## 2026-11-15 — Noise-biased exit rule for all cells (cell #3 fix)

- **Done (Sai approved):** `trace_single_cell_dp_ridge` gained `exit_rule="noise_biased"`: exit frame = argmax(total DP score − p75 noise × duration), no averaging, no left-exit bonus. The re-trace step now applies to every main-pass track (neighbours masked ±14 px / ±4 frames, exit funnel steps 1–5 with target 4.0→1.5→1.2, two rounds). Main pass unchanged; gap-recovered tracks untouched.
- **Result:** 15000: cell #3 now follows its cell to 15686 (was 15658, ~70 px ahead before). 1450: cell #2 ends at 2672 when its cell leaves (was trailing 40–50 px to 2714). 18000: exits 18325 / 19028, both closer to the faint cells. Same cells accepted everywhere. Details and residuals in `lab/failures.md`.
- **Evidence:** `lab/evidence/20261115_gap_dp_fix/` (`check_*`, `cell3_*`, `run_*_exit_rule.log`). Final outputs `C:\Users\amara\rbc_tracking\output\stage3_*_exitrule`; `output/stage3` = 15000.
- **Open:** Box #4 (15000) slightly ahead at the end; ~5–19 px trailing in funnels; nothing committed (git sync blocked).
- **Confirmed (2026-11-15):** Sai reviewed the 15000–16500 video: cell #1 and cell #3 look right.

## 2026-11-16 — Whole-video runner built; test window shows window-dependence

- **Done (Sai approved):** `run_full_video.py` runs `tracking.py`'s Stage 3 over the whole AVI in windows: 1500-frame cores, plus 450 frames decoded before and after each. A cell belongs to the window whose core contains its entrance. Each window writes its own files atomically (`per_window/`) and the merged CSVs are rebuilt from them, so `--resume` can't duplicate rows. Outputs: trajectories, biomechanics, windows, summary, one video and one kymograph PNG per window. Default output folder is `OneDrive\Desktop\RBC_full_video_results`. Added `.vscode/launch.json` (test window + full video) and `.vscode/settings.json` (Python 3.14 with cv2); `.vscode/` is gitignored. `tracking.py` unchanged.
- **Test (core 15000–16500, decoded 14550–16950):** 5 cells kept vs yesterday's 6. Explained by debug run: (a) the "15006" cell actually enters at 14985 → owned by the previous window (correct); (b) the 15034 cell failed `raw_check` by a hair (10.61 < 10.88) and was gap-recovered at 15075; (c) the faint 15459 cell was missed: the gap pass only searches gaps >1.25× median spacing (here 266 frames; this gap is 168); (d) a NEW cell entering ~15776 was found because this window's `r_thresh` is lower (19.8 vs 23.8). It is real (round cell visible in the entrance funnel at 15790), so yesterday's 15000–16500 result missed a 7th cell (the rendered video only reached 15890).
- **Evidence:** `lab/evidence/20261115_gap_dp_fix/real_cell_15776_missed_by_1500_window.png`; test outputs in `OneDrive\Desktop\RBC_test_window_15000-16500`.
- **Open:** `tracking.py` results depend on window placement: per-window `r_thresh`, borderline `raw_check`, and the gap-pass "long gap" rule. Proposed before the full run: gap pass searches every interval between accepted cells (and before the first / after the last), relying on the fit-ratio validation. Not yet approved.

## 2026-11-16 — Gap pass searches all intervals; alternate-frame flicker bug fixed

- **Done (Sai approved):** `tracking.py` gap pass searches every interval between accepted entrances (plus before first / after last), splitting an interval when a cell is recovered. Energy is cached per window; at most 60 searches. Fixed the gap DP's every-2nd-frame baseline (alternate-frame flicker made the noise reference parity-dependent).
- **Result:** 1450 / 18000 unchanged; 15000 unchanged except cell #5 now ends 15789 (was 15768). Runner test window still 5 cells (misses faint 15459; finds the real 15776 cell). Details in `lab/failures.md`.
- **Next:** Sai to decide between running the full video now (QA with per-window kymograph PNGs) and first trying two offset window grids merged.

## 2026-11-16 — First whole-video run: 106 cells

- **Done (Sai: "run full video now"):** `py run_full_video.py --resume` over frames 0–37,274: 25 windows, ~40 min, no errors. Results in `C:\Users\amara\OneDrive\Desktop\RBC_full_video_results\` (trajectories.csv, biomechanics.csv, windows.csv, summary.json, 25 videos, 25 kymograph PNGs; ~1 GB).
- **Result:** 106 cells; 1–8 per 1500-frame window (median 4). Track length 295–461 frames (mean 363); every track reaches x ≤ 30 (none stops mid-channel).
- **Cross-check vs hand-validated windows:** 1450 window's 3 cells identical. 18000 window's 2 cells found (the 18001 cell enters at 17989, so it is owned by the previous window). 15000 window: all confirmed cells except the faint 15459 one (known placement-dependent miss), plus the real 15776 cell. #35/#36 end ~20 frames later in the exit funnel than the boxes Sai confirmed.
- **QA:** contact sheet of all 25 kymographs shows nearly every bright ridge tracked. The last window has a cell entering at the end of the recording (untrackable by design). Faint ridges need full-resolution review.
- **Open:** faint-cell placement dependence (≈ ±1 cell in busy windows); exit-funnel residual lag; nothing committed (git sync blocked).

## 2026-11-16 — Pipeline write-up (PDF)

- **Done (Sai asked):** Compiled `writeup/draft.md` from `lab/`, the code and the whole-video results, with 9 figures in `writeup/figs/`. Rendered to `writeup/RBC_tracking_pipeline.pdf` (14 pages; copy in `Desktop\RBC_full_video_results\`). Covers data/geometry, every stage with its parameters, the whole-video runner, validation, results with caveats (model-shaped speeds/durations, box-model shape numbers, unverified 10,000 fps vs CIHX frame range), the failure history, the sketch assessment, limitations, files and output columns. No physical speeds reported.
- **Evidence:** the PDF itself; build scripts kept in the session scratchpad (not project code).
- **2026-11-16 — Comparison video:** `C:\Users\amara\rbc_tracking\comparison_raw_vs_tracked_frames_14986-15900.mp4` (raw AVI frames on top, tracked video below, same 915 frames, H.264, 20 fps). Frame alignment verified by pixel comparison (offset 0 best at 4 sample frames).
- **2026-11-16 — Three-panel comparison video:** `C:\Users\amara\rbc_tracking\comparison_raw_enhanced_tracked_frames_14986-15900.mp4` — raw / enhanced / tracked, same 915 frames, H.264 1280×1512, 20 fps, 68 MB. The enhanced panel (viewing aid only, not used by the tracker) shows |frame − local median background (every 8th frame, 14700–16200)| after removing each frame's median offset (alternate-frame flicker), minus 1.5× each pixel's typical |diff| (suppresses tube-jitter walls), ×12 with a square-root curve, the bath dimmed to 12 %, speckle removed, and a faint channel outline. Bore cells visible as bright blobs; the faint 15459 cell barely shows (its per-frame signal is at noise level).

## 2026-11-18 — Offset grid (grid B) added to the whole-video runner

- **Done (Sai: "add the offset-grid pass"):** `run_full_video.py` now runs a second window grid with cores shifted by core/2 (750 frames), stored in `per_window_b/`, and merges it into grid A. Grid A is never changed. A grid-B track is a duplicate of an A track (median bore distance < 20 px over ≥ 20 shared frames → the A cell is labelled `A+B`), a conflict (shares > 30 % of its bore frames within 14 px, or crosses it → reported, not added), or a new cell (added, `source = B`). Cells are renumbered by entrance; `merge_report.csv` lists every grid-B track. Videos/kymographs re-render only where the drawn cells changed (`render_state.json`). `--no-offset-grid` turns it off. `tracking.py` unchanged.
- **Test (frames 13500–18000, scratchpad):** 13 cells; 10 found by both grids with identical paths (median distance 0.0 px), 0 conflicts, 0 additions. Grid B missed 2 cells grid A found (17118, 17989), so placement-dependent misses go both ways and the union helps. **Grid B did not recover the faint 15459 cell either:** both placements put the cell ahead at entrance 15350, and its entry-flare path is blanked ±20 px over the faint cell's start. This miss is a masking problem, not a placement one.
- **Next:** whole-video grid-B run in progress (reuses grid A). Candidate fix for the masking miss (needs approval, `tracking.py`): in the gap pass, blank accepted tracks only in the bore, not in the entry flare where DP paths are approximate.

## 2026-11-18 — Offset grid finished with ridge-evidence gate: 109 cells

- **Done:** Whole-video grid B (25 windows, 62 min incl. merge and rendering) merged into grid A. The first merge added 14 cells, 11 of them noise (see `lab/failures.md`), so grid-B additions now need ridge evidence ≥ 2.5 (`RIDGE_MIN_RATIO` in `run_full_video.py`). Re-merge: 3 added (entering 6628, 7260, 23898; evidence 3.11 / 3.43 / 2.82; round cell visible in the entry funnel), 11 rejected (1.07–1.91), 1 conflict (crosses cell 56), 92 duplicates.
- **Result:** **109 cells** = 92 found by both grids (`A+B`) + 14 grid A only + 3 added by grid B. IDs 1–109 contiguous, sorted by entrance; 25 videos + 25 kymographs (7 re-rendered).
- **Prototype rejected:** blanking accepted tracks only in the bore (for the faint 15459 cell) made the gap pass latch onto the next cell's entrance and lose recoveries (`lab/failures.md`).
- **Open:** (1) about 8–11 of the original grid-A tracks look like noise (empty entry funnel; clustered at 27700–29100 and 35200–35800); the ridge score cannot remove them because some bright real cells also score low. (2) The faint 15459 cell is still missed in both grid placements. (3) `writeup/` still describes the 106-cell single-grid run.

## 2026-11-18 — Entry-funnel cell check (flags noise tracks)

- **Done (Sai approved):** `run_full_video.py` now checks every final track for a real cell in the ENTRY funnel while the track is there (`EntryChecker`). Per grid-A window it builds a local background (median of every 8th frame) and jitter floor on an entry-funnel crop; per track it takes the strongest cell-sized blob (|diff| − frame median − 1.5 × jitter, blurred σ = 4 px) within ±25 px of the track on up to 8 funnel frames; score = median. Flags only: `entry_cell_score` and `entry_check` (`cell_seen` ≥ 2.0, `uncertain` 1.0–2.0, `no_cell_seen` < 1.0) in `biomechanics.csv`, counts in `summary.json`, crops of flagged tracks in `entry_check_review.png`. Nothing is removed.
- **Calibration (labels by eye from earlier crops):** 21 empty-funnel tracks 0.32–1.28 (median 0.59); 30 real cells 1.61–50 (median 5.4), including the bright cells the ridge score had rated 1.2–1.8.
- **Result on the 109 cells:** 90 `cell_seen`, 7 `uncertain`, 12 `no_cell_seen`; 0 contradictions with eye labels. All 12 `no_cell_seen` show an empty funnel on the review sheet (2 newly found: entering 26094, 28856). Uncertain by eye: real-looking 30487, 23898, 24084; noise-looking 7830, 12872, 11189; unclear 29610. Best estimate ≈ 93–94 real passages.
- **Evidence:** `lab/evidence/20261118_offset_grid/entry_check_review_109cells.png`, prototype scores.
- **Open:** flagged tracks are kept in all outputs; decide whether to exclude `no_cell_seen` from analyses. `writeup/` not yet updated.

## 2026-11-18 — Write-up updated (offset grid + entry check)

- **Done (Sai asked):** `writeup/draft.md` and `writeup/RBC_tracking_pipeline.pdf` (17 pages; copy in `Desktop\RBC_full_video_results\`) now cover the two window grids, the merge rules and ridge-evidence gate (§9.1), the entry-funnel noise flag (§9.2), the whole-video checks (§10.1, new Figures 10–11), the 109-cell result with statistics on the 90 `cell_seen` cells (§11), four new rows in the failure table (§12), and updated limitations and output columns. Pipeline diagram and statistics charts regenerated. No physical speeds reported.

## 2026-11-21 — Pushed to Asheesh10/RBC_tracking_final

- **Done (Sai asked):** README gained a "Run the tracker" section and file-table rows for the scripts. Created branch `rbc-final` from local `master` and committed all code (`tracking.py`, `run_full_video.py`, `zone13_tracking.py`, `scripts/sync.ps1`), lab notes and evidence, and `writeup/`. Pushed as `main` to the new repository https://github.com/Asheesh10/RBC_tracking_final (remote `final`). Nothing was pushed to Sai2125/iisc-rbc-tracking.
- **Left out on purpose:** the AVI (git-ignored), `kymo_dump.npy`, `kymo_crop.png`, `tracking_pre_fix_backup.py` (scratch), `rbc-tracking/` (a nested old clone), `.vscode/` (git-ignored, machine-specific).
- **Open:** local `master` is unchanged and still diverged from `origin/master` (1 ahead, 22 behind). The new repository's history is local `master` + this commit, so Sai2125's 22 intermediate commits are not in it as separate commits (their final code is, inside `tracking.py`).

## 2026-11-23 — Literature review Section 6 corrected

- **Done (Sai asked):** Rewrote Section 6 of `Desktop\RBC_Tracking_Literature_Review.txt` (outside the repo; original kept as `RBC_Tracking_Literature_Review.before_sec6_rewrite.txt`). The old text described the planned per-frame design (bilateral filter, watershed, Kalman + Hungarian matching), which lives only in `zone13_tracking.py`. The new text describes the working pipeline: whole-passage kymograph DP, zone speed limits, motion × contrast energy, random-line and entry-funnel checks against noise, Kalman + RTS smoothing. It adds an honest "not yet better" list (no run-twice test, hand-fitted settings, box-model shape numbers, no watershed/bilateral in the tracker, no ground-truth accuracy) and a corrected summary table. Sections 1–5 unchanged.

## 2026-11-23 — Cell outlines (`cell_outlines.py`)

- **Done (Sai: "do the cell outlines"):** new `cell_outlines.py` traces each cell's real outline in the entry and exit funnels and measures area, perimeter, circularity, equivalent-ellipse axes, aspect ratio, deformation index, equivalent diameter and solidity per frame, with quality flags (cut off, wall contact, touching neighbour, low contrast, irregular, size outlier). Method: per-block background built from frames aligned by phase correlation (the tube image wobbles up to ±2 px, which caused the wall flicker), each frame aligned before subtracting, Hungarian assignment of blobs to tracks, outline = darkest line of the cell's rim found by a closed polar dynamic programme. Exit outlines are found without the tracks and matched to tracks in passage order. `tracking.py` and `run_full_video.py` unchanged. README updated.
- **Checks:** synthetic cells with the measured rim profile (`--self-test`): area independent of contrast (±3 %), circularity 0.996 (0.978 faint) for a circle, AR 1.42 / 1.82 for true 1.5 / 2.0. Window 14900–16500: 6/6 cells with 24–76 good entry outlines, review sheet checked by eye.
- **Whole video (33 min):** 82/90 `cell_seen` tracks have ≥ 5 good entry outlines vs 1/12 `no_cell_seen`. Areas split into small particles (34–73 px², 24 tracks, incl. all 12 `no_cell_seen`) and RBCs (135–299 px², 85 tracks) → best RBC count 85. Resting RBC shape: circularity 0.990, AR 1.05, deformation index 0.024, equivalent diameter 16.4 px (medians, 78 RBC-sized cells with ≥ 5 good frames); entry-taper AR p90 1.36. Absolute area depends on focus (area vs contrast Spearman −0.87; faint cells read up to ~2× larger). Exit outlines for 34 tracks. Outputs `cell_outlines_per_frame.csv`, `cell_outlines_per_cell.csv`, `cell_outline_contours.npz`, `cell_outline_review_01–10.png` in `Desktop\RBC_full_video_results`.
- **Found:** tracker bug — cell #37's track slides onto #38's ridge in the lower bore (frames 15200–15330), so exit frames/transit times of both are off by ~10–25 frames (`lab/failures.md`).
- **Evidence:** `lab/evidence/20261123_cell_outlines/`.
- **Git:** Sai approved; `cell_outlines.py`, README, lab notes and evidence committed on branch `rbc-final` and pushed to remote `final` (`main`). Sync pull still fails (`rbc-final` has no upstream); `master`/`origin` untouched.
