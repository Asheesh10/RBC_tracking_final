# Current

Updated: 2026-11-23

## Now

Stage 3 tracker (`tracking.py`: kymograph DP + gap pass + re-trace + RTS smoother) tracks **6/6 cells** on frames 15000–16500, 3 on 1450–2950 and 2 on 18000–19500. Sai confirmed cells #1 and #3 visually on 15000–16500 (2026-11-15).

Today's changes (all Sai-approved, all in `tracking.py`, **not committed**):

1. Gap pass is now a kymograph DP with fit-based validation. It recovers the faint 15459 cell and removes a ghost track stuck on debris at x≈731.
2. Video is written as H.264 via ffmpeg (fallback OpenCV `avc1` → `mp4v`).
3. Every main-pass track is re-traced with neighbours masked, a faster exit-funnel model and a noise-biased exit rule. This fixes cell #1 riding cell #2, cell #3 ending early, and the 1450 cell #2 lag.

Repo: private GitHub [Sai2125/iisc-rbc-tracking](https://github.com/Sai2125/iisc-rbc-tracking). AVI is gitignored and lives at `C:\Users\amara\Downloads\`. `AGENTS.md` is the rule file.

**Whole video (2026-11-18):** `run_full_video.py` with two offset window grids tracked **109 cells** over frames 0–37,274 (92 found by both grids, 14 grid A only, 3 added by grid B after a ridge-evidence check). Results on the Desktop in `RBC_full_video_results`. An entry-funnel check flags each cell: 90 `cell_seen`, 7 `uncertain`, 12 `no_cell_seen` (likely noise; kept but flagged). The faint 15459 cell is still missed.

**Cell outlines (2026-11-23):** `cell_outlines.py` traces real outlines in the entry/exit funnels (dark-rim tracing on wobble-aligned background differences) and writes shape numbers per frame and per cell into `RBC_full_video_results`. 85 of the 109 tracks are RBC-sized; 24 are small particles (34–73 px²; includes all 12 `no_cell_seen`). 78 RBCs have ≥ 5 good entry outlines: resting circularity 0.990, AR 1.05, equivalent diameter 16.4 px (medians). Absolute area depends on focus. Exit outlines for 34 tracks. Committed on `rbc-final` and pushed to `final` (Asheesh10/RBC_tracking_final, `main`) on 2026-11-23.

## Next

1. Tracker bug from the outlines work: a faint cell's track can slide onto the next cell's ridge in the lower bore (#37 → #38, frames 15200–15330), corrupting exit frames/transit times. Measure how often it happens before deciding on a fix in `tracking.py` (needs approval).
2. Bore shape (length/elongation in the 7-px bore) is not measured yet; the outlines cover only the funnels.
3. Residual tracking issues: box #4 (15000) runs ~10 px ahead in its last frames; boxes trail the cell centre by ~5–19 px in the exit funnel; cell #3 ends ~20 frames before the cell fully leaves.
4. Structural (see `lab/failures.md`): faint cells can still be hidden by a neighbour's entry-funnel blanking (15459). Cells already in the bore at a window's start are untracked by design. Longer term: exit-side detection so an exit with no matching track is flagged.

## Blockers / open

- **Git:** `scripts/sync.ps1 pull` fails. Local `master` is 1 ahead / 22 behind `origin/master`, and `tracking.py` has large uncommitted changes (earlier work plus today's). Needs Sai's decision; agents must not rebase/merge ad hoc.
- `tracking.py` header comment says 1450 cell #2 exits at 2730; frames show ~2675.
- Exit-funnel walls flicker in background-difference images (tube jitter); per-frame funnel detection is not usable as is.
- AVI header count 37,274 still not reconciled with CIHX saved range 0–17,471.
- No µm/pixel calibration. No physical speeds.

## Key paths

| What | Path |
|---|---|
| Agent rules (all harnesses) | `AGENTS.md` |
| Human setup + map | `README.md` |
| Git sync | `scripts/sync.ps1`, `scripts/sync.sh` |
| Stage 3 tracker | `tracking.py` (run with `py`; `RBC_DEBUG_ENTRANCES=1` for decisions) |
| Zone 1/3 prototype | `zone13_tracking.py` |
| Related work | `lab/related-work.md` |
| Failures | `lab/failures.md` |
| Evidence (2026-11-15 tracker fixes) | `lab/evidence/20261115_gap_dp_fix/` |
| Stage 1 outputs + Stage 3 runs (local) | `C:\Users\amara\rbc_tracking\output\` (`stage3/` = latest run; `stage3_*` = saved comparisons) |
| AVI (local only) | `C:\Users\amara\Downloads\2A002 2_C001H001S0001-001.avi` |
| Tube-mask script | `lab/tools/probe_tube_mask.py` |
| Annotated screenshot (do not use for detection) | `image (3).png` |
