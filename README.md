# IIsc RBC tracking

Track red blood cells in a high-speed channel video. Algorithmic detect-and-measure — not dumping the video into an LLM.

The ~23 GB AVI is **not** in git. Copy it next to this clone as
`2A002 2_C001H001S0001-001.avi`.

## Setup

1. `git clone` this repo and open the folder in Cursor, VS Code, or another harness.
2. Place the AVI (keep `.cihx` / `.mii` if they are not already in the clone).
3. At the start of a session, run the pull script (agents should do this without being asked):

```powershell
.\scripts\sync.ps1 pull
```

```bash
./scripts/sync.sh pull
```

Agent rules for every harness: read [`AGENTS.md`](AGENTS.md). Cursor skills under `.cursor/skills/` are overlays only.

## Run the tracker

Needs Python with `opencv-python`, `numpy` and `scipy` (on the original laptop: the `py` launcher, Python 3.14). `ffmpeg` on PATH is optional but gives small H.264 videos.

1. **Stage 1 (once):** `py preprocessing.py` writes `calibration.json` and `median_bg.npy` (median background, channel line, widths).
2. **One window (for checking):** `py tracking.py --start-frame 15000 --max-frames 1500` writes `output/stage3/` (trajectories, biomechanics, tracked video). Add `RBC_DEBUG_ENTRANCES=1` to print every accept/reject decision.
3. **Whole video:** `py run_full_video.py --resume` runs both window grids, merges them, flags likely-noise tracks and writes everything to `Desktop\RBC_full_video_results\`. Re-running continues where it stopped.
4. **Cell outlines (after step 3):** `py cell_outlines.py` traces each cell's outline in the entry and exit funnels and writes shape numbers (area, perimeter, circularity, aspect ratio, deformation index) per frame and per cell, plus review sheets, into the same results folder. `--from-frame/--to-frame` limit it to cells entering in a range; `--self-test` runs the synthetic accuracy check.

In VS Code: open this folder, then Run and Debug → "RBC - test window 15000-16500" or "RBC - full video (resumes if stopped)" → Ctrl+F5. Those two configurations live in `.vscode/`, which is git-ignored because it holds machine-specific paths, so on a new machine create them again or use the commands above.

How it all works, with results and known limits: [`writeup/RBC_tracking_pipeline.pdf`](writeup/RBC_tracking_pipeline.pdf).

## What each piece is for

| Path | Role |
|------|------|
| `AGENTS.md` | Rules for agents (permission, no full-video, logging, git pull). Source of truth. |
| `README.md` | Humans: setup and this map. |
| `preprocessing.py` | Stage 1: median background, channel centreline, width profile (`calibration.json`). |
| `tracking.py` | Stage 3 tracker for one window: kymograph, main pass, gap pass, re-trace, Kalman + RTS smoother, H.264 video. |
| `run_full_video.py` | Whole video: two offset window grids, merge with ridge-evidence check, entry-funnel noise flag, resume. |
| `cell_outlines.py` | Real cell outlines in the entry/exit funnels (dark-rim tracing on wobble-aligned background differences) and shape numbers from them; replaces the box-model shape numbers for analysis. |
| `zone13_tracking.py` | Prototype: funnel detection + Kalman + Hungarian matching (not used by the main pipeline). |
| `lab/CURRENT.md` | One page: now / next / blockers. Read this first. |
| `lab/log.md` | Dated history of what we did. |
| `lab/decisions.md` | Choices that should stick. |
| `lab/failures.md` | What went wrong, why, what we changed. Source of the “why this method”. |
| `lab/gameplan.md` | Current detect/measure plan. |
| `lab/related-work.md` | Papers/tools we actually use or reject. |
| `lab/evidence/` | Small proof images/notes. Never the AVI. |
| `lab/tools/` | Experiment scripts (e.g. tube mask). Not auto-pushed. |
| `.cursor/skills/rbc-lab-log/` | Cursor: keep `lab/` updated. |
| `.cursor/skills/rbc-git-sync/` | Cursor: call the sync script. |
| `.cursor/skills/rbc-writeup/` | Cursor: on request, draft `writeup/draft.md` from `lab/`. |
| `scripts/sync.ps1` / `scripts/sync.sh` | Deterministic `git pull --ff-only`; optional lab-only push. |
| `writeup/` | Pipeline description: `draft.md` (source), `figs/`, `RBC_tracking_pipeline.pdf`. |
| `.github/copilot-instructions.md`, `CLAUDE.md`, `.clinerules` | One-liners: follow `AGENTS.md`. |
