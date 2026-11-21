# Decisions

Add a line when a choice should outlive a chat recap. Do not delete old items; mark them superseded.

- **2026-09-16 — Method class:** Algorithmic detect and measure. Do not dump the video into an LLM/VLM to explain what is happening.
- **2026-09-16 — Next scientific step:** Single-frame identification of the cell **in the tube**, not whole-field RBC detection.
- **2026-09-16 — Scope of this lab notebook:** Lightweight files under `lab/`. No extra platform. Raw AVI stays in the project root.
- **2026-09-16 — Analysis data:** Original AVI / decoded pixels / lossless PNG. Lossy previews are for viewing only. Do not overwrite raw files.
- **2026-09-16 — Speeds:** No physical speed until capture timing and spatial calibration are verified. Pixel motion can be discussed after tracking works.
- **2026-09-16 — Annotation:** Red circle on `image (3).png` is not part of the experiment image. Do not use it as a detection feature. Screenshot coordinates may not match AVI coordinates.
- **2026-09-16 — Timing fields:** Prefer CIHX `recordInfo/recordRate` / MII film speed (10000). Do not use `deviceInfo/recordRate` (125) or AVI playback FPS as acquisition timing without verification.
- **2026-09-16 — Code permission:** Agents must not create or modify code until Sai gives permission. Propose first, wait, then implement. Rule lives in `AGENTS.md`.
- **2026-09-16 — Tube-mask prototype:** Sai approved seeking one random AVI still (bypassing the red-circle PNG) and building a lumen mask. Seed 20260916 → frame 24642 / 37274. Script: `lab/tools/probe_tube_mask.py`.
- **2026-09-16 — Log failures:** Misses go in `lab/failures.md` (what we saw, why, what changed, what remains). Do not keep only the successful overlay.
- **2026-10-05 — Git-shared notebook:** Private GitHub repo; AVI never committed. `AGENTS.md` is canonical for Cursor and other harnesses. Lab notes push via `scripts/sync.ps1 push-lab` / `sync.sh`; analysis code still needs Sai’s permission.
- **2026-11-15 — Gap-pass recovery:** Faint cells missed by the main kymograph pass are recovered with a kymograph DP (time-integrated), not a per-frame follower. Per-frame signal of such cells is at noise level; the per-frame follower drifted and locked onto debris (x≈731). Accept a recovery only if a 3-parameter fit along it scores ≥2.8× random-line noise. Sai approved.
- **2026-11-15 — Video encoding:** `tracking.py` writes H.264 via ffmpeg (CRF 18) so outputs play in browsers and editor viewers; OpenCV `mp4v` only as a last-resort fallback. Sai approved.
- **2026-11-15 — Shared-ridge re-trace:** Only tracks whose DP paths share a ridge (within ±14 px for ≥5 frames) are re-traced, with neighbours masked and a faster exit-funnel model. Main-pass extraction keeps the original speed model so accept/reject decisions do not move. Re-tracing every track was rejected because it moved correct exits.
- **2026-11-15 — Exit frame rule:** In the re-trace, a passage's exit frame is chosen by total evidence above the noise level, not by average score + left-exit bonus (which favoured short paths and neighbours' exits). Every main-pass track is re-traced this way; the main pass keeps the original rule so acceptance does not change. Supersedes the shared-ridge-only trigger from earlier the same day. Sai approved.
- **2026-11-18 — Offset grid:** Whole-video runs use two window grids (cores at 0, 1500, … and at 750, 2250, …). Grid A is the base; grid-B tracks are only added when they match no grid-A track and do not conflict with one. Cells found by both grids are labelled `A+B`. Sai approved.
- **2026-11-18 — Noise-track flag:** Each whole-video track gets an entry-funnel cell check (`cell_seen` / `uncertain` / `no_cell_seen`). Tracks are flagged, never deleted, so analyses can filter on `entry_check`. Sai approved.
