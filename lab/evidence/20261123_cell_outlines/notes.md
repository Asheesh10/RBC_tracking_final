# 2026-11-23 — Cell outlines (`cell_outlines.py`)

## Radial profile of entry-funnel cells (background difference, grey levels)

Mean over 1-px rings around the outline centre, r = 0.5, 1.5, … px:

| cell / frame | profile (r = 0.5 … 11.5 px) |
|---|---|
| #38 f15097 (bright) | +64 +70 +75 +57 +39 +17 −14 −33 −33 −20 −8 −2 |
| #41 f15545 | +37 +38 +42 +39 +26 +7 −11 −19 −18 −12 −5 −2 |
| #40 f15420 | +71 +67 +59 +43 +26 +2 −19 −28 −24 −14 −4 +1 |
| #37 f15013 | +13 +17 +24 +29 +25 +19 +5 −10 −14 −13 −7 −4 |
| #42 f15799 (faint) | +4 +4 +8 +9 +9 +7 +2 −2 −4 −4 −4 −2 |

Every cell: bright core to r ≈ 5, zero at r ≈ 6, darkest rim at r ≈ 7.5–8.5, background at r ≈ 11.
A fixed noise threshold traced the faint outer tail (equivalent diameter 20 px for bright cells,
15 px for faint #42), so the outline follows the rim minimum instead.

## Synthetic check (`py cell_outlines.py --self-test`)

Measured profile, noise σ = 3.5, ellipses at 0/30/60°, 9 runs each:

| true AR | contrast | area (rim) | AR measured | circularity true / measured |
|---|---|---|---|---|
| 1.0 | full | 203 ± 2 | 1.02 | 1.000 / 0.996 |
| 1.2 | full | 204 ± 2 | 1.17 | 0.988 / 0.987 |
| 1.5 | full | 205 ± 1 | 1.42 | 0.941 / 0.951 |
| 2.0 | full | 205 ± 1 | 1.82 | 0.841 / 0.875 |
| 1.0 | 1/4 | 205 ± 5 | 1.06 | 1.000 / 0.978 |
| 1.2 | 1/4 | 205 ± 4 | 1.17 | 0.988 / 0.966 |
| 1.5 | 1/4 | 204 ± 7 | 1.40 | 0.941 / 0.933 |
| 2.0 | 1/4 | 209 ± 7 | 1.76 | 0.841 / 0.868 |

Area does not depend on contrast (±3 %). Elongated cells read a little rounder than they are
(AR 2.0 → 1.8, circularity 0.84 → 0.87) because the rim has a fixed width. Faint round cells read
slightly less circular (0.978) because noise roughens the rim. Smoothing settings were chosen from a
sweep (more smoothing fixes faint-cell circularity but lowers AR of elongated cells).

## Tube wobble

Phase correlation of funnel crops against the background, frames 15330–15387: the image shifts by up
to ±1.4 px (exit), ±2.2 px (entry), ±0.8 px (bore), mostly along x, period ≈ 60 frames. Different
amounts per region, so it is not one rigid shift. Aligning each crop to its background before
subtracting removes the bright/dark wall flicker (`exit_funnel_plain_vs_aligned_difference.png`).
Backgrounds are now built from aligned frames too (sharper walls).

## Exit funnel identity (cells 37 / 38)

`kymo_15150_15345_track37_slides_onto_38.png`: smoothed energy kymograph (bottom) and with tracks
(top, red #37, green #38). The faint #37 ridge continues to the exit at ~15290 (`exit_funnel_15228_15333_cell_emerging.png`:
a cell emerges at ~15285 and becomes the left ring). Track #37 leaves its ridge at x ≈ 510 and ends
on #38's bright ridge; #38's track runs 10–20 frames late. So the "9"-shaped ring that track #37
points at in the exit funnel is really #38. `cell_outlines.py` matches exit tracklets to tracks in
passage order and gives the "9" ring to #38.

## Window result (cells entering 14900–16500)

`window_15000_16500_outline_review.png` (green good, yellow usable, red poor):

| cell | good entry frames | rest area px² | eq. diam px | circularity | AR | taper AR p90 | exit |
|---|---|---|---|---|---|---|---|
| 37 | 45 | 245 | 17.7 | 0.992 | 1.03 | 1.84 | none found |
| 38 | 42 | 199 | 15.9 | 0.994 | 1.08 | 1.65 | ring, A 42, n 47 |
| 39 | 24 | 219 | 16.7 | 0.991 | 1.03 | 1.53 | ring, A 44, n 49 |
| 40 | 76 | 188 | 15.5 | 0.990 | 1.07 | 1.46 | none found |
| 41 | 43 | 201 | 16.0 | 0.997 | 1.03 | 1.30 | none found |
| 42 | 47 | 265 | 18.4 | 0.970 | 1.08 | 1.25 | none found |

Faint cells (#37 contrast 8, #42 contrast 5) read larger than bright ones (#38 35, #40 27):
probably focus (a defocused cell has a wider, fainter rim), not the tracer (synthetic area is
contrast-independent). Absolute area therefore depends on focus.

## Whole video (109 tracks, `Desktop\RBC_full_video_results`, 33 min)

- 7,030 entry-funnel outlines (3,364 good) and 1,891 exit-funnel outlines.
- Agreement with the entry check: 82 of 90 `cell_seen` tracks have ≥ 5 good entry outlines; only 1 of 12 `no_cell_seen` does (and it is a small particle, below).
- **Size classes** (`size_class`): entry outline areas split with a gap between 73 and 135 px².
  24 tracks are `small_particle` (34–73 px²: a bright speck 4–5 px across with no dark rim,
  `small_particle_tracks_review.png`): all 12 `no_cell_seen`, 5 of 7 `uncertain`, 7 `cell_seen`.
  85 tracks are `rbc_sized` (83 `cell_seen` + `uncertain` #63, #87, which are very defocused but real).
  So the best count of RBC passages is 85.
- **Resting shape** (78 `rbc_sized` cells with ≥ 5 good frames, all `cell_seen`): circularity median
  0.990 (IQR 0.983–0.996, range 0.952–0.999), aspect ratio 1.048 (1.032–1.067, max 1.40), deformation
  index 0.024 (0.016–0.032), equivalent diameter 16.4 px (15.6–17.7), area 211 px² (190–247, range
  135–299). Frame-to-frame area spread within a cell: median IQR 2.7 % of area. The other 7 RBC-sized
  tracks have < 5 good frames (faint/defocused or crowded: #9, #45, #59, #63, #66, #87, #100).
- **Entry taper** (70 cells): aspect ratio p90 median 1.36 (1.28–1.51), deformation index p90 0.15 (0.12–0.20).
- **Focus:** area vs contrast Spearman −0.87; area by contrast tercile 264 / 209 / 183 px² (faint → large).
  Sharply focused cells (#6, #8, #109) read 135–144 px²; faint defocused discs 250–300 px²
  (`focus_bright_vs_faint_review.png`). Elongated in-focus cells (#69, #99, AR 1.4) look genuinely elliptical.
- **Box model for the same 78 cells:** `area_min` 66.76 for 44 cells, `circularity_min` 0.50 for 44,
  `circularity_mean` has 11 distinct values. The outline numbers have 78 distinct values.
- **Exit funnel:** 36 cell-like exit tracklets; 34 tracks matched (30 of them RBC-sized; 22 clear,
  12 with another track reaching the exit within 40 frames). Exit ring: equivalent diameter 8.0 px (7.5–8.9), circularity 0.93,
  AR 1.26, median 40 good frames per matched cell. 75 tracks have no exit outline (faint).
