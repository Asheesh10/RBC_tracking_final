"""
Cell outlines in the entry and exit funnels, and shape numbers measured from them.

Why this exists
---------------
The tracker (tracking.py) measures shape from a box model, which is not a real outline:
the box area takes the same few values for most cells. This script traces each cell's real
outline in the two wide funnels, where a cell is free (or nearly free) of the walls.

Method (per frame, per cell)
----------------------------
1. Background difference. Per block of 1500 frames, a local background (median of every 8th
   frame over the block +-300 frames) and a per-pixel jitter floor (median |difference|,
   which is high on the flickering tube walls) are built for an entry crop (x >= 1100) and an
   exit crop (x <= 230). d = frame - background, minus its median inside the channel
   (alternate-frame flicker), smoothed with a 0.8 px Gaussian.
2. Cell centre. Blob map = |d| minus 1.5 x jitter floor, blurred to cell size (sigma 3 px).
   Local maxima inside the channel are candidates. Each cell present in that frame is assigned
   at most one candidate (Hungarian assignment on distance to its predicted position: the
   previous outline centre if there is one, else the track point), within a gate.
3. Outline. A bright-field RBC here is a bright core with a dark rim (radial profile: bright
   to r~5 px, zero at r~6, darkest at r~7.5, background at r~11). The outline is traced along
   the darkest line of the rim with a closed dynamic programme in polar coordinates around the
   centre (72 angles, radius 3-15 px, max 1 px radial change per 5 degrees, curvature
   penalty), refined to sub-pixel, then lightly smoothed. The rim may not cross the channel
   walls. The centre is moved to the outline's centroid and the outline traced again.
   Why the rim and not a threshold: a noise threshold lands on the faint outer tail of the
   rim, so the outline (and area) shrinks for faint cells. The rim minimum does not move with
   contrast (synthetic test: area changes ~3% between full and 1/4 contrast).
4. Numbers from the outline polygon: area, perimeter, circularity 4 pi A / P^2, equivalent
   ellipse from second moments (major, minor, angle), aspect ratio major/minor, deformation
   index (major-minor)/(major+minor), equivalent diameter, solidity (area / convex area).
   All in pixels: there is no um/pixel calibration yet.

Synthetic check (run `py cell_outlines.py --self-test`): with the measured rim profile and
the real noise level, circularity of a circle reads 0.996 at full contrast and 0.978 at 1/4
contrast; aspect ratio reads 1.02 for a circle, 1.42 for 1.5 and 1.82 for 2.0 (the fixed-width
rim makes elongated cells read a little rounder: circularity 0.875 for a true 0.841). Area does
not change with contrast (within ~3%).

Quality flags (per frame)
-------------------------
cut_off            outline reaches the image edge (cell entering at the right edge)
wall_contact       fraction of the outline within 1.5 px of a wall (> 0.15 flagged)
touching_neighbour another cell's outline centre closer than the two radii + 5 px
low_contrast       (core - rim) / noise < MIN_CONTRAST
irregular          solidity < 0.90
size_outlier       area differs > 35% from the cell's median entry area (set afterwards)
quality = good     wide part of the funnel, none of the flags
          usable   contrast ok, not cut off, not touching a neighbour (taper / wall contact allowed)
          poor     everything else

Per cell
--------
rest_*      medians over `good` entry-funnel frames (wide part): the cell's resting shape
taper_*     90th percentile over good/usable frames in the entry taper (cell squeezing in)
exit_*      medians over good (else usable) exit-funnel frames; see "Exit funnel" below
size_class  small_particle if the entry outline area is < 100 px^2. In the whole video the areas
            split cleanly: 34-73 px^2 are small bright specks with no dark rim (debris, possibly
            platelets, or noise), 135-299 px^2 are RBCs. Not an RBC if small_particle.
Absolute area depends on focus: sharply focused cells read ~135-180 px^2, defocused (faint)
cells ~250-300 px^2 (Spearman area vs contrast -0.87 over 78 RBCs), while the synthetic check
shows the tracer itself is contrast-independent. Prefer circularity / aspect ratio /
deformation index for comparisons between cells until focus is controlled or calibrated.

Exit funnel
-----------
The exit end of a track is not always on the right cell (a faint cell's track can slide onto
the next cell in the lower bore: cells 37/38 at frames 15200-15330). Exit outlines are
therefore found without the tracks: blobs in the exit funnel are linked over frames into exit
tracklets, and tracklets are matched to tracks in order (cells cannot overtake in the bore)
by the time each reaches the bore exit. Exit numbers carry their own quality and an
`exit_match` label. In the exit funnel a cell shows as a small dark ring (~8 px across even for
sharply focused cells that are ~13 px at entry), so exit outlines are a different visible
feature from entry outlines: compare exit numbers with each other, not with entry numbers.

Outputs (in --results)
----------------------
cell_outlines_per_frame.csv   one row per cell per funnel frame
cell_outlines_per_cell.csv    one row per cell: resting shape (entry funnel, wide part),
                              entry-taper deformation, exit shape
cell_outline_contours.npz     outline points (72 per outline) for overlays
cell_outline_review_NN.png    review sheets: outlines drawn on background-difference crops
"""

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

import cv2
import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment

import tracking as T

# ------------------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------------------

BLOCK = 1500                 # background block length (frames)
BG_PAD = 300                 # background block padding (frames)
BG_STEP = 8                  # every 8th frame goes into the background median
ENTRY_X0 = 1100              # entry crop: x >= this
EXIT_X1 = 230                # exit crop: x < this
FUNNEL_MIN_WIDTH = 13.0      # px between wall lines; narrower = taper end / bore: no outline
WIDE_MIN_WIDTH = 26.0        # wide part of a funnel (cell can be free of the walls)
ENTRY_GATE = 16.0            # px, first association of a track with a blob (entry)
CONT_GATE = 7.0              # px, extra gate around the previous outline centre
MIN_CONTRAST = 4.0           # rim depth / noise for a usable outline
EXIT_LINK_CONTRAST = 2.0     # exit funnel: weakest detection that may extend a tracklet
EXIT_TRACKLET_CONTRAST = 3.0 # exit funnel: median contrast a tracklet needs to count as a cell
EXIT_BLOB_MIN = 0.7          # exit funnel: blob strength / noise worth tracing (speed)
WALL_FLAG = 0.15
TOUCH_MARGIN = 5.0           # px beyond the two rim radii: the rims' dark halos overlap inside this
SOLIDITY_FLAG = 0.90
SIZE_FLAG = 0.35
SMALL_OBJECT_AREA = 100.0    # px^2: whole video has a gap between 73 (small bright specks, no rim)
                             # and 135 (smallest in-focus RBC); below this = small_particle

N_ANG = 72
R_MIN, R_MAX, R_STEP = 3.0, 15.0, 0.5
MAX_DR = 2                   # radial samples per angular step (1 px per 5 degrees)
SMOOTH_PEN = 0.5             # curvature penalty per (radial sample)^2
R_SMOOTH = 1.5               # gaussian smoothing of r(theta), angular samples

FRAME_FIELDS = ["track_id", "frame", "funnel", "zone", "cx", "cy", "track_cx", "track_cy",
                "offset_px", "area_px2", "perimeter_px", "circularity", "major_px", "minor_px",
                "aspect_ratio", "deformation_index", "eq_diameter_px", "solidity", "angle_deg",
                "contrast", "wall_contact", "cut_off", "touching_neighbour", "quality", "flags",
                "exit_tracklet"]
CELL_FIELDS = ["track_id", "entrance_frame", "entry_check",
               "n_entry", "n_entry_good",
               "rest_area_px2", "rest_area_iqr_px2", "rest_eq_diameter_px", "rest_perimeter_px",
               "rest_circularity", "rest_aspect_ratio", "rest_deformation_index", "rest_contrast",
               "taper_n", "taper_aspect_ratio_p90", "taper_deformation_index_p90",
               "exit_match", "exit_tracklet", "exit_first_frame", "n_exit", "n_exit_good",
               "exit_area_px2", "exit_eq_diameter_px", "exit_circularity", "exit_aspect_ratio",
               "exit_deformation_index", "exit_contrast", "outline_status", "size_class"]


# ------------------------------------------------------------------------------------------
# outline tracing
# ------------------------------------------------------------------------------------------

def _polar_sample(img, cx, cy, cval=0.0, order=1):
    th = np.arange(N_ANG) * 2 * np.pi / N_ANG
    rr = np.arange(R_MIN, R_MAX + 1e-6, R_STEP)
    X = cx + rr[None, :] * np.cos(th)[:, None]
    Y = cy + rr[None, :] * np.sin(th)[:, None]
    P = ndi.map_coordinates(img, [Y.ravel(), X.ravel()], order=order, mode="constant", cval=cval)
    return P.reshape(X.shape), th, X, Y


def _closed_dp(C):
    """Minimum-cost closed path j(i), i = 0..n-1, |j(i) - j(i-1)| <= MAX_DR, j(n-1) ~ j(0).
    One DP per possible start radius, vectorised over start radii."""
    n, M = C.shape
    INF = 1e18
    cost = np.full((M, M), INF)
    cost[np.arange(M), np.arange(M)] = C[0]
    back = np.zeros((n, M, M), np.int8)
    for i in range(1, n):
        best = np.full((M, M), INF)
        arg = np.zeros((M, M), np.int8)
        for s in range(-MAX_DR, MAX_DR + 1):
            sh = np.full((M, M), INF)          # sh[:, j] = cost[:, j - s]
            if s > 0:
                sh[:, s:] = cost[:, :M - s]
            elif s < 0:
                sh[:, :M + s] = cost[:, -s:]
            else:
                sh[:] = cost
            sh += SMOOTH_PEN * s * s
            better = sh < best
            best[better] = sh[better]
            arg[better] = s
        cost = best + C[i][None, :]
        back[i] = arg
    tot = np.full(M, INF)
    endj = np.zeros(M, int)
    for s0 in range(M):
        lo, hi = max(0, s0 - MAX_DR), min(M, s0 + MAX_DR + 1)
        k = lo + int(np.argmin(cost[s0, lo:hi]))
        tot[s0] = cost[s0, k] + SMOOTH_PEN * (k - s0) ** 2
        endj[s0] = k
    s0 = int(np.argmin(tot))
    path = np.zeros(n, int)
    path[-1] = endj[s0]
    for i in range(n - 1, 0, -1):
        path[i - 1] = path[i] - back[i, s0, path[i]]
    return path


def trace_rim(ds, cx, cy, forbid):
    """Outline along the dark rim around (cx, cy). ds: smoothed signed background difference.
    forbid: True where the outline may not go (outside the walls / outside the image)."""
    P, th, X, Y = _polar_sample(ds, cx, cy)
    F = _polar_sample(forbid.astype(np.float32), cx, cy, cval=1.0, order=0)[0] > 0.5
    F = np.maximum.accumulate(F, axis=1)     # beyond a wall hit, everything is forbidden
    C = P.copy()
    C[F] = 1e6
    path = _closed_dp(C)
    j = path.astype(float)
    for i, k in enumerate(path):             # sub-pixel: parabola through the minimum
        if 0 < k < C.shape[1] - 1:
            a, b, c = C[i, k - 1], C[i, k], C[i, k + 1]
            den = a - 2 * b + c
            if den > 1e-6 and max(a, b, c) < 1e5:
                j[i] = k + float(np.clip(0.5 * (a - c) / den, -0.5, 0.5))
    r = ndi.gaussian_filter1d(R_MIN + j * R_STEP, R_SMOOTH, mode="wrap")
    xs, ys = cx + r * np.cos(th), cy + r * np.sin(th)
    rim_vals = P[np.arange(N_ANG), path]
    # wall contact: outline point with a forbidden pixel within 1.5 px outward
    probe = ndi.map_coordinates(forbid.astype(np.float32), [ys + 1.5 * np.sin(th), xs + 1.5 * np.cos(th)],
                                order=0, mode="constant", cval=1.0) > 0.5
    return xs, ys, r, rim_vals, probe


def shape_measures(xs, ys):
    x, y = np.asarray(xs, float), np.asarray(ys, float)
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A = 0.5 * cr.sum()
    sgn = 1.0 if A >= 0 else -1.0
    A = abs(A)
    if A < 1e-6:
        return None
    P = float(np.sum(np.hypot(x1 - x, y1 - y)))
    cx = sgn * np.sum((x + x1) * cr) / (6 * A)
    cy = sgn * np.sum((y + y1) * cr) / (6 * A)
    Ixx = sgn * np.sum((x * x + x * x1 + x1 * x1) * cr) / 12 - A * cx * cx
    Iyy = sgn * np.sum((y * y + y * y1 + y1 * y1) * cr) / 12 - A * cy * cy
    Ixy = sgn * np.sum((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * cr) / 24 - A * cx * cy
    ev, evec = np.linalg.eigh(np.array([[Ixx, Ixy], [Ixy, Iyy]]) / A)
    ev = np.maximum(ev, 1e-9)
    L, Wd = 4 * np.sqrt(ev[1]), 4 * np.sqrt(ev[0])
    hull = cv2.convexHull(np.stack([x, y], 1).astype(np.float32))
    return dict(cx=float(cx), cy=float(cy), area_px2=float(A), perimeter_px=P,
                circularity=float(4 * np.pi * A / P ** 2), major_px=float(L), minor_px=float(Wd),
                aspect_ratio=float(L / Wd), deformation_index=float((L - Wd) / (L + Wd)),
                eq_diameter_px=float(2 * np.sqrt(A / np.pi)),
                solidity=float(A / max(cv2.contourArea(hull), 1e-6)),
                angle_deg=float(np.degrees(np.arctan2(evec[1, 1], evec[0, 1]))))


def outline_at(ds, forbid, sig, x, y):
    """Trace, recentre, trace again. Returns (measures, xs, ys, wall_fraction) or None."""
    xs, ys, r, rim, probe = trace_rim(ds, x, y, forbid)
    m = shape_measures(xs, ys)
    if m is None:
        return None
    xs, ys, r, rim, probe = trace_rim(ds, m["cx"], m["cy"], forbid)
    m2 = shape_measures(xs, ys)
    if m2 is None:
        return None
    # contrast: how much darker the rim is than the brighter of (a) the core (inside half the rim
    # radius; entry cells have a bright core) and (b) the background just outside (r+2..r+4 px;
    # exit-funnel cells are a dark ring with no bright core), over the noise
    rr = np.median(r)
    core = ds[max(0, int(m2["cy"] - 0.5 * rr)):int(m2["cy"] + 0.5 * rr) + 1,
              max(0, int(m2["cx"] - 0.5 * rr)):int(m2["cx"] + 0.5 * rr) + 1]
    core_v = float(np.median(core)) if core.size else 0.0
    th = np.arange(N_ANG) * 2 * np.pi / N_ANG
    ring = []
    for dr in (2.0, 3.0, 4.0):
        ring.append(ndi.map_coordinates(ds, [m2["cy"] + (r + dr) * np.sin(th), m2["cx"] + (r + dr) * np.cos(th)],
                                        order=1, mode="nearest"))
    ok = ~(ndi.map_coordinates(forbid.astype(np.float32), [m2["cy"] + (r + 3) * np.sin(th), m2["cx"] + (r + 3) * np.cos(th)],
                               order=0, mode="constant", cval=1.0) > 0.5)
    outer_v = float(np.median(np.stack(ring)[:, ok])) if ok.any() else 0.0
    m2["contrast"] = float((max(core_v, outer_v) - np.median(rim)) / max(sig, 1e-3))
    return m2, xs, ys, float(probe.mean())


# ------------------------------------------------------------------------------------------
# background models and frame access
# ------------------------------------------------------------------------------------------

class Frames:
    def __init__(self, video):
        self.cap = cv2.VideoCapture(video)
        self.total = int(self.cap.get(cv2.CAP_PROP_FRAME_COUNT))
        self.pos = -1

    def get(self, f):
        if f != self.pos + 1:
            self.cap.set(cv2.CAP_PROP_POS_FRAMES, int(f))
        ok, im = self.cap.read()
        self.pos = f if ok else -1
        return cv2.cvtColor(im, cv2.COLOR_BGR2GRAY).astype(np.float32) if ok else None


class Side:
    """Geometry of one funnel crop (entry or exit)."""

    def __init__(self, name, lumen, top, bot, x0, x1):
        H, W = lumen.shape
        self.name, self.x0, self.x1 = name, x0, x1
        rows = np.where(lumen[:, x0:x1].any(axis=1))[0]
        self.y0, self.y1 = max(0, int(rows.min()) - 10), min(H, int(rows.max()) + 11)
        Y = np.arange(self.y0, self.y1)[:, None]
        X = np.arange(x0, x1)[None, :]
        t, b = top[x0:x1][None, :], bot[x0:x1][None, :]
        self.inside = (Y >= t - 1.0) & (Y <= b + 1.0)
        self.core = (Y >= t + 2.5) & (Y <= b - 2.5)      # away from the walls: noise estimate
        self.forbid = ~self.inside
        self.forbid[:, :1] |= (x0 == 0)
        self.forbid[:, -1:] |= (x1 == W)
        self.width = (bot - top)[x0:x1]
        self.win = cv2.createHanningWindow((x1 - x0, self.y1 - self.y0), cv2.CV_32F)

    def crop(self, g):
        return g[self.y0:self.y1, self.x0:self.x1]


def _shift(img, dx, dy):
    M = np.float32([[1, 0, dx], [0, 1, dy]])
    return cv2.warpAffine(img, M, (img.shape[1], img.shape[0]), flags=cv2.INTER_LINEAR,
                          borderMode=cv2.BORDER_REFLECT)


def align_to(ref, img, win):
    """Shift (dx, dy) of img relative to ref (phase correlation). The tube image wobbles by up
    to ~2 px between frames (different amounts at entry, bore and exit), which is what makes the
    walls flicker in a plain background difference."""
    (dx, dy), _ = cv2.phaseCorrelate(ref, img, win)
    if abs(dx) > 4 or abs(dy) > 4:
        return 0.0, 0.0
    return dx, dy


class Backgrounds:
    """Per block of BLOCK frames and per funnel crop: a background built from frames aligned to
    each other (so the walls stay sharp) and a per-pixel jitter floor of aligned differences."""

    def __init__(self, frames, sides):
        self.frames, self.sides, self.cache = frames, sides, {}

    def get(self, f):
        b = f // BLOCK
        if b not in self.cache:
            lo = max(0, b * BLOCK - BG_PAD)
            hi = min(self.frames.total, (b + 1) * BLOCK + BG_PAD)
            stacks = {s.name: [] for s in self.sides}
            for k in range(lo, hi, BG_STEP):
                g = self.frames.get(k)
                if g is None:
                    continue
                for s in self.sides:
                    stacks[s.name].append(s.crop(g))
            model = {}
            for s in self.sides:
                st = np.stack(stacks[s.name])
                ref = np.median(st, axis=0)
                for _ in range(2):            # align every sample to the current reference
                    al = np.stack([_shift(x, *(-np.array(align_to(ref, x, s.win)))) for x in st])
                    ref = np.median(al, axis=0)
                dd = []
                for x in al[::2]:
                    d = x - ref
                    d -= np.median(d)
                    dd.append(np.abs(cv2.GaussianBlur(d, (3, 3), 0)))
                model[s.name] = (ref, np.median(np.stack(dd), axis=0))
            self.cache = {b: model}            # keep one block in memory
        return self.cache[b]


def prepare(side, g, bgm, floor):
    c = side.crop(g)
    dx, dy = align_to(bgm, c, side.win)
    d = c - _shift(bgm, dx, dy)
    d -= np.median(d[side.inside])
    ds = cv2.GaussianBlur(d, (0, 0), 0.8)
    p16, p84 = np.percentile(ds[side.core], [16, 84])
    sig = 0.5 * float(p84 - p16)
    blob = cv2.GaussianBlur(np.maximum(0.0, np.abs(ds) - 1.5 * floor), (0, 0), 3.0)
    blob[~side.inside] = 0
    return d, ds, sig, blob


def candidates(blob, min_val):
    pk = (blob == ndi.maximum_filter(blob, size=9)) & (blob > min_val)
    ys, xs = np.nonzero(pk)
    return np.stack([xs, ys], 1).astype(float) if len(xs) else np.zeros((0, 2))


def zone_of(side, x_crop):
    xi = int(np.clip(round(x_crop), 0, len(side.width) - 1))
    return "wide" if side.width[xi] >= WIDE_MIN_WIDTH else "taper"


# ------------------------------------------------------------------------------------------
# entry funnel: outlines from the tracks
# ------------------------------------------------------------------------------------------

def funnel_limits(top, bot):
    w = bot - top
    W = len(w)
    ent = next(x for x in range(W // 2, W) if (w[x:] >= FUNNEL_MIN_WIDTH).all())
    ext = [x for x in range(0, W // 2) if (w[:x + 1] >= FUNNEL_MIN_WIDTH).all()][-1]
    return ent, ext


def flags_for(m, wall, cut, touching, zone):
    fl = []
    if cut:
        fl.append("cut_off")
    if wall > WALL_FLAG:
        fl.append("wall_contact")
    if touching:
        fl.append("touching_neighbour")
    if m["contrast"] < MIN_CONTRAST:
        fl.append("low_contrast")
    if m["solidity"] < SOLIDITY_FLAG:
        fl.append("irregular")
    if not fl and zone == "wide":
        q = "good"
    elif not (cut or touching or m["contrast"] < MIN_CONTRAST):
        q = "usable"
    else:
        q = "poor"
    return fl, q


def entry_pass(tracks, side, frames, bgs, entry_lo, W, log):
    """tracks: {tid: [(frame, cx, cy), ...]}. Returns list of per-frame result dicts."""
    by_frame = {}
    for tid, pts in tracks.items():
        for i, (f, cx, cy) in enumerate(pts):
            if cx >= entry_lo:
                vx = pts[min(i + 1, len(pts) - 1)][1] - pts[max(i - 1, 0)][1]
                by_frame.setdefault(f, []).append((tid, cx, cy, 0.5 * vx))
    last = {}
    out = []
    t0 = time.time()
    for n, f in enumerate(sorted(by_frame)):
        g = frames.get(f)
        if g is None:
            continue
        bgm, floor = bgs.get(f)[side.name]
        d, ds, sig, blob = prepare(side, g, bgm, floor)
        cand = candidates(blob, 0.5 * sig)
        act = by_frame[f]
        if len(cand) == 0:
            continue
        preds, gates = [], []
        for tid, cx, cy, vx in act:
            lf = last.get(tid)
            if lf is not None and f - lf[0] <= 3:
                p = (lf[1] + vx * (f - lf[0]), lf[2])
                gates.append(CONT_GATE + abs(vx))
            else:
                p = (cx - side.x0, cy - side.y0)
                gates.append(ENTRY_GATE)
            preds.append(p)
        preds = np.array(preds)
        D = np.hypot(preds[:, None, 0] - cand[None, :, 0], preds[:, None, 1] - cand[None, :, 1])
        Dg = np.where(D <= np.array(gates)[:, None], D, 1e6)
        ri, ci = linear_sum_assignment(Dg)
        found = []
        for a, c in zip(ri, ci):
            if Dg[a, c] >= 1e6:
                continue
            res = outline_at(ds, side.forbid, sig, cand[c, 0], cand[c, 1])
            if res is None:
                continue
            found.append((a, res))
        centres = [(res[0]["cx"], res[0]["cy"], np.sqrt(res[0]["area_px2"] / np.pi)) for _, res in found]
        peak = {a: blob[int(cand[c, 1]), int(cand[c, 0])] for a, c in zip(ri, ci) if Dg[a, c] < 1e6}
        for k, (a, (m, xs, ys, wall)) in enumerate(found):
            tid, tcx, tcy, vx = act[a]
            gx, gy = m["cx"] + side.x0, m["cy"] + side.y0
            cut = bool((xs + side.x0).max() > W - 3.0 or (xs + side.x0).min() < 2.0)
            rad = np.sqrt(m["area_px2"] / np.pi)
            touching = any(np.hypot(m["cx"] - c[0], m["cy"] - c[1]) < rad + c[2] + TOUCH_MARGIN
                           for j, c in enumerate(centres) if j != k)
            # a strong blob with no track yet (e.g. a cell entering at the frame edge) counts too
            for cx_, cy_ in cand:
                dd_ = np.hypot(cx_ - m["cx"], cy_ - m["cy"])
                if rad + 3 < dd_ < 2 * rad + TOUCH_MARGIN and blob[int(cy_), int(cx_)] >= 0.5 * peak[a]:
                    touching = True
            zone = zone_of(side, m["cx"])
            fl, q = flags_for(m, wall, cut, touching, zone)
            last[tid] = (f, m["cx"], m["cy"])
            row = dict(track_id=tid, frame=f, funnel="entry", zone=zone, track_cx=tcx, track_cy=tcy,
                       offset_px=float(np.hypot(gx - tcx, gy - tcy)), wall_contact=wall, cut_off=int(cut),
                       touching_neighbour=int(touching), quality=q, flags=";".join(fl), exit_tracklet="")
            row.update({k2: v for k2, v in m.items()})
            row["cx"], row["cy"] = gx, gy
            row["_xs"], row["_ys"] = xs + side.x0, ys + side.y0
            out.append(row)
        if log and n % 500 == 0:
            log(f"  entry: {n}/{len(by_frame)} frames, {len(out)} outlines, {time.time() - t0:.0f}s")
    return out


# ------------------------------------------------------------------------------------------
# exit funnel: blobs -> tracklets -> order-preserving match to tracks
# ------------------------------------------------------------------------------------------

def exit_pass(tracks, side, frames, bgs, exit_hi, log):
    """Detect blobs in the exit funnel in every frame around the tracks' exit times, link them
    into tracklets, and match tracklets to tracks in passage order."""
    # time each track reaches the bore exit (first frame with cx <= exit_hi + 10)
    t_exit = {}
    for tid, pts in tracks.items():
        hit = [f for f, cx, cy in pts if cx <= exit_hi + 10]
        if hit:
            t_exit[tid] = hit[0]
    if not t_exit:
        return [], {}
    spans = sorted((t - 50, t + 130) for t in t_exit.values())
    fr_set = set()
    for a, b in spans:
        fr_set.update(range(max(0, a), min(frames.total, b)))
    dets = {}
    t0 = time.time()
    for n, f in enumerate(sorted(fr_set)):
        g = frames.get(f)
        if g is None:
            continue
        bgm, floor = bgs.get(f)[side.name]
        d, ds, sig, blob = prepare(side, g, bgm, floor)
        cand = candidates(blob, EXIT_BLOB_MIN * sig)   # keeps ~86% of real exit cells, ~44% of noise
        keep = []
        for x, y in cand:
            if x > exit_hi + 2:      # still in the taper / bore
                continue
            res = outline_at(ds, side.forbid, sig, x, y)
            if (res is None or res[0]["contrast"] < EXIT_LINK_CONTRAST
                    or res[3] > 0.3 or res[0]["aspect_ratio"] > 2.5):
                continue                  # weak, hugging a wall, or a wall streak
            keep.append(res)
        # one outline per place (two candidates can converge on the same cell)
        uniq = []
        for r in sorted(keep, key=lambda r: -r[0]["contrast"]):
            if all(np.hypot(r[0]["cx"] - u[0]["cx"], r[0]["cy"] - u[0]["cy"]) > 6 for u in uniq):
                uniq.append(r)
        dets[f] = uniq
        if log and n % 1000 == 0:
            log(f"  exit: {n}/{len(fr_set)} frames, {time.time() - t0:.0f}s")
    # link: a detection continues a tracklet if within 5 px (x may only decrease by <= 6 px/frame
    # or increase by <= 2 px), allowing gaps of up to 3 frames
    tracklets = []
    active = []
    for f in sorted(dets):
        active = [t for t in active if f - t["frames"][-1] <= 4]
        used = set()
        if active and dets[f]:
            P = np.array([[t["pts"][-1][0], t["pts"][-1][1]] for t in active])
            Q = np.array([[r[0]["cx"], r[0]["cy"]] for r in dets[f]])
            dx = Q[None, :, 0] - P[:, None, 0]
            dy = Q[None, :, 1] - P[:, None, 1]
            ok = (dx >= -7) & (dx <= 2.5) & (np.abs(dy) <= 4)
            D = np.where(ok, np.hypot(dx, dy), 1e6)
            ri, ci = linear_sum_assignment(D)
            for a, c in zip(ri, ci):
                if D[a, c] < 1e6:
                    t = active[a]
                    t["frames"].append(f)
                    t["pts"].append((Q[c, 0], Q[c, 1]))
                    t["res"].append(dets[f][c])
                    used.add(c)
        for c, r in enumerate(dets[f]):
            if c not in used:
                t = dict(frames=[f], pts=[(r[0]["cx"], r[0]["cy"])], res=[r])
                tracklets.append(t)
                active.append(t)
    # keep tracklets that look like a cell leaving the bore: >= 8 detections, starting in the
    # narrow part of the funnel (x > 40 in crop coordinates = near the taper), moving left
    good = []
    for t in tracklets:
        if len(t["frames"]) < 8 or np.median([r[0]["contrast"] for r in t["res"]]) < EXIT_TRACKLET_CONTRAST:
            continue
        x_first, x_last = t["pts"][0][0], t["pts"][-1][0]
        if x_first < exit_hi - 45 or x_first - x_last < 10:
            continue
        good.append(t)
    good.sort(key=lambda t: t["frames"][0])
    for i, t in enumerate(good):
        t["id"] = i + 1
    # order-preserving match (both lists sorted by time); cost |dt|, skip cost 40 frames
    tids = sorted(t_exit, key=lambda k: t_exit[k])
    a = [t_exit[k] for k in tids]
    b = [t["frames"][0] for t in good]
    SK, MAXDT = 40.0, 60
    na, nb = len(a), len(b)
    Dp = np.full((na + 1, nb + 1), np.inf)
    Dp[0, :] = SK * np.arange(nb + 1)
    Dp[:, 0] = SK * np.arange(na + 1)
    mv = np.zeros((na + 1, nb + 1), np.int8)
    for i in range(1, na + 1):
        for j in range(1, nb + 1):
            dt = abs(a[i - 1] - b[j - 1])
            opts = [Dp[i - 1, j - 1] + (dt if dt <= MAXDT else np.inf), Dp[i - 1, j] + SK, Dp[i, j - 1] + SK]
            k = int(np.argmin(opts))
            Dp[i, j], mv[i, j] = opts[k], k
    match = {}
    i, j = na, nb
    while i > 0 and j > 0:
        if mv[i, j] == 0:
            match[tids[i - 1]] = good[j - 1]
            i, j = i - 1, j - 1
        elif mv[i, j] == 1:
            i -= 1
        else:
            j -= 1
    # label each match: clear if no other track reached the exit within 40 frames
    info = {}
    for tid in tids:
        near = [o for o in tids if o != tid and abs(t_exit[o] - t_exit[tid]) < 40]
        if tid in match:
            t = match[tid]
            info[tid] = dict(exit_match="matched_close_neighbour" if near else "matched",
                             exit_tracklet=t["id"], dt=t["frames"][0] - t_exit[tid])
        else:
            info[tid] = dict(exit_match="no_exit_outline", exit_tracklet="", dt=None)
    rows = []
    for tid, t in match.items():
        pts = dict((f, p) for f, p, _ in tracks_xy(tracks[tid]))
        for f, (m, xs, ys, wall), (px, py) in zip(t["frames"], t["res"], t["pts"]):
            gx, gy = m["cx"] + side.x0, m["cy"] + side.y0
            cut = bool((xs + side.x0).min() < 2.0)
            others = [r for r in dets[f] if abs(r[0]["cx"] - m["cx"]) + abs(r[0]["cy"] - m["cy"]) > 1e-6]
            touching = any(np.hypot(m["cx"] - r[0]["cx"], m["cy"] - r[0]["cy"])
                           < np.sqrt(m["area_px2"] / np.pi) + np.sqrt(r[0]["area_px2"] / np.pi) + TOUCH_MARGIN for r in others)
            zone = zone_of(side, m["cx"])
            fl, q = flags_for(m, wall, cut, touching, zone)
            tp = pts.get(f)
            row = dict(track_id=tid, frame=f, funnel="exit", zone=zone,
                       track_cx=tp[0] if tp else "", track_cy=tp[1] if tp else "",
                       offset_px=float(np.hypot(gx - tp[0], gy - tp[1])) if tp else "",
                       wall_contact=wall, cut_off=int(cut), touching_neighbour=int(touching),
                       quality=q, flags=";".join(fl), exit_tracklet=t["id"])
            row.update(m)
            row["cx"], row["cy"] = gx, gy
            row["_xs"], row["_ys"] = xs + side.x0, ys + side.y0
            rows.append(row)
    if log:
        log(f"  exit: {len(tracklets)} tracklets, {len(good)} cell-like, {len(match)}/{len(tids)} tracks matched")
    return rows, info


def tracks_xy(pts):
    for f, cx, cy in pts:
        yield f, (cx, cy), None


# ------------------------------------------------------------------------------------------
# per-cell summary
# ------------------------------------------------------------------------------------------

def _med(rows, key):
    v = [float(r[key]) for r in rows]
    return float(np.median(v)) if v else None


def summarise(tracks, bios, rows, exit_info):
    by = {}
    for r in rows:
        by.setdefault((r["track_id"], r["funnel"]), []).append(r)
    # size outliers relative to the cell's own median entry area
    for (tid, fun), rs in by.items():
        ref = _med([r for r in rs if r["quality"] in ("good", "usable")], "area_px2")
        if not ref:
            continue
        for r in rs:
            if abs(r["area_px2"] - ref) > SIZE_FLAG * ref and "size_outlier" not in r["flags"]:
                r["flags"] = ";".join([x for x in r["flags"].split(";") if x] + ["size_outlier"])
                if r["quality"] == "good":
                    r["quality"] = "usable"
    out = []
    for tid in sorted(tracks):
        b = bios.get(tid, {})
        ent = by.get((tid, "entry"), [])
        ext = by.get((tid, "exit"), [])
        rest = [r for r in ent if r["quality"] == "good"]
        taper = [r for r in ent if r["zone"] == "taper" and r["quality"] in ("good", "usable")]
        exg = [r for r in ext if r["quality"] == "good"]
        exu = exg or [r for r in ext if r["quality"] == "usable"]
        a = [r["area_px2"] for r in rest]
        ei = exit_info.get(tid, {})
        row = dict(track_id=tid, entrance_frame=b.get("entrance_frame", tracks[tid][0][0]),
                   entry_check=b.get("entry_check", ""),
                   n_entry=len(ent), n_entry_good=len(rest),
                   rest_area_px2=_med(rest, "area_px2"),
                   rest_area_iqr_px2=float(np.subtract(*np.percentile(a, [75, 25]))) if len(a) >= 2 else None,
                   rest_eq_diameter_px=_med(rest, "eq_diameter_px"), rest_perimeter_px=_med(rest, "perimeter_px"),
                   rest_circularity=_med(rest, "circularity"), rest_aspect_ratio=_med(rest, "aspect_ratio"),
                   rest_deformation_index=_med(rest, "deformation_index"), rest_contrast=_med(rest, "contrast"),
                   taper_n=len(taper),
                   taper_aspect_ratio_p90=float(np.percentile([r["aspect_ratio"] for r in taper], 90)) if len(taper) >= 3 else None,
                   taper_deformation_index_p90=float(np.percentile([r["deformation_index"] for r in taper], 90)) if len(taper) >= 3 else None,
                   exit_match=ei.get("exit_match", "no_exit_outline"), exit_tracklet=ei.get("exit_tracklet", ""),
                   exit_first_frame=min((r["frame"] for r in ext), default=None),
                   n_exit=len(ext), n_exit_good=len(exg),
                   exit_area_px2=_med(exu, "area_px2"), exit_eq_diameter_px=_med(exu, "eq_diameter_px"),
                   exit_circularity=_med(exu, "circularity"), exit_aspect_ratio=_med(exu, "aspect_ratio"),
                   exit_deformation_index=_med(exu, "deformation_index"), exit_contrast=_med(exu, "contrast"))
        row["outline_status"] = ("ok" if len(rest) >= 5 else "few_good_frames" if rest else "no_good_outline")
        size = row["rest_area_px2"] or _med([r for r in ent if r["quality"] != "poor"], "area_px2") or _med(ent, "area_px2")
        row["size_class"] = "unknown" if size is None else ("small_particle" if size < SMALL_OBJECT_AREA else "rbc_sized")
        out.append(row)
    return out


def load_frame_rows(path):
    """Read cell_outlines_per_frame.csv back (for --summary-only)."""
    rows = []
    ints = {"track_id", "frame", "cut_off", "touching_neighbour"}
    text = {"funnel", "zone", "quality", "flags", "exit_tracklet"}
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            for k, v in r.items():
                if k in ints:
                    r[k] = int(v)
                elif k not in text:
                    r[k] = _num(v)
            rows.append(r)
    return rows


# ------------------------------------------------------------------------------------------
# review sheets
# ------------------------------------------------------------------------------------------

QCOL = {"good": (0, 220, 0), "usable": (0, 220, 255), "poor": (0, 0, 255)}


def review_tile(frames, bgs, sides, row, R=28, S=3):
    g = frames.get(row["frame"])
    side = sides[row["funnel"]]
    bgm, floor = bgs.get(row["frame"])[side.name]
    d = prepare(side, g, bgm, floor)[0]
    cx, cy = row["cx"] - side.x0, row["cy"] - side.y0
    x0, y0 = int(round(cx)) - R, int(round(cy)) - R
    t = np.zeros((2 * R, 2 * R), np.float32)
    xs0, ys0 = max(0, x0), max(0, y0)
    xs1, ys1 = min(d.shape[1], x0 + 2 * R), min(d.shape[0], y0 + 2 * R)
    t[ys0 - y0:ys1 - y0, xs0 - x0:xs1 - x0] = d[ys0:ys1, xs0:xs1]
    im = cv2.cvtColor(cv2.resize(np.clip(128 + 4 * t, 0, 255).astype(np.uint8), None, fx=S, fy=S,
                                 interpolation=cv2.INTER_NEAREST), cv2.COLOR_GRAY2BGR)
    pts = np.stack([(row["_xs"] - side.x0 - x0 + 0.5) * S, (row["_ys"] - side.y0 - y0 + 0.5) * S], 1)
    cv2.polylines(im, [pts.astype(np.int32)], True, QCOL[row["quality"]], 1, cv2.LINE_AA)
    cv2.putText(im, f"f{row['frame']} {row['funnel'][:2]}", (2, 11), cv2.FONT_HERSHEY_SIMPLEX, 0.36, (0, 255, 255), 1, cv2.LINE_AA)
    cv2.putText(im, f"A{row['area_px2']:.0f} C{row['circularity']:.2f} AR{row['aspect_ratio']:.2f}",
                (2, 2 * R * S - 4), cv2.FONT_HERSHEY_SIMPLEX, 0.34, (0, 255, 255), 1, cv2.LINE_AA)
    return im


def review_sheets(out_dir, frames, bgs, sides, rows, cells, per_page=12, n_entry=6, n_exit=3):
    by = {}
    for r in rows:
        by.setdefault(r["track_id"], []).append(r)
    files = []
    lines = []
    for c in cells:
        rs = sorted(by.get(c["track_id"], []), key=lambda r: r["frame"])
        ent = [r for r in rs if r["funnel"] == "entry" and r["quality"] != "poor"] or [r for r in rs if r["funnel"] == "entry"]
        ext = [r for r in rs if r["funnel"] == "exit"]
        pick = [ent[i] for i in np.linspace(0, len(ent) - 1, min(n_entry, len(ent))).astype(int)] if ent else []
        pick += [ext[i] for i in np.linspace(0, len(ext) - 1, min(n_exit, len(ext))).astype(int)] if ext else []
        tiles = [review_tile(frames, bgs, sides, r) for r in pick]
        th, tw = 168, 168
        tiles += [np.zeros((th, tw, 3), np.uint8)] * (n_entry + n_exit - len(tiles))
        strip = np.hstack([np.hstack([t, np.full((th, 2, 3), 255, np.uint8)]) for t in tiles])
        bar = np.full((20, strip.shape[1], 3), 30, np.uint8)
        f = lambda v, fmt: ("-" if v is None else fmt.format(v))
        txt = (f"cell #{c['track_id']} entr {c['entrance_frame']} [{c['entry_check']}]  rest: A {f(c['rest_area_px2'], '{:.0f}')} "
               f"d {f(c['rest_eq_diameter_px'], '{:.1f}')} C {f(c['rest_circularity'], '{:.3f}')} AR {f(c['rest_aspect_ratio'], '{:.2f}')} "
               f"(n {c['n_entry_good']})  taper AR90 {f(c['taper_aspect_ratio_p90'], '{:.2f}')}  exit: {c['exit_match']} "
               f"A {f(c['exit_area_px2'], '{:.0f}')} C {f(c['exit_circularity'], '{:.2f}')} (n {c['n_exit']})")
        cv2.putText(bar, txt, (4, 14), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (0, 255, 255), 1, cv2.LINE_AA)
        lines.append(np.vstack([bar, strip, np.full((4, strip.shape[1], 3), 255, np.uint8)]))
    for p in range(0, len(lines), per_page):
        page = np.vstack(lines[p:p + per_page])
        fn = out_dir / f"cell_outline_review_{p // per_page + 1:02d}.png"
        cv2.imwrite(str(fn), page)
        files.append(fn)
    return files


# ------------------------------------------------------------------------------------------
# self-test on synthetic cells
# ------------------------------------------------------------------------------------------

def self_test():
    prof_r = np.arange(0.5, 18, 1.0)
    prof_v = np.array([37, 38, 42, 39, 26, 7, -11, -19, -18, -12, -5, -2, 1, 1, -1, 1, 0, 0], float)
    rng = np.random.default_rng(1)
    sz, up = 60, 4
    print("true AR  contrast | area (rim)  | AR meas     | circularity true / meas")
    for scale in (1.0, 0.25):
        for ar in (1.0, 1.2, 1.5, 2.0):
            res = []
            a, b = 7.5 * np.sqrt(ar), 7.5 / np.sqrt(ar)
            for ang in (0, 30, 60):
                for _ in range(3):
                    c = sz / 2 + rng.uniform(-0.5, 0.5, 2)
                    t = np.radians(ang)
                    yy, xx = np.mgrid[0:sz * up, 0:sz * up].astype(float) / up
                    u = (xx - c[0]) * np.cos(t) + (yy - c[1]) * np.sin(t)
                    v = -(xx - c[0]) * np.sin(t) + (yy - c[1]) * np.cos(t)
                    ins = (u / a) ** 2 + (v / b) ** 2 <= 1
                    sd = (ndi.distance_transform_edt(ins) - ndi.distance_transform_edt(~ins)) / up
                    sd = cv2.resize(sd.astype(np.float32), (sz, sz), interpolation=cv2.INTER_AREA)
                    img = (np.interp(7.5 - sd, prof_r, prof_v) * scale + rng.normal(0, 3.5, (sz, sz))).astype(np.float32)
                    ds = cv2.GaussianBlur(img, (0, 0), 0.8)
                    g0 = c + rng.uniform(-2, 2, 2)
                    r = outline_at(ds, np.zeros((sz, sz), bool), 3.0, g0[0], g0[1])
                    res.append((r[0]["area_px2"], r[0]["aspect_ratio"], r[0]["circularity"]))
            res = np.array(res)
            h = ((a - b) / (a + b)) ** 2
            P = np.pi * (a + b) * (1 + 3 * h / (10 + np.sqrt(4 - 3 * h)))
            print(f"  {ar:4.1f}    {scale:4.2f}   | {res[:, 0].mean():6.1f}+-{res[:, 0].std():4.1f} | "
                  f"{res[:, 1].mean():.2f}+-{res[:, 1].std():.2f} | {4 * np.pi * (np.pi * a * b) / P ** 2:.3f} /{res[:, 2].mean():.3f}+-{res[:, 2].std():.3f}")


# ------------------------------------------------------------------------------------------
# main
# ------------------------------------------------------------------------------------------

def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def write_csv(path, fields, rows):
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else ("" if v is None else v)) for k, v in r.items() if k in fields})
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description="Trace RBC outlines in the entry and exit funnels.")
    ap.add_argument("--results", type=str, default=None,
                    help="Folder with trajectories.csv (+ biomechanics.csv) from run_full_video.py "
                         "(default: Desktop/RBC_full_video_results)")
    ap.add_argument("--video", type=str, default=None)
    ap.add_argument("--from-frame", type=int, default=None, help="Only cells entering at or after this frame")
    ap.add_argument("--to-frame", type=int, default=None, help="Only cells entering before this frame")
    ap.add_argument("--out", type=str, default=None, help="Write outputs here instead of --results")
    ap.add_argument("--no-exit", action="store_true", help="Entry funnel only")
    ap.add_argument("--no-review", action="store_true", help="Skip the review PNG sheets")
    ap.add_argument("--self-test", action="store_true", help="Run the synthetic-cell check and exit")
    ap.add_argument("--summary-only", action="store_true",
                    help="Rebuild cell_outlines_per_cell.csv from an existing cell_outlines_per_frame.csv (no video)")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    import run_full_video as RF
    res_dir = Path(args.results) if args.results else RF.default_out_dir()
    out_dir = Path(args.out) if args.out else res_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    traj = res_dir / "trajectories.csv"
    if not traj.is_file():
        sys.exit(f"No trajectories.csv in {res_dir}")

    tracks = {}
    with open(traj, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            tracks.setdefault(int(r["track_id"]), []).append((int(r["frame"]), float(r["cx"]), float(r["cy"])))
    for t in tracks.values():
        t.sort()
    bios = {}
    if (res_dir / "biomechanics.csv").is_file():
        with open(res_dir / "biomechanics.csv", newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                bios[int(r["track_id"])] = {"entrance_frame": int(float(r["entrance_frame"])),
                                            "entry_check": r.get("entry_check", "")}
    ent_of = {tid: bios.get(tid, {}).get("entrance_frame", p[0][0]) for tid, p in tracks.items()}
    sel = {tid: p for tid, p in tracks.items()
           if (args.from_frame is None or ent_of[tid] >= args.from_frame)
           and (args.to_frame is None or ent_of[tid] < args.to_frame)}
    if not sel:
        sys.exit("No cells in the selected range.")

    if args.summary_only:
        pf = out_dir / "cell_outlines_per_frame.csv"
        if not pf.is_file():
            sys.exit(f"No {pf.name} in {out_dir}; run without --summary-only first.")
        rows = [r for r in load_frame_rows(pf) if r["track_id"] in sel]
        exit_info = {}
        for r in rows:
            if r["funnel"] == "exit" and r["exit_tracklet"]:
                exit_info[r["track_id"]] = dict(exit_match="matched", exit_tracklet=r["exit_tracklet"])
        old = out_dir / "cell_outlines_per_cell.csv"
        if old.is_file():                      # keep the exit match labels from the full run
            with open(old, newline="", encoding="utf-8") as f:
                for r in csv.DictReader(f):
                    if int(r["track_id"]) in exit_info:
                        exit_info[int(r["track_id"])]["exit_match"] = r["exit_match"]
        cells = summarise(sel, bios, rows, exit_info)
        write_csv(out_dir / "cell_outlines_per_frame.csv", FRAME_FIELDS, rows)
        write_csv(out_dir / "cell_outlines_per_cell.csv", CELL_FIELDS, cells)
        print(f"per-cell summary rebuilt for {len(cells)} cells in {out_dir}")
        return

    calib = T.load_calibration()
    bg = T.load_background()
    video = args.video
    if not video:
        sp = res_dir / "settings.json"
        video = json.loads(sp.read_text(encoding="utf-8")).get("video") if sp.is_file() else None
        if not video or not os.path.isfile(video):
            video = T._resolve_video_path(calib.get("video"))
    if not video or not os.path.isfile(video):
        sys.exit(f"Video not found: {video!r}. Pass --video.")
    lumen, top, bot, _, _ = T.build_exact_lumen_geometry(bg)
    H, W = bg.shape[:2]
    entry_lo, exit_hi = funnel_limits(top, bot)
    sides = {"entry": Side("entry", lumen, top, bot, ENTRY_X0, W), "exit": Side("exit", lumen, top, bot, 0, EXIT_X1)}
    frames = Frames(video)
    bgs = Backgrounds(frames, list(sides.values()))
    log = lambda s: print(s, flush=True)
    log(f"cells: {len(sel)}  entry funnel x >= {entry_lo} (wide x >= "
        f"{next(x for x in range(W // 2, W) if ((bot - top)[x:] >= WIDE_MIN_WIDTH).all())}), exit funnel x <= {exit_hi}")
    t0 = time.time()
    rows = entry_pass(sel, sides["entry"], frames, bgs, entry_lo, W, log)
    log(f"entry funnel: {len(rows)} outlines in {time.time() - t0:.0f}s")
    exit_info = {}
    if not args.no_exit:
        t1 = time.time()
        erows, exit_info = exit_pass(sel, sides["exit"], frames, bgs, exit_hi - sides["exit"].x0, log)
        rows += erows
        log(f"exit funnel: {len(erows)} outlines in {time.time() - t1:.0f}s")
    cells = summarise(sel, bios, rows, exit_info)
    rows.sort(key=lambda r: (r["track_id"], r["funnel"] != "entry", r["frame"]))
    write_csv(out_dir / "cell_outlines_per_frame.csv", FRAME_FIELDS, rows)
    write_csv(out_dir / "cell_outlines_per_cell.csv", CELL_FIELDS, cells)
    np.savez_compressed(out_dir / "cell_outline_contours.npz",
                        track_id=np.array([r["track_id"] for r in rows], np.int32),
                        frame=np.array([r["frame"] for r in rows], np.int32),
                        funnel=np.array([r["funnel"] for r in rows]),
                        xs=np.array([r["_xs"] for r in rows], np.float32).reshape(-1, N_ANG),
                        ys=np.array([r["_ys"] for r in rows], np.float32).reshape(-1, N_ANG))
    if not args.no_review:
        files = review_sheets(out_dir, frames, bgs, sides, rows, cells)
        log(f"review sheets: {', '.join(p.name for p in files)}")
    n_ok = sum(c["outline_status"] == "ok" for c in cells)
    n_ex = sum(c["exit_match"].startswith("matched") for c in cells)
    log(f"done in {time.time() - t0:.0f}s: {len(cells)} cells, {n_ok} with >= 5 good entry outlines, "
        f"{n_ex} with an exit outline match. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
