"""
ZONE 1 / ZONE 3 TRACKING — entry-funnel / exit-funnel detection + tracking.

New architecture per the theory spec worked out in chat (2026-11-15):
  - Detection: background-subtraction + motion fusion, bilateral smoothing,
    Otsu threshold, size filter, watershed split for touching cells.
  - Tracking: per-track Kalman filter (constant-velocity) predicts next
    position; Hungarian algorithm (scipy linear_sum_assignment) does the
    frame-to-frame assignment, gated to a max-plausible-displacement radius.

Deliberately separate from tracking.py / preprocessing.py (do not touch
those - this is a fresh module, reusing Stage-1's calibration.json and
median_bg.npy as read-only inputs).

Zone definitions (pixel coords, this video's calibration):
    X_IN_END     = 110   -> entry funnel / capillary boundary
    X_OUT_START  = 1180  -> capillary / exit funnel boundary
  Flow direction is right -> left (matches tracking.py's own convention):
    Zone 1 "entry funnel"  = x > X_OUT_START  (cells arrive from the bath)
    Zone 3 "exit funnel"   = x < X_IN_END     (cells leave towards the bath)
  The capillary itself (X_IN_END..X_OUT_START) is NOT handled here - that
  is the kymograph + ADI_proxy zone, a separate module.

This file only implements Zone 1 / Zone 3. Run as a script to test on a
handful of real seeked frames (no full-video pass).
"""
import json
import cv2
import numpy as np
from scipy.optimize import linear_sum_assignment

# ====================================================================
# PATHS (read-only reuse of existing Stage-1 output)
# ====================================================================
CALIB_FILE = r"C:\Users\amara\rbc_tracking\output\calibration.json"
BG_FILE = r"C:\Users\amara\rbc_tracking\output\median_bg.npy"
VIDEO = r"C:\Users\amara\Downloads\2A002 2_C001H001S0001-001.avi"

# ====================================================================
# ZONE BOUNDARIES (pixels; this video's calibration)
# ====================================================================
X_IN_END = 110
X_OUT_START = 1180

# Measured from this footage's free-floating (undeformed) RBCs (vision-
# assisted bounding-box check, 2026-11-15): diameter ~18 px.
RESTING_DIAM_PX = 18.0
_RESTING_AREA = np.pi * (RESTING_DIAM_PX / 2.0) ** 2
MIN_AREA = 0.25 * _RESTING_AREA          # reject dust/noise speckles
MAX_AREA = 2.5 * _RESTING_AREA           # above this: probably 2+ touching cells

GATE_RADIUS_PX = 20.0   # max plausible per-frame displacement for Hungarian gating
MAX_MISSES = 5          # frames a track can go unmatched before death

# Measured on real frames (1500-1649): fused-signal max() stays ~35-55 when
# no real cell is present in a flare ROI (pure background/wall-jitter noise),
# vs ~70-470 when a real cell is in view. These floors reject noise-only
# frames instead of letting Otsu invent a false split on flat noise.
MIN_SIGNAL_FLOOR = 60.0
MIN_OTSU_FLOOR = 25.0


def load_calibration():
    with open(CALIB_FILE, "r", encoding="utf-8") as f:
        return json.load(f)


def load_background():
    return np.load(BG_FILE)


# ====================================================================
# ROI: generous band around the centerline, FLARE REGIONS ONLY
# (the wide bath is full of unrelated free-floating cells - restricting
# to the flares is what keeps detection scoped to Zone 1 / Zone 3)
# ====================================================================
def build_flare_roi_mask(shape, slope, intercept, margin=45):
    h, w = shape
    mask = np.zeros((h, w), dtype=np.uint8)
    for x in range(w):
        if x <= X_IN_END or x >= X_OUT_START:
            yc = slope * x + intercept
            y0 = max(0, int(yc - margin))
            y1 = min(h, int(yc + margin))
            mask[y0:y1, x] = 255
    return mask


# ====================================================================
# DETECTION
# ====================================================================
def _watershed_split(binary, labels, label_id, stats):
    x, y, w, h, area = stats[label_id]
    pad = 5
    x0, y0 = max(0, x - pad), max(0, y - pad)
    x1, y1 = min(binary.shape[1], x + w + pad), min(binary.shape[0], y + h + pad)
    blob_mask = (labels[y0:y1, x0:x1] == label_id).astype(np.uint8) * 255

    dist = cv2.distanceTransform(blob_mask, cv2.DIST_L2, 5)
    if dist.max() <= 0:
        ys, xs = np.where(blob_mask > 0)
        return [{"cx": float(x0 + xs.mean()), "cy": float(y0 + ys.mean()), "area": float(area)}]

    _, sure_fg = cv2.threshold(dist, 0.5 * dist.max(), 255, 0)
    sure_fg = sure_fg.astype(np.uint8)
    n2, markers = cv2.connectedComponents(sure_fg)
    if n2 <= 2:
        ys, xs = np.where(blob_mask > 0)
        return [{"cx": float(x0 + xs.mean()), "cy": float(y0 + ys.mean()), "area": float(area)}]

    markers = markers + 1
    markers[blob_mask == 0] = 0
    blob_bgr = cv2.cvtColor(blob_mask, cv2.COLOR_GRAY2BGR)
    cv2.watershed(blob_bgr, markers)

    out = []
    for lbl in range(2, n2 + 1):
        ys, xs = np.where(markers == lbl)
        if xs.size < MIN_AREA * 0.3:
            continue
        out.append({"cx": float(x0 + xs.mean()), "cy": float(y0 + ys.mean()), "area": float(xs.size)})
    return out if out else [{"cx": float(x + w / 2.0), "cy": float(y + h / 2.0), "area": float(area)}]


def detect_blobs(gray, prev_gray, bg, roi_mask):
    """Fused background-subtraction (signal A) + motion (signal B) detection.

    Both signals are needed: static background-subtracted contrast alone was
    found (earlier optical-flow check) too subtle for reliable detection on
    this footage; motion evidence recovers cells that are faint in contrast.
    """
    gray_f = gray.astype(np.float32)
    bg_f = bg.astype(np.float32)

    # Darkening-only contrast: RBCs absorb light and appear DARKER than
    # background in this brightfield footage (confirmed in Stage-1 calibration:
    # cell pixels ~0-40 vs background ~100-115). A plain abs(gray-bg) also
    # fires on anything that gets BRIGHTER than background - found in testing
    # to include a transient bright artifact at the extreme frame corner
    # (debris/reflection/vignetting, not a cell) that dragged the per-frame
    # Otsu threshold up ~10x while active and caused severe track flicker.
    # Restricting to darkening removes that whole class of false positives.
    dark_diff = np.maximum(0.0, bg_f - gray_f)
    motion = np.abs(gray_f - prev_gray.astype(np.float32)) if prev_gray is not None else np.zeros_like(dark_diff)

    fused = dark_diff * np.sqrt(np.maximum(1.0, motion))
    fused[roi_mask == 0] = 0

    fused_u8 = np.clip(fused, 0, 255).astype(np.uint8)
    # bilateral: edge-preserving, does not blur cell boundaries the way a
    # Gaussian/box blur would (that was Kumar 2020's flaw - do not repeat it)
    fused_smooth = cv2.bilateralFilter(fused_u8, d=7, sigmaColor=50, sigmaSpace=4)

    roi_vals = fused_smooth[roi_mask > 0]
    # MIN_SIGNAL_FLOOR: measured on real frames (1500-1649) that when no real
    # cell is present in a flare, fused.max() stays in the ~35-55 noise band,
    # while frames with a real cell reach ~70-470. Below this floor, Otsu
    # still finds *some* split in pure noise and reports fake blobs (observed
    # directly: frame with fused.max()=34.9 produced 20 "detections" on a
    # frame with zero visible cells) - worse than reporting nothing.
    if roi_vals.size == 0 or roi_vals.max() < MIN_SIGNAL_FLOOR:
        return [], fused_smooth

    # Otsu: threshold computed from this frame's own histogram, not a fixed
    # hand-picked constant - keeps the pipeline deterministic/self-calibrating.
    # Clamped to a minimum so a mostly-flat histogram (low true contrast) does
    # not collapse to a near-zero split.
    thresh_val, _ = cv2.threshold(roi_vals, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    thresh_val = max(thresh_val, MIN_OTSU_FLOOR)
    binary = (fused_smooth > thresh_val).astype(np.uint8) * 255
    binary[roi_mask == 0] = 0

    n, labels, stats, centroids = cv2.connectedComponentsWithStats(binary, connectivity=8)

    detections = []
    for i in range(1, n):
        area = stats[i, cv2.CC_STAT_AREA]
        if area < MIN_AREA:
            continue
        if area > MAX_AREA:
            detections.extend(_watershed_split(binary, labels, i, stats))
        else:
            cx, cy = centroids[i]
            detections.append({"cx": float(cx), "cy": float(cy), "area": float(area)})
    return detections, fused_smooth


# ====================================================================
# TRACKING: per-track Kalman filter + Hungarian assignment
# ====================================================================
class Track:
    _next_id = 1

    def __init__(self, cx, cy):
        self.id = Track._next_id
        Track._next_id += 1
        self.kf = cv2.KalmanFilter(4, 2)
        self.kf.transitionMatrix = np.array(
            [[1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0], [0, 0, 0, 1]], np.float32
        )
        self.kf.measurementMatrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0]], np.float32)
        self.kf.processNoiseCov = np.eye(4, dtype=np.float32) * 0.3
        self.kf.measurementNoiseCov = np.eye(2, dtype=np.float32) * 4.0
        self.kf.statePost = np.array([[cx], [cy], [0], [0]], np.float32)
        self.kf.errorCovPost = np.eye(4, dtype=np.float32) * 10.0
        self.history = [(cx, cy)]
        self.misses = 0
        self.hits = 1

    def predict(self):
        pred = self.kf.predict()
        return float(pred[0, 0]), float(pred[1, 0])

    def update(self, cx, cy):
        self.kf.correct(np.array([[cx], [cy]], np.float32))
        self.history.append((cx, cy))
        self.misses = 0
        self.hits += 1


def associate(tracks, detections):
    if not tracks or not detections:
        return [], list(range(len(detections))), list(range(len(tracks)))

    preds = [t.predict() for t in tracks]
    cost = np.full((len(tracks), len(detections)), 1e6)
    for i, (px, py) in enumerate(preds):
        for j, d in enumerate(detections):
            dist = float(np.hypot(px - d["cx"], py - d["cy"]))
            if dist <= GATE_RADIUS_PX:
                cost[i, j] = dist

    row_ind, col_ind = linear_sum_assignment(cost)
    matches, assigned_t, assigned_d = [], set(), set()
    for r, c in zip(row_ind, col_ind):
        if cost[r, c] < 1e6:
            matches.append((r, c))
            assigned_t.add(r)
            assigned_d.add(c)

    unmatched_tracks = [i for i in range(len(tracks)) if i not in assigned_t]
    unmatched_dets = [j for j in range(len(detections)) if j not in assigned_d]
    return matches, unmatched_dets, unmatched_tracks


def run_tracking(video_path, start_frame, n_frames, bg, roi_mask, verbose=True):
    cap = cv2.VideoCapture(video_path)
    cap.set(cv2.CAP_PROP_POS_FRAMES, start_frame)

    tracks = []
    dead_tracks = []
    prev_gray = None
    log = []

    for offset in range(n_frames):
        ok, frame = cap.read()
        if not ok:
            break
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        dets, _ = detect_blobs(gray, prev_gray, bg, roi_mask)
        matches, unmatched_dets, unmatched_tracks = associate(tracks, dets)

        for r, c in matches:
            tracks[r].update(dets[c]["cx"], dets[c]["cy"])
        for i in unmatched_tracks:
            tracks[i].misses += 1
        for j in unmatched_dets:
            tracks.append(Track(dets[j]["cx"], dets[j]["cy"]))

        alive, died = [], []
        for t in tracks:
            (died if t.misses > MAX_MISSES else alive).append(t)
        tracks = alive
        dead_tracks.extend(died)

        log.append({
            "frame": start_frame + offset,
            "n_detections": len(dets),
            "n_tracks_alive": len(tracks),
        })
        if verbose and (offset + 1) % 20 == 0:
            print(f"  frame {start_frame + offset}: {len(dets)} detections, {len(tracks)} live tracks")

        prev_gray = gray

    cap.release()
    return tracks + dead_tracks, log


# ====================================================================
# TEST (no full-video pass - a handful of seeked/consecutive frames only)
# ====================================================================
if __name__ == "__main__":
    calib = load_calibration()
    bg = load_background()
    slope = calib["channel_line"]["slope"]
    intercept = calib["channel_line"]["intercept"]

    roi_mask = build_flare_roi_mask(bg.shape, slope, intercept, margin=45)
    print(f"ROI mask pixel count: {int((roi_mask > 0).sum())} "
          f"({100.0 * (roi_mask > 0).sum() / roi_mask.size:.2f}% of frame)")

    # Test range: a known-active window from tracking.py's own debug output
    # (Cell #1 observed crossing Frames 1538 -> 1894 in the existing Stage 3
    # kymograph tracker, i.e. the entrance flare is active around here).
    START = 1500
    N_FRAMES = 150

    print(f"\nRunning Zone 1/3 detect+track on frames {START}..{START + N_FRAMES - 1}")
    all_tracks, log = run_tracking(VIDEO, START, N_FRAMES, bg, roi_mask)

    print(f"\nTotal tracks created: {len(all_tracks)}")
    for t in sorted(all_tracks, key=lambda t: -t.hits)[:10]:
        xs = [p[0] for p in t.history]
        ys = [p[1] for p in t.history]
        print(f"  track #{t.id}: {t.hits} hits, x {xs[0]:.1f}->{xs[-1]:.1f}, "
              f"y {ys[0]:.1f}->{ys[-1]:.1f}")

    # quick overlay dump on the last frame processed, for visual sanity check
    cap = cv2.VideoCapture(VIDEO)
    cap.set(cv2.CAP_PROP_POS_FRAMES, START + N_FRAMES - 1)
    ok, last_frame = cap.read()
    cap.release()
    if ok:
        vis = last_frame.copy()
        for t in all_tracks:
            if len(t.history) >= 2:
                pts = np.array(t.history, dtype=np.int32)
                cv2.polylines(vis, [pts], False, (0, 255, 0), 1, cv2.LINE_AA)
            cx, cy = t.history[-1]
            cv2.circle(vis, (int(cx), int(cy)), 4, (0, 0, 255), 1)
        out_path = r"C:\Users\amara\iisc-rbc-tracking\lab\evidence\zone13_test_overlay.png"
        cv2.imwrite(out_path, vis)
        print(f"\nOverlay saved: {out_path}")
