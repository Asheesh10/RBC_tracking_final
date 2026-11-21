"""
Whole-video RBC passage tracking.

Runs the Stage 3 pipeline from tracking.py (kymograph DP + gap pass + neighbour-masked
re-trace + RTS smoother + exit-flare extension) over the full AVI in overlapping windows and
stitches the results. tracking.py itself is not modified.

Why windows: tracking.py keeps every decoded ROI strip in RAM (~326 KB per frame), so the
whole 37,274-frame video would need ~12 GB. Each window here holds ~2,400 frames (~0.8 GB).

Window layout (defaults):
    core      = 1500 frames  the window size validated on 1450, 15000 and 18000
    pre-roll  = 450 frames   decoded before the core, so a cell entering just before the core
                             is detected and masked during the re-trace
    post-roll = 450 frames   decoded after the core; longer than the longest passage (435
                             frames), so every cell that ENTERS in the core finishes in-window
A cell belongs to the window whose core contains its entrance frame, so no cell is counted
twice and none falls between windows. Cores start at 0, 1500, 3000, ... so the core
15000-16500 is exactly the window checked by eye on 2026-11-15.

Not trackable by design: cells already inside the channel at frame 0, and cells entering in
roughly the last 260 frames of the recording (no room to finish a passage).

Offset grid (default on; --no-offset-grid to skip), added 2026-11-18:
    Whether tracking.py finds a FAINT cell depends on where its window starts (thresholds and
    gap intervals are computed from the decoded frames; measured on 15000-16500: the faint
    15459 cell is found when the window starts at 15000 but not when it starts at 14550). So a
    second grid B is run with cores shifted by core/2 (750 frames). Grid A (cores 0, 1500, ...)
    is the base and is never changed. A grid-B track is
      - a duplicate of an A track if, over >= 20 shared frames in the bore, the median distance
        between the two paths is < 20 px (the A track is then marked "A+B": found twice);
      - a conflict if it shares part of a ridge with an A track (>30 % of shared bore frames
        within 14 px) or crosses it (cells cannot overtake in the bore) -- reported, not added;
      - otherwise a candidate new cell. It is added (source "B") only if its ridge evidence
        (energy along the track vs random lines) is >= RIDGE_MIN_RATIO (2.5); else "rejected".
        Without this check the first whole-video merge added 14 tracks of which 11 were noise
        (empty entry funnel, evidence 1.1-1.9); the 3 real ones scored 2.75-3.34.
    Final cells are renumbered by entrance frame. merge_report.csv lists every grid-B track and
    its verdict so additions can be checked on the kymographs.

Outputs (default folder: <Desktop>/RBC_full_video_results):
    trajectories.csv   every cell, every frame (same columns as tracking.py)
    biomechanics.csv   one row per cell (+ window_core_start, source = A / A+B / B,
                       entry_cell_score, entry_check = cell_seen / uncertain / no_cell_seen)
    entry_check_review.png  entry-funnel crops of every track not flagged cell_seen, for review
                       (the entry check only FLAGS tracks; nothing is removed -- see EntryChecker)
    windows.csv        per window and grid: frame ranges, cells detected / kept, runtime
    merge_report.csv   every grid-B track: added / duplicate / conflict, and which A track
    summary.json       totals + settings
    videos/            one tracked H.264 video per grid-A core (every cell present is drawn,
                       including cells from the previous window still passing)
    kymographs/        one kymograph PNG per grid-A core with the tracks drawn on it -- every
                       bright diagonal ridge should have a coloured line on it
    per_window/        grid A raw results per window (written atomically; used for --resume)
    per_window_b/      grid B raw results per window
    run_log.txt        console output

Usage:
    py run_full_video.py --resume                       # whole video; continues if stopped
    py run_full_video.py --start 15000 --end 16500 --out <folder> --overwrite   # one window
"""
import argparse
import csv
import json
import os
import shutil
import sys
import time
from pathlib import Path

import cv2
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tracking as T  # noqa: E402

TRAJ_FIELDS = ["frame", "track_id", "cx", "cy", "bbox_x", "bbox_y", "bbox_w", "bbox_h",
               "velocity_px_per_frame", "area", "circularity", "elongation", "ADI_running"]
BIO_FIELDS = ["track_id", "entrance_frame", "exit_frame", "track_length_frames", "duration_ms",
              "x_start_px", "x_end_px", "net_travel_px", "velocity_mean", "velocity_max",
              "area_mean", "area_min", "circularity_mean", "circularity_min", "elongation_max",
              "ADI", "shape_recovery", "window_core_start"]
FINAL_BIO_FIELDS = BIO_FIELDS + ["source", "entry_cell_score", "entry_check"]
# Entry-funnel cell check thresholds. Calibrated 2026-11-18 against tracks labelled by eye from
# entry-funnel crops: 21 empty-funnel (noise) tracks scored 0.32-1.28 (median 0.59), 30 real
# cells 1.61-50 (median 5.4). >= ENTRY_SEEN: cell_seen; < ENTRY_NONE: no_cell_seen; between:
# uncertain (look at entry_check_review.png).
ENTRY_SEEN, ENTRY_NONE = 2.0, 1.0
WIN_FIELDS = ["core_start", "core_end", "decode_start", "decode_end", "cells_detected_in_window",
              "cells_kept", "kept_track_ids", "seconds"]
WIN_OUT_FIELDS = ["grid"] + WIN_FIELDS
REPORT_FIELDS = ["grid_b_track", "entrance_frame", "exit_frame", "verdict", "final_track_id",
                 "matched_grid_a_track", "shared_bore_frames", "median_distance_px", "ridge_ratio", "note"]
OUR_FILES = ["trajectories.csv", "biomechanics.csv", "windows.csv", "summary.json",
             "settings.json", "run_log.txt", "merge_report.csv", "render_state.json", "ridge_scores.json",
             "entry_scores.json", "entry_check_review.png"]
# A grid-B track that matches no grid-A track is added only if its ridge evidence
# (ridge_evidence) reaches this. Calibrated 2026-11-18 on the first whole-video merge: the 14
# grid-B-only tracks split into 3 real cells (round cell visible in the entry funnel; ratios
# 2.75-3.34, like confirmed faint cells at 3.0-3.1) and 11 noise tracks (empty funnel; 1.11-1.92).
RIDGE_MIN_RATIO = 2.5
OUR_DIRS = ["videos", "kymographs", "per_window", "per_window_b"]
B_ID_BASE = 10_000_000     # grid-B raw track ids: B_ID_BASE + core_start * 100 + k (never clash with A)
SETTING_KEYS = ("video", "start", "end", "core", "pre", "post")


def default_out_dir():
    home = Path.home()
    for desk in (home / "OneDrive" / "Desktop", home / "Desktop"):
        if desk.is_dir():
            return desk / "RBC_full_video_results"
    return home / "RBC_full_video_results"


class Tee:
    """Mirror console output into run_log.txt."""

    def __init__(self, path, stream):
        self.f = open(path, "a", encoding="utf-8")
        self.out = stream

    def write(self, s):
        self.out.write(s)
        self.f.write(s)

    def flush(self):
        self.out.flush()
        self.f.flush()


def plan_windows(total, start, end, core, pre, post):
    wins = []
    cs = start
    while cs < end:
        ce = min(cs + core, end)
        wins.append({"core_start": cs, "core_end": ce,
                     "decode_start": max(0, cs - pre), "decode_end": min(total, ce + post)})
        cs = ce
    return wins


# --------------------------------------------------------------------------------------
# persistence: one atomic set of files per window; merged CSVs are rebuilt from them
# --------------------------------------------------------------------------------------

def write_csv_atomic(path, fields, rows):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def read_csv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def traj_rows(gid, path):
    rows = []
    for p in path:
        bx, by, bw, bh = p["bbox"]
        rows.append({
            "frame": p["frame"], "track_id": gid,
            "cx": round(p["cx"], 2), "cy": round(p["cy"], 2),
            "bbox_x": bx, "bbox_y": by, "bbox_w": bw, "bbox_h": bh,
            "velocity_px_per_frame": round(p["velocity"], 3), "area": round(p["area"], 2),
            "circularity": round(p["circularity"], 4), "elongation": round(p["elongation"], 3),
            "ADI_running": round(p["ADI_running"], 4),
        })
    return rows


def rows_to_tracks(rows):
    tracks = {}
    for r in rows:
        gid = int(r["track_id"])
        tracks.setdefault(gid, []).append({
            "frame": int(r["frame"]), "track_id": gid, "cx": float(r["cx"]), "cy": float(r["cy"]),
            "bbox": (int(r["bbox_x"]), int(r["bbox_y"]), int(r["bbox_w"]), int(r["bbox_h"])),
            "velocity": float(r["velocity_px_per_frame"]), "area": float(r["area"]),
            "circularity": float(r["circularity"]), "elongation": float(r["elongation"]),
            "ADI_running": float(r["ADI_running"]),
        })
    for path in tracks.values():
        path.sort(key=lambda p: p["frame"])
    return tracks


def load_done(out, wins, sub="per_window"):
    """Return {core_start: info} for finished windows plus all their tracks / biomarkers."""
    done, tracks, bios = {}, {}, {}
    pw = out / sub
    for win in wins:
        cs = win["core_start"]
        info_p = pw / f"info_{cs:05d}.json"
        if not info_p.exists():
            continue
        done[cs] = json.loads(info_p.read_text(encoding="utf-8"))
        tracks.update(rows_to_tracks(read_csv(pw / f"traj_{cs:05d}.csv")))
        for r in read_csv(pw / f"bio_{cs:05d}.csv"):
            bios[int(r["track_id"])] = r
    return done, tracks, bios


# --------------------------------------------------------------------------------------
# tracking one grid of windows
# --------------------------------------------------------------------------------------

def run_grid(name, out, sub, wins, video, bg, lumen, top, bot, roi, b_ids):
    """Track every not-yet-done window of one grid. Grid A numbers cells sequentially in
    entrance order; grid B uses B_ID_BASE + core_start*100 + k (raw ids, replaced at merge)."""
    pw = out / sub
    pw.mkdir(exist_ok=True)
    done, tracks, _ = load_done(out, wins, sub)
    next_id = max(tracks) + 1 if (tracks and not b_ids) else 1
    todo = [w for w in wins if w["core_start"] not in done]
    print(f"\nGRID {name}: {len(wins)} windows, {len(wins) - len(todo)} already done", flush=True)
    t0_grid, n_new = time.time(), 0
    for wi, win in enumerate(wins, 1):
        cs, ce = win["core_start"], win["core_end"]
        if cs in done:
            continue
        tag = f"[{name} {wi}/{len(wins)}] core {cs}-{ce}"
        t0 = time.time()
        print(f"{tag}: tracking (decoding frames {win['decode_start']}-{win['decode_end']})...", flush=True)
        found = T.detect_and_track_passages(
            video_path=video, start_frame=win["decode_start"],
            frames_to_run=win["decode_end"] - win["decode_start"],
            background=bg, lumen_mask=lumen, top_wall=top, bot_wall=bot, roi_y_bounds=roi)
        found = T.extend_exit_flare(video, found, bg, lumen, top, bot, roi[0], roi[1], verbose=False)
        kept = sorted((p for p in found.values() if cs <= p[0]["frame"] < ce), key=lambda p: p[0]["frame"])
        t_rows, b_rows, ids = [], [], []
        for k, path in enumerate(kept):
            if b_ids:
                gid = B_ID_BASE + cs * 100 + k
            else:
                gid = next_id
                next_id += 1
            for p in path:
                p["track_id"] = gid
            tracks[gid] = path
            bio = T.compute_biomarker(gid, path)
            bio["window_core_start"] = cs
            t_rows += traj_rows(gid, path)
            b_rows.append(bio)
            ids.append(gid)
        secs = round(time.time() - t0, 1)
        info = {"core_start": cs, "core_end": ce, "decode_start": win["decode_start"],
                "decode_end": win["decode_end"], "cells_detected_in_window": len(found),
                "cells_kept": len(kept), "kept_track_ids": " ".join(map(str, ids)), "seconds": secs}
        write_csv_atomic(pw / f"traj_{cs:05d}.csv", TRAJ_FIELDS, t_rows)
        write_csv_atomic(pw / f"bio_{cs:05d}.csv", BIO_FIELDS, b_rows)
        tmp = pw / f"info_{cs:05d}.json.tmp"
        tmp.write_text(json.dumps(info, indent=2), encoding="utf-8")
        os.replace(tmp, pw / f"info_{cs:05d}.json")       # written last = window done
        done[cs] = info
        n_new += 1
        left = len(wins) - len(done)
        eta = (time.time() - t0_grid) / n_new * left / 60.0
        print(f"{tag}: {len(kept)} cells kept (of {len(found)} seen incl. pre/post-roll)  [{secs}s]"
              + (f"  ~{eta:.0f} min left in grid {name}" if left else ""), flush=True)
    return done, tracks


# --------------------------------------------------------------------------------------
# merging grid B into grid A
# --------------------------------------------------------------------------------------

def _bore_pairs(pa, pb):
    """(x_a, x_b) at frames both tracks share, where either is inside the bore."""
    ma = {p["frame"]: p["cx"] for p in pa}
    pairs = []
    for p in pb:
        xa = ma.get(p["frame"])
        if xa is not None and (110 <= xa <= 1180 or 110 <= p["cx"] <= 1180):
            pairs.append((xa, p["cx"]))
    return pairs


def merge_grids(tracks_a, tracks_b, evidence_fn=None):
    """Grid A is kept as is. Each grid-B track is a duplicate of a kept track (median bore
    distance < 20 px over >= 20 shared frames), a conflict (shares > 30 % of its shared bore
    frames within 14 px, or crosses it), rejected (matches nothing but its ridge evidence is
    below RIDGE_MIN_RATIO), or a new cell (added)."""
    def ent(p):
        return p[0]["frame"]
    kept = [(aid, pa, "A") for aid, pa in sorted(tracks_a.items(), key=lambda kv: ent(kv[1]))]
    confirmed, report = set(), []
    for bid, pb in sorted(tracks_b.items(), key=lambda kv: ent(kv[1])):
        row = {"bid": bid, "grid_b_entrance_frame": ent(pb), "grid_b_exit_frame": pb[-1]["frame"],
               "grid_b_window_core": (bid - B_ID_BASE) // 100, "verdict": "added", "partner": None,
               "shared_bore_frames": "", "median_distance_px": "", "ridge_ratio": "",
               "note": "not found by grid A"}
        for oid, po, src in kept:
            pairs = _bore_pairs(po, pb)
            if len(pairs) < 20:
                continue
            d = np.array([a - b for a, b in pairs], dtype=np.float64)
            ad = np.abs(d)
            med = float(np.median(ad))
            if med < 20.0:
                row.update(verdict="duplicate", partner=oid, shared_bore_frames=len(pairs),
                           median_distance_px=round(med, 1), note="same cell as a grid-A track")
                if src == "A":
                    confirmed.add(oid)
                break
            close = float(np.mean(ad < 14.0))
            signs = {int(np.sign(v)) for v in d if abs(v) > 5.0}
            if close > 0.3 or len(signs) > 1:
                row.update(verdict="conflict", partner=oid, shared_bore_frames=len(pairs),
                           median_distance_px=round(med, 1),
                           note=("shares part of a ridge with" if close > 0.3 else "crosses (overtakes)")
                           + " a kept track; not added")
                break
        if row["verdict"] == "added" and evidence_fn is not None:
            ratio = evidence_fn(bid, pb)
            row["ridge_ratio"] = round(ratio, 2)
            if ratio < RIDGE_MIN_RATIO:
                row.update(verdict="rejected",
                           note=f"not found by grid A, but ridge evidence {ratio:.2f} < {RIDGE_MIN_RATIO} (noise level)")
            else:
                row["note"] = f"not found by grid A; ridge evidence {ratio:.2f}"
        if row["verdict"] == "added":
            kept.append((bid, pb, "B"))
        report.append(row)
    kept.sort(key=lambda t: ent(t[1]))
    return kept, confirmed, report


def finalize(out, settings, wins_a, wins_b, done_a, done_b, kept, confirmed, report, entry_fn=None):
    """Renumber by entrance frame, recompute biomarkers, run the entry-funnel check (flag only,
    nothing is removed), write the merged outputs."""
    tracks, bios, new_id = {}, {}, {}
    traj, bio_rows = [], []
    for k, (oid, path, src) in enumerate(kept, 1):
        pts = [dict(p, track_id=k) for p in path]
        bio = T.compute_biomarker(k, pts)
        e = pts[0]["frame"]
        owner = next((w["core_start"] for w in wins_a if w["core_start"] <= e < w["core_end"]), wins_a[-1]["core_start"])
        bio["window_core_start"] = owner
        bio["source"] = "A+B" if (src == "A" and oid in confirmed) else src
        if entry_fn is not None:
            sc = entry_fn(pts, owner)
            bio["entry_cell_score"] = round(sc, 2)
            bio["entry_check"] = entry_label(sc)
        else:
            bio["entry_cell_score"] = ""
            bio["entry_check"] = "not_checked"
        new_id[oid] = k
        tracks[k] = pts
        bios[k] = {kk: str(v) for kk, v in bio.items()}
        traj += traj_rows(k, pts)
        bio_rows.append(bio)
    rep_rows = []
    for r in report:
        rep_rows.append({
            "grid_b_track": r["bid"], "entrance_frame": r["grid_b_entrance_frame"],
            "exit_frame": r["grid_b_exit_frame"], "verdict": r["verdict"],
            "final_track_id": new_id.get(r["bid"], "") if r["verdict"] == "added" else "",
            "matched_grid_a_track": new_id.get(r["partner"], "") if r["partner"] is not None else "",
            "shared_bore_frames": r["shared_bore_frames"], "median_distance_px": r["median_distance_px"],
            "ridge_ratio": r["ridge_ratio"], "note": r["note"]})
    win_rows = [dict(grid="A", **{k: done_a[w["core_start"]][k] for k in WIN_FIELDS}) for w in wins_a if w["core_start"] in done_a]
    win_rows += [dict(grid="B", **{k: done_b[w["core_start"]][k] for k in WIN_FIELDS}) for w in wins_b if w["core_start"] in done_b]
    write_csv_atomic(out / "trajectories.csv", TRAJ_FIELDS, traj)
    write_csv_atomic(out / "biomechanics.csv", FINAL_BIO_FIELDS, bio_rows)
    write_csv_atomic(out / "windows.csv", WIN_OUT_FIELDS, win_rows)
    write_csv_atomic(out / "merge_report.csv", REPORT_FIELDS, rep_rows)

    vels = [b["velocity_mean"] for b in bio_rows if b["velocity_mean"]]
    adis = [b["ADI"] for b in bio_rows if b["ADI"] is not None]
    n_src = {s: sum(1 for b in bio_rows if b["source"] == s) for s in ("A", "A+B", "B")}
    summary = {
        "total_cells": len(bio_rows),
        "cells_by_source": {"grid A only": n_src["A"], "grid A, confirmed by grid B": n_src["A+B"],
                            "added by grid B (missed by grid A)": n_src["B"]},
        "grid_b_tracks": {v: sum(1 for r in report if r["verdict"] == v)
                          for v in ("added", "duplicate", "conflict", "rejected")},
        "ridge_min_ratio_for_grid_b_additions": RIDGE_MIN_RATIO,
        "entry_funnel_check": {
            "cell_seen": sum(1 for b in bio_rows if b["entry_check"] == "cell_seen"),
            "uncertain": sum(1 for b in bio_rows if b["entry_check"] == "uncertain"),
            "no_cell_seen": sum(1 for b in bio_rows if b["entry_check"] == "no_cell_seen"),
            "not_checked": sum(1 for b in bio_rows if b["entry_check"] == "not_checked"),
            "thresholds": {"cell_seen_at_or_above": ENTRY_SEEN, "no_cell_seen_below": ENTRY_NONE},
            "note": "flags only: no track is removed; see entry_check_review.png for flagged tracks",
        },
        "windows_done": {"A": len(done_a), "B": len(done_b)},
        "windows_planned": {"A": len(wins_a), "B": len(wins_b)},
        "cells_per_window": {str(w["core_start"]): sum(1 for b in bio_rows if b["window_core_start"] == w["core_start"])
                             for w in wins_a},
        "settings": settings,
        "metrics": {
            "velocity_mean_px_per_frame": float(np.mean(vels)) if vels else None,
            "ADI_mean": float(np.mean(adis)) if adis else None,
        },
        "not_trackable_by_design": [
            "cells already inside the channel at the first decoded frame",
            "cells entering in roughly the last 260 frames of the recording",
        ],
    }
    tmp = out / "summary.json.tmp"
    tmp.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    os.replace(tmp, out / "summary.json")
    return tracks, bios, summary


def render_all(out, wins_a, tracks, bios, video, bg, lumen, top, bot, roi, no_video, no_kymo):
    """Per grid-A core: video + kymograph. Re-rendered only when the cells drawn in that core
    changed (ids, frame spans, end points), recorded in render_state.json."""
    state_p = out / "render_state.json"
    state = json.loads(state_p.read_text(encoding="utf-8")) if state_p.exists() else {}
    for wi, win in enumerate(wins_a, 1):
        cs, ce = win["core_start"], win["core_end"]
        clipped = clip_tracks(tracks, cs, ce)
        sig = json.dumps(sorted([gid, p[0]["frame"], p[-1]["frame"], round(p[0]["cx"], 1), round(p[-1]["cx"], 1)]
                                for gid, p in clipped.items()))
        own = {gid for gid, b in bios.items() if str(b.get("window_core_start")) == str(cs)}
        vid = out / "videos" / f"cells_frames_{cs:05d}-{ce:05d}.mp4"
        kym = out / "kymographs" / f"kymo_frames_{cs:05d}-{ce:05d}.png"
        tag = f"[render {wi}/{len(wins_a)}] core {cs}-{ce}"
        if not no_video and (state.get(f"video_{cs}") != sig or (clipped and not vid.exists())):
            v = render_window_video(video, top, bot, tracks, bios, win, out)
            state[f"video_{cs}"] = sig
            state_p.write_text(json.dumps(state), encoding="utf-8")
            if v:
                print(f"{tag}: video -> videos/{v.name}", flush=True)
        if not no_kymo and (state.get(f"kymo_{cs}") != sig or not kym.exists()):
            k = kymograph_png(video, bg, lumen, roi, tracks, win, out, own)
            state[f"kymo_{cs}"] = sig
            state_p.write_text(json.dumps(state), encoding="utf-8")
            if k:
                print(f"{tag}: kymograph -> kymographs/{k.name}", flush=True)


# --------------------------------------------------------------------------------------
# per-window QA outputs
# --------------------------------------------------------------------------------------

def clip_tracks(tracks, f0, f1):
    out = {}
    for gid, path in tracks.items():
        pts = [p for p in path if f0 <= p["frame"] < f1]
        if pts:
            out[gid] = pts
    return out


def render_window_video(video, top, bot, tracks, bios, win, out):
    """Tracked video of the window core with every cell present drawn. Uses
    tracking.render_tracked_video, which renders only frames with a cell (+-20)."""
    final = out / "videos" / f"cells_frames_{win['core_start']:05d}-{win['core_end']:05d}.mp4"
    clipped = clip_tracks(tracks, win["core_start"], win["core_end"])
    if not clipped:
        return None
    tmp = out / "videos" / f"_tmp_{win['core_start']}"
    tmp.mkdir(parents=True, exist_ok=True)
    adi = {gid: {"ADI": float(b["ADI"])} for gid, b in bios.items() if b.get("ADI") not in (None, "", "None")}
    mp4, _ = T.render_tracked_video(video, top, bot, clipped, adi, str(tmp), pad_frames=20, playback_fps=20.0)
    os.replace(mp4, final)
    shutil.rmtree(tmp, ignore_errors=True)
    return final


def energy_rows(video, bg, lumen, roi, f0, f1):
    """Per frame f0..f1-1: column maximum of motion x contrast energy inside the lumen (the same
    signal the tracker's gap pass uses), minus each column's median over these frames."""
    y0, y1 = roi
    bg_roi = bg[y0:y1].astype(np.float32)
    m = lumen[y0:y1]
    cap = cv2.VideoCapture(video)
    prev = None
    if f0 > 0:
        cap.set(cv2.CAP_PROP_POS_FRAMES, f0 - 1)
        ok, fr = cap.read()
        if ok:
            prev = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY)[y0:y1].astype(np.float32)
    else:
        cap.set(cv2.CAP_PROP_POS_FRAMES, 0)
    rows = []
    for _ in range(f0, f1):
        ok, fr = cap.read()
        if not ok:
            break
        g = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY)[y0:y1].astype(np.float32)
        if prev is None:
            prev = g
        mot = np.abs(g - prev)
        absd = np.abs(g - bg_roi)
        mot[m == 0] = 0
        absd[m == 0] = 0
        rows.append(cv2.GaussianBlur(mot * np.sqrt(np.maximum(1.0, absd)), (5, 5), 0).max(axis=0))
        prev = g
    cap.release()
    if not rows:
        return None
    K = np.stack(rows)
    return np.maximum(0.0, K - np.median(K, axis=0, keepdims=True))


def ridge_evidence(video, bg, lumen, roi, path, pad=300, n_null=300):
    """Ridge evidence of one track: median 9-column-max energy along the track's bore part
    (150 < x < 1100) divided by the same statistic on random straight lines at bore speeds
    (3.5-6.5 px/frame) over the same frames. Noise ~1; real cells >= ~2.7 (see RIDGE_MIN_RATIO)."""
    f0 = max(1, path[0]["frame"] - pad)
    f1 = path[-1]["frame"] + pad
    K = energy_rows(video, bg, lumen, roi, f0, f1)
    if K is None:
        return 0.0
    K9 = cv2.dilate(K.astype(np.float32), np.ones((1, 9), np.uint8))
    W = K9.shape[1]
    vals = [K9[p["frame"] - f0, int(round(p["cx"]))] for p in path
            if 150 < p["cx"] < 1100 and 0 <= p["frame"] - f0 < K9.shape[0]]
    if not vals:
        return 0.0
    rng = np.random.default_rng(12345)
    nulls = []
    for _ in range(n_null):
        v = rng.uniform(3.5, 6.5)
        x0 = rng.uniform(1000.0, 1100.0)
        n = int((x0 - 150.0) / v)
        if K9.shape[0] - n <= 1:
            break
        t0 = int(rng.integers(0, K9.shape[0] - n))
        tt = np.arange(t0, t0 + n)
        nulls.append(float(np.median(K9[tt, np.clip(np.round(x0 - v * (tt - t0)).astype(int), 0, W - 1)])))
    ref = float(np.median(nulls)) if nulls else 1.0
    return float(np.median(vals)) / max(ref, 1e-3)


class EntryChecker:
    """Is a real cell visible in the ENTRY funnel while the track is there?

    In the entry funnel (x > 1180) cells are round rings, clearly visible in a background
    difference, while noise tracks show an empty funnel. Per grid-A window a local background
    (median of every 8th frame from core_start-300 to core_start+1800) and a per-pixel jitter
    floor (median |difference|) are built on an entry-funnel crop. For up to 8 frames while the
    track is in the funnel: |frame - background| minus the frame's own median (alternate-frame
    flicker), minus 1.5 x jitter floor (tube jitter), blurred to cell size (sigma 4 px); the
    strongest value inside the lumen within +-25 px of the track's x. Score = median over those
    frames. Tracks with < 3 funnel points are checked over the whole funnel 8-40 frames before
    they start."""

    X0 = 1100

    def __init__(self, video, lumen, total):
        self.video, self.total = video, total
        rows = np.where(lumen[:, self.X0:].any(axis=1))[0]
        self.y0, self.y1 = max(0, int(rows.min()) - 8), min(lumen.shape[0], int(rows.max()) + 9)
        self.lum = cv2.dilate((lumen[self.y0:self.y1, self.X0:] > 0).astype(np.uint8), np.ones((5, 5), np.uint8)) > 0
        self.xs = np.arange(self.lum.shape[1]) + self.X0
        self.cap = cv2.VideoCapture(video)
        self.models = {}

    def crop(self, f):
        self.cap.set(cv2.CAP_PROP_POS_FRAMES, int(f))
        ok, im = self.cap.read()
        if not ok:
            return None
        return cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)[self.y0:self.y1, self.X0:].astype(np.float32)

    def model(self, cs):
        if cs not in self.models:
            fr = [self.crop(f) for f in range(max(0, cs - 300), min(self.total - 1, cs + 1800), 8)]
            st = np.stack([x for x in fr if x is not None])
            bg = np.median(st, axis=0)
            dd = []
            for x in st[::2]:
                d = x - bg
                d -= np.median(d)
                dd.append(np.abs(cv2.GaussianBlur(d, (3, 3), 0)))
            self.models = {cs: (bg, np.median(np.stack(dd), axis=0))}   # keep only the current window
        return self.models[cs]

    def blob(self, g, bg, floor):
        d = g - bg
        d -= np.median(d)
        e = np.maximum(0.0, np.abs(cv2.GaussianBlur(d, (3, 3), 0)) - 1.5 * floor)
        return cv2.GaussianBlur(e, (0, 0), 4.0)

    def score(self, path, core_start):
        bg, floor = self.model(core_start)
        funnel = [p for p in path if 1185 < p["cx"] < 1276]
        if len(funnel) >= 3:
            idx = np.linspace(0, len(funnel) - 1, min(8, len(funnel))).astype(int)
            checks = [(funnel[i]["frame"], funnel[i]["cx"]) for i in idx]
        else:
            checks = [(path[0]["frame"] - back, None) for back in (8, 16, 24, 32, 40)]
        vals = []
        for f, x in checks:
            g = self.crop(f)
            if g is None:
                continue
            b = self.blob(g, bg, floor)
            region = self.lum & ((np.abs(self.xs - x) <= 25) if x is not None else (self.xs > 1185))[None, :]
            if region.any():
                vals.append(float(b[region].max()))
        return float(np.median(vals)) if vals else 0.0

    def review_tile(self, path, core_start, label):
        """Two background-difference crops of the entry funnel (track position ticked) + label."""
        bg, _ = self.model(core_start)
        funnel = [p for p in path if 1185 < p["cx"] < 1276] or path[:1]
        tiles = []
        for p in (funnel[0], funnel[len(funnel) // 2]):
            g = self.crop(p["frame"])
            if g is None:
                g = bg.copy()
            d = np.clip(128 + 3 * (g - bg - np.median(g - bg)), 0, 255).astype(np.uint8)
            d = cv2.cvtColor(cv2.resize(d, None, fx=2, fy=2, interpolation=cv2.INTER_NEAREST), cv2.COLOR_GRAY2BGR)
            px = int(round((p["cx"] - self.X0) * 2))
            cv2.line(d, (px, 0), (px, 8), (0, 0, 255), 3)
            tiles.append(d)
        t = np.hstack([tiles[0], np.full((tiles[0].shape[0], 4, 3), 255, np.uint8), tiles[1]])
        bar = np.full((18, t.shape[1], 3), 30, np.uint8)
        cv2.putText(bar, label, (4, 13), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (0, 255, 255), 1, cv2.LINE_AA)
        return np.vstack([t, bar])


def entry_label(score):
    if score >= ENTRY_SEEN:
        return "cell_seen"
    return "no_cell_seen" if score < ENTRY_NONE else "uncertain"


def kymograph_png(video, bg, lumen, roi, tracks, win, out, own_ids):
    """Column-max motion x contrast energy over the core frames (time down, x across) minus
    each column's median, with every cell's track drawn: this window's cells coloured and
    labelled, cells from the previous window grey."""
    f0, f1 = win["core_start"], win["core_end"]
    K = energy_rows(video, bg, lumen, roi, f0, f1)
    if K is None:
        return None
    rows = K
    img = cv2.cvtColor(np.clip(K / 25.0 * 255.0, 0, 255).astype(np.uint8), cv2.COLOR_GRAY2BGR)
    for gid, path in clip_tracks(tracks, f0, f0 + len(rows)).items():
        col = T._cell_color_bgr(gid) if gid in own_ids else (150, 150, 150)
        pts = np.array([[int(round(p["cx"])), p["frame"] - f0] for p in path], np.int32)
        cv2.polylines(img, [pts], False, col, 2, cv2.LINE_AA)
        x0, yy = pts[0]
        cv2.putText(img, f"#{gid}", (int(min(max(x0 - 40, 2), img.shape[1] - 60)), int(min(yy + 14, img.shape[0] - 4))),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.5, col, 1, cv2.LINE_AA)
    cv2.putText(img, f"frames {f0}-{f0 + len(rows)}  time down, x across (entrance right, exit left)",
                (8, 18), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 255), 1, cv2.LINE_AA)
    final = out / "kymographs" / f"kymo_frames_{f0:05d}-{f1:05d}.png"
    cv2.imwrite(str(final), img)
    return final


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Track every RBC passage in the whole video (windowed).")
    ap.add_argument("--out", type=str, default=None, help="Output folder (default: Desktop/RBC_full_video_results)")
    ap.add_argument("--video", type=str, default=None, help="Override AVI path")
    ap.add_argument("--start", type=int, default=0, help="First core frame (default 0)")
    ap.add_argument("--end", type=int, default=None, help="End of last core, exclusive (default: last frame)")
    ap.add_argument("--core", type=int, default=1500)
    ap.add_argument("--pre", type=int, default=450)
    ap.add_argument("--post", type=int, default=450)
    ap.add_argument("--resume", action="store_true", help="Continue an interrupted run in the same folder")
    ap.add_argument("--overwrite", action="store_true", help="Delete this script's previous outputs in --out first")
    ap.add_argument("--no-video", action="store_true", help="Skip per-window videos")
    ap.add_argument("--no-kymo", action="store_true", help="Skip per-window kymograph PNGs")
    ap.add_argument("--no-offset-grid", action="store_true",
                    help="Skip the second window grid (offset by core/2) that catches faint cells grid A misses")
    args = ap.parse_args()

    out = Path(args.out) if args.out else default_out_dir()
    out.mkdir(parents=True, exist_ok=True)
    if args.overwrite:
        for name in OUR_FILES:
            if (out / name).exists():
                (out / name).unlink()
        for name in OUR_DIRS:
            shutil.rmtree(out / name, ignore_errors=True)
    has_previous = (out / "settings.json").exists()
    if has_previous and not args.resume:
        sys.exit(f"{out} already has results. Use --resume to continue that run, "
                 f"--overwrite to start again, or --out for a new folder.")
    for d in OUR_DIRS:
        (out / d).mkdir(exist_ok=True)
    sys.stdout = Tee(out / "run_log.txt", sys.stdout)

    if not os.path.isfile(T.CALIB_FILE) or not os.path.isfile(T.BG_FILE):
        sys.exit(f"Stage 1 outputs not found ({T.CALIB_FILE}, {T.BG_FILE}). Run preprocessing.py first.")
    calib = T.load_calibration()
    bg = T.load_background()
    video = args.video or T._resolve_video_path(calib.get("video"))
    if not video or not os.path.isfile(video):
        sys.exit(f"Video not found: {video!r}. Pass --video <path to the .avi>.")
    cap = cv2.VideoCapture(video)
    total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    cap.release()
    end = total if args.end is None else min(args.end, total)

    settings = {"video": video, "start": args.start, "end": end, "core": args.core,
                "pre": args.pre, "post": args.post, "total_frames_in_video": total,
                "offset_grid": 0 if args.no_offset_grid else args.core // 2}
    if has_previous:
        old = json.loads((out / "settings.json").read_text(encoding="utf-8"))
        if any(old.get(k) != settings[k] for k in SETTING_KEYS):
            sys.exit(f"--resume: settings differ from the saved run in {out}:\n  saved {old}\n  now   {settings}")
    else:
        (out / "settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")

    lumen, top, bot, _, _ = T.build_exact_lumen_geometry(bg)
    ys = np.where(lumen > 0)[0]
    roi = (max(0, int(ys.min()) - 6), min(bg.shape[0], int(ys.max()) + 6))
    wins_a = plan_windows(total, args.start, end, args.core, args.pre, args.post)
    offset = 0 if args.no_offset_grid else args.core // 2
    wins_b = (plan_windows(total, args.start + offset, end, args.core, args.pre, args.post)
              if offset and args.start + offset < end else [])

    print("=" * 76)
    print(f"WHOLE-VIDEO RBC TRACKING | frames {args.start}-{end} of {total}")
    print(f"grid A: {len(wins_a)} windows (cores from {args.start})"
          + (f" | grid B: {len(wins_b)} windows (cores offset by {offset})" if wins_b else " | offset grid off"))
    print(f"Video : {video}\nOutput: {out}")
    print("=" * 76, flush=True)

    t_run = time.time()
    done_a, tracks_a = run_grid("A", out, "per_window", wins_a, video, bg, lumen, top, bot, roi, b_ids=False)
    done_b, tracks_b = ({}, {})
    if wins_b:
        done_b, tracks_b = run_grid("B", out, "per_window_b", wins_b, video, bg, lumen, top, bot, roi, b_ids=True)

    score_p = out / "ridge_scores.json"
    scores = json.loads(score_p.read_text(encoding="utf-8")) if score_p.exists() else {}

    def evidence(bid, path):
        key = f"{bid}:{path[0]['frame']}:{path[-1]['frame']}"
        if key not in scores:
            scores[key] = ridge_evidence(video, bg, lumen, roi, path)
            score_p.write_text(json.dumps(scores, indent=1), encoding="utf-8")
            print(f"  ridge evidence for grid-B candidate entering {path[0]['frame']}: {scores[key]:.2f}", flush=True)
        return scores[key]

    kept, confirmed, report = merge_grids(tracks_a, tracks_b, evidence_fn=evidence)

    checker = EntryChecker(video, lumen, total)
    entry_p = out / "entry_scores.json"
    entry_cache = json.loads(entry_p.read_text(encoding="utf-8")) if entry_p.exists() else {}
    print("\nENTRY-FUNNEL CHECK (is a real cell visible where each track enters?)", flush=True)

    def entry_fn(path, owner):
        key = f"{path[0]['frame']}:{path[-1]['frame']}:{round(path[0]['cx'], 1)}"
        if key not in entry_cache:
            entry_cache[key] = checker.score(path, owner)
            entry_p.write_text(json.dumps(entry_cache, indent=1), encoding="utf-8")
        return entry_cache[key]

    tracks, bios, summary = finalize(out, settings, wins_a, wins_b, done_a, done_b, kept, confirmed, report,
                                     entry_fn=entry_fn)
    ec = summary["entry_funnel_check"]
    print(f"ENTRY CHECK: {ec['cell_seen']} cell seen | {ec['uncertain']} uncertain | "
          f"{ec['no_cell_seen']} no cell seen (flags only; nothing removed)", flush=True)
    flagged = sorted((float(b["entry_cell_score"]), gid) for gid, b in bios.items() if b["entry_check"] != "cell_seen")
    review_p = out / "entry_check_review.png"
    if flagged:
        tiles = []
        for sc, gid in flagged:
            b = bios[gid]
            tiles.append(checker.review_tile(
                tracks[gid], int(b["window_core_start"]),
                f"#{gid} enters {b['entrance_frame']} | score {sc:.2f} | {b['entry_check']} | {b['source']}"))
        while len(tiles) % 3:
            tiles.append(np.full_like(tiles[0], 40))
        sheet = np.vstack([np.hstack([np.pad(t, ((3, 3), (3, 3), (0, 0)), constant_values=255) for t in tiles[i:i + 3]])
                           for i in range(0, len(tiles), 3)])
        cv2.imwrite(str(review_p), sheet)
        print(f"review sheet of {len(flagged)} flagged tracks -> {review_p.name}", flush=True)
    elif review_p.exists():
        review_p.unlink()
    src = summary["cells_by_source"]
    print(f"\nMERGE: {summary['total_cells']} cells = {src['grid A only']} grid A only + "
          f"{src['grid A, confirmed by grid B']} confirmed by grid B + "
          f"{src['added by grid B (missed by grid A)']} added by grid B "
          f"| grid-B tracks: {summary['grid_b_tracks']}  (details: merge_report.csv)", flush=True)

    render_all(out, wins_a, tracks, bios, video, bg, lumen, top, bot, roi, args.no_video, args.no_kymo)

    print("\n" + "=" * 76)
    print(f"FINISHED: {summary['total_cells']} cells over frames {args.start}-{end} "
          f"in {(time.time() - t_run) / 60:.1f} min. Results in {out}")
    print("=" * 76, flush=True)


if __name__ == "__main__":
    main()
