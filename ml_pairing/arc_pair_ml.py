import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np


@dataclass
class ArcSegment:
    points: np.ndarray  # [N, 2], (x, y)
    center: np.ndarray
    length: float
    curvature: float
    angle_span: float
    gt_id: int


def read_manifest(manifest_path: str) -> List[Dict]:
    with open(manifest_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, dict) and "samples" in data:
        return data["samples"]
    if isinstance(data, list):
        return data
    raise ValueError("Manifest must be a list or object with 'samples'.")


def point_to_ellipse_residual(points: np.ndarray, ellipse: Sequence[float]) -> np.ndarray:
    cx, cy, a, b, phi_deg = ellipse
    phi = math.radians(phi_deg)
    cp, sp = math.cos(phi), math.sin(phi)
    x = points[:, 0] - cx
    y = points[:, 1] - cy
    xr = cp * x + sp * y
    yr = -sp * x + cp * y
    val = np.sqrt((xr / max(a, 1e-6)) ** 2 + (yr / max(b, 1e-6)) ** 2)
    return np.abs(val - 1.0)


def extract_arc_segments(
    gray: np.ndarray,
    gt_ellipses: Sequence[Sequence[float]],
    min_points: int = 18,
    residual_thresh: float = 0.12,
) -> List[ArcSegment]:
    edges = cv2.Canny(gray, 70, 150)
    num, labels = cv2.connectedComponents((edges > 0).astype(np.uint8), connectivity=8)

    gx = cv2.Sobel(gray, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(gray, cv2.CV_32F, 0, 1, ksize=3)

    arcs: List[ArcSegment] = []
    for cc in range(1, num):
        ys, xs = np.where(labels == cc)
        if len(xs) < min_points:
            continue
        pts = np.stack([xs, ys], axis=1).astype(np.float32)
        center = pts.mean(axis=0)

        order = np.argsort(np.arctan2(pts[:, 1] - center[1], pts[:, 0] - center[0]))
        pts = pts[order]

        diffs = np.diff(pts, axis=0)
        seg_len = np.sqrt((diffs ** 2).sum(axis=1)).sum()
        if seg_len < min_points:
            continue

        vectors = pts - center
        ang = np.arctan2(vectors[:, 1], vectors[:, 0])
        angle_span = float(np.unwrap(ang).ptp())

        curvature = float(angle_span / max(seg_len, 1e-6))

        best_gt, best_res = -1, 1e9
        for gt_i, ell in enumerate(gt_ellipses):
            res = point_to_ellipse_residual(pts, ell).mean()
            if res < best_res:
                best_res = res
                best_gt = gt_i

        if best_res > residual_thresh:
            best_gt = -1

        arcs.append(
            ArcSegment(
                points=pts,
                center=center,
                length=float(seg_len),
                curvature=curvature,
                angle_span=angle_span,
                gt_id=best_gt,
            )
        )

    return arcs


def arc_local_image_features(gray: np.ndarray, arc: ArcSegment) -> np.ndarray:
    gx = cv2.Sobel(gray, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(gray, cv2.CV_32F, 0, 1, ksize=3)
    grad_mag = np.sqrt(gx * gx + gy * gy)

    h, w = gray.shape[:2]
    xs = np.clip(arc.points[:, 0].astype(np.int32), 0, w - 1)
    ys = np.clip(arc.points[:, 1].astype(np.int32), 0, h - 1)

    vals = grad_mag[ys, xs]
    p90 = np.percentile(vals, 90) if len(vals) else 0.0
    p50 = np.percentile(vals, 50) if len(vals) else 0.0

    edge_conf = float((vals > 25).mean()) if len(vals) else 0.0
    return np.array([float(vals.mean()), float(vals.std()), float(p50), float(p90), edge_conf], dtype=np.float32)


def pair_features(a: ArcSegment, b: ArcSegment, gray: np.ndarray) -> np.ndarray:
    dist = float(np.linalg.norm(a.center - b.center))
    geom = np.array(
        [
            dist,
            abs(a.length - b.length),
            a.length + b.length,
            abs(a.curvature - b.curvature),
            abs(a.angle_span - b.angle_span),
        ],
        dtype=np.float32,
    )
    img_a = arc_local_image_features(gray, a)
    img_b = arc_local_image_features(gray, b)
    img = np.concatenate([np.abs(img_a - img_b), 0.5 * (img_a + img_b)], axis=0)
    return np.concatenate([geom, img], axis=0)


def make_training_samples(
    manifest_path: str,
    max_neighbors: int = 10,
    neg_dist_scale: float = 0.25,
) -> Tuple[np.ndarray, np.ndarray, List[Dict]]:
    samples = read_manifest(manifest_path)
    feat_list, label_list = [], []
    meta = []

    for item in samples:
        image_path = Path(item["image"])
        gt_ellipses = item["gt_ellipses"]
        gray = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)
        if gray is None:
            raise FileNotFoundError(f"Cannot read image: {image_path}")

        arcs = extract_arc_segments(gray, gt_ellipses)
        if len(arcs) < 2:
            continue

        centers = np.stack([a.center for a in arcs], axis=0)
        diag = float(np.linalg.norm([gray.shape[0], gray.shape[1]]))
        max_dist = diag * neg_dist_scale

        for i, arc_i in enumerate(arcs):
            d = np.linalg.norm(centers - arc_i.center[None, :], axis=1)
            nn = np.argsort(d)[1 : max_neighbors + 1]
            for j in nn:
                if j <= i:
                    continue
                arc_j = arcs[j]
                if d[j] > max_dist:
                    continue

                y = int(arc_i.gt_id >= 0 and arc_i.gt_id == arc_j.gt_id)
                if y == 0 and (arc_i.gt_id < 0 or arc_j.gt_id < 0):
                    continue

                x = pair_features(arc_i, arc_j, gray)
                feat_list.append(x)
                label_list.append(y)
                meta.append({"image": str(image_path), "arc_i": i, "arc_j": int(j), "label": y})

    if not feat_list:
        return np.zeros((0, 15), dtype=np.float32), np.zeros((0,), dtype=np.int64), meta

    return np.stack(feat_list, axis=0), np.array(label_list, dtype=np.int64), meta


class TinyMLP:
    def __init__(self, in_dim: int, hidden_dim: int = 32, seed: int = 0):
        rng = np.random.default_rng(seed)
        self.w1 = (rng.standard_normal((in_dim, hidden_dim)) * 0.1).astype(np.float32)
        self.b1 = np.zeros((hidden_dim,), dtype=np.float32)
        self.w2 = (rng.standard_normal((hidden_dim, 1)) * 0.1).astype(np.float32)
        self.b2 = np.zeros((1,), dtype=np.float32)

    def forward_logits(self, x: np.ndarray) -> np.ndarray:
        h = np.maximum(0.0, x @ self.w1 + self.b1)
        return h @ self.w2 + self.b2

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        z = self.forward_logits(x)
        return 1.0 / (1.0 + np.exp(-z))

    def fit(self, x: np.ndarray, y: np.ndarray, epochs: int = 60, lr: float = 1e-2) -> Dict[str, float]:
        y = y.astype(np.float32).reshape(-1, 1)
        n = max(len(x), 1)
        for _ in range(epochs):
            h_pre = x @ self.w1 + self.b1
            h = np.maximum(0.0, h_pre)
            z = h @ self.w2 + self.b2
            p = 1.0 / (1.0 + np.exp(-z))

            dz = (p - y) / n
            dw2 = h.T @ dz
            db2 = dz.sum(axis=0)
            dh = dz @ self.w2.T
            dh[h_pre <= 0] = 0
            dw1 = x.T @ dh
            db1 = dh.sum(axis=0)

            self.w2 -= lr * dw2
            self.b2 -= lr * db2
            self.w1 -= lr * dw1
            self.b1 -= lr * db1

        p = self.predict_proba(x)
        pred = (p[:, 0] >= 0.5).astype(np.int64)
        acc = float((pred == y[:, 0]).mean())
        return {"train_acc": acc}

    def to_json(self, path: str, mean: np.ndarray, std: np.ndarray) -> None:
        obj = {
            "w1": self.w1.tolist(),
            "b1": self.b1.tolist(),
            "w2": self.w2.tolist(),
            "b2": self.b2.tolist(),
            "mean": mean.tolist(),
            "std": std.tolist(),
        }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(obj, f)

    @staticmethod
    def from_json(path: str) -> Tuple["TinyMLP", np.ndarray, np.ndarray]:
        with open(path, "r", encoding="utf-8") as f:
            obj = json.load(f)
        m = TinyMLP(len(obj["mean"]), len(obj["b1"]))
        m.w1 = np.array(obj["w1"], dtype=np.float32)
        m.b1 = np.array(obj["b1"], dtype=np.float32)
        m.w2 = np.array(obj["w2"], dtype=np.float32)
        m.b2 = np.array(obj["b2"], dtype=np.float32)
        mean = np.array(obj["mean"], dtype=np.float32)
        std = np.array(obj["std"], dtype=np.float32)
        return m, mean, std


def compute_metrics(y_true: np.ndarray, y_score: np.ndarray, threshold: float = 0.5) -> Dict[str, float]:
    y_pred = (y_score >= threshold).astype(np.int64)
    tp = int(((y_pred == 1) & (y_true == 1)).sum())
    fp = int(((y_pred == 1) & (y_true == 0)).sum())
    fn = int(((y_pred == 0) & (y_true == 1)).sum())

    recall = tp / max(tp + fn, 1)
    false_positive_rate = fp / max((y_true == 0).sum(), 1)
    precision = tp / max(tp + fp, 1)
    return {
        "recall": float(recall),
        "false_positive_rate": float(false_positive_rate),
        "precision": float(precision),
    }


def benchmark_rule_vs_mlp(x: np.ndarray, y: np.ndarray, y_score: np.ndarray) -> Dict[str, Dict[str, float]]:
    rule_score = 1.0 / (1.0 + x[:, 0])

    t0 = time.perf_counter()
    rule_metrics = compute_metrics(y, rule_score, threshold=0.2)
    t1 = time.perf_counter()

    mlp_metrics = compute_metrics(y, y_score, threshold=0.5)
    t2 = time.perf_counter()

    rule_metrics["avg_infer_ms"] = float((t1 - t0) * 1000)
    mlp_metrics["avg_infer_ms"] = float((t2 - t1) * 1000)

    return {"rule": rule_metrics, "mlp": mlp_metrics}
