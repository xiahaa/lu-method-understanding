import argparse
import json
from pathlib import Path

import cv2
import numpy as np


def gen_sample(path: Path, n_ellipses: int, noise: int, seed: int):
    rng = np.random.default_rng(seed)
    h, w = 512, 512
    img = np.zeros((h, w), dtype=np.uint8)
    gts = []
    for _ in range(n_ellipses):
        cx = int(rng.integers(80, w - 80))
        cy = int(rng.integers(80, h - 80))
        a = int(rng.integers(30, 110))
        b = int(rng.integers(20, min(90, a)))
        phi = float(rng.integers(0, 180))
        cv2.ellipse(img, (cx, cy), (a, b), phi, 0, 360, 255, 2)
        gts.append([float(cx), float(cy), float(a), float(b), float(phi)])

    noise_pts = rng.integers(0, min(h, w), size=(noise, 2))
    for y, x in noise_pts:
        img[y, x] = 255

    cv2.imwrite(str(path), img)
    return {"image": str(path), "gt_ellipses": gts}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", default="ml_pairing/synth_data")
    ap.add_argument("--train_n", type=int, default=24)
    ap.add_argument("--test_n", type=int, default=10)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    train, test = [], []
    for i in range(args.train_n):
        train.append(gen_sample(out / f"train_{i:03d}.png", n_ellipses=3, noise=300, seed=10 + i))
    for i in range(args.test_n):
        test.append(gen_sample(out / f"test_{i:03d}.png", n_ellipses=4, noise=500, seed=100 + i))

    with open(out / "train_manifest.json", "w", encoding="utf-8") as f:
        json.dump(train, f, indent=2)
    with open(out / "test_manifest.json", "w", encoding="utf-8") as f:
        json.dump(test, f, indent=2)

    print(out / "train_manifest.json")
    print(out / "test_manifest.json")


if __name__ == "__main__":
    main()
