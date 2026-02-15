import argparse
import json
from pathlib import Path

import cv2
import numpy as np

from arc_pair_ml import TinyMLP, extract_arc_segments, pair_features


def main() -> None:
    ap = argparse.ArgumentParser(description="Use trained MLP to rerank arc-pair candidates.")
    ap.add_argument("--image", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--gt_ellipses", default="", help="Optional GT ellipses json list for arc assignment/extraction stability.")
    ap.add_argument("--top_k", type=int, default=100)
    ap.add_argument("--out", default="ml_pairing/reranked_pairs.json")
    args = ap.parse_args()

    gray = cv2.imread(args.image, cv2.IMREAD_GRAYSCALE)
    if gray is None:
        raise FileNotFoundError(args.image)

    gt_ellipses = []
    if args.gt_ellipses:
        with open(args.gt_ellipses, "r", encoding="utf-8") as f:
            gt_ellipses = json.load(f)

    arcs = extract_arc_segments(gray, gt_ellipses, residual_thresh=10.0)
    model, mean, std = TinyMLP.from_json(args.model)

    cand = []
    for i in range(len(arcs)):
        for j in range(i + 1, len(arcs)):
            x = pair_features(arcs[i], arcs[j], gray)
            xn = (x - mean) / std
            score = float(model.predict_proba(xn[None, :])[0, 0])
            cand.append({"arc_i": i, "arc_j": j, "score": score})

    cand.sort(key=lambda z: z["score"], reverse=True)
    result = {
        "image": args.image,
        "num_arcs": len(arcs),
        "top_pairs": cand[: args.top_k],
        "note": "MLP scores are for candidate reranking only; final ellipse geometry fitting stays unchanged.",
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    print(f"Saved {len(result['top_pairs'])} pairs to {args.out}")


if __name__ == "__main__":
    main()
