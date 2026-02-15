import argparse
import json
from pathlib import Path

import numpy as np

from arc_pair_ml import TinyMLP, benchmark_rule_vs_mlp, make_training_samples


def main() -> None:
    ap = argparse.ArgumentParser(description="Train lightweight MLP arc-pair classifier.")
    ap.add_argument("--train_manifest", required=True)
    ap.add_argument("--test_manifest", required=True)
    ap.add_argument("--out_model", default="ml_pairing/arc_pair_mlp.json")
    ap.add_argument("--out_report", default="ml_pairing/benchmark_report.json")
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--hidden", type=int, default=32)
    args = ap.parse_args()

    x_train, y_train, _ = make_training_samples(args.train_manifest)
    x_test, y_test, _ = make_training_samples(args.test_manifest)

    if len(x_train) == 0 or len(x_test) == 0:
        raise RuntimeError("Empty dataset after sample generation. Check manifests and GT ellipses.")

    mean = x_train.mean(axis=0)
    std = x_train.std(axis=0) + 1e-6
    x_train_n = (x_train - mean) / std
    x_test_n = (x_test - mean) / std

    model = TinyMLP(in_dim=x_train_n.shape[1], hidden_dim=args.hidden)
    train_stats = model.fit(x_train_n, y_train, epochs=args.epochs)

    y_score = model.predict_proba(x_test_n)[:, 0]
    report = benchmark_rule_vs_mlp(x_test, y_test, y_score)
    report["train"] = train_stats
    report["train_size"] = int(len(x_train))
    report["test_size"] = int(len(x_test))

    Path(args.out_model).parent.mkdir(parents=True, exist_ok=True)
    model.to_json(args.out_model, mean, std)

    Path(args.out_report).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_report, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
