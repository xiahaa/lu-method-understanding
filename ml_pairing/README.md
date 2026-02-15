# Arc-Pair 轻量学习器（MLP）

该目录提供一个**不改动最终几何拟合**的轻量模型流程，用于预测“弧段对是否同属一椭圆”，并在推理阶段用于候选重排（或替代部分阈值打分规则）。

## 能力覆盖

1. **自动样本生成**：从图像 + GT 椭圆自动提取弧段并构造弧段对。
   - 正样本：同一 GT 椭圆上的两段弧。
   - 负样本：空间邻近但 GT 不一致的弧段对。
2. **特征**：
   - 几何特征：中心距离、长度差/和、曲率差、角跨度差。
   - 局部图像特征：梯度均值/方差/分位数、边缘置信度（强梯度占比）。
3. **推理定位**：仅输出弧段对 score 做重排，不替代后续 ellipse fitting。
4. **对比评估**：内置 `rule vs mlp` 指标输出（速度、召回、误检率、精度）。

## 数据清单格式

`manifest.json` 支持两种形式：

- 列表：`[{"image": "...", "gt_ellipses": [[cx,cy,a,b,phi_deg], ...]}, ...]`
- 字典：`{"samples": [...同上...]}`

## 快速开始

```bash
# 1) 生成可运行的合成数据（可替换为真实数据清单）
python3 ml_pairing/generate_synthetic_manifests.py

# 2) 训练并输出基线对比报告
python3 ml_pairing/train_arc_pair_mlp.py \
  --train_manifest ml_pairing/synth_data/train_manifest.json \
  --test_manifest ml_pairing/synth_data/test_manifest.json \
  --out_model ml_pairing/arc_pair_mlp.json \
  --out_report ml_pairing/benchmark_report.json

# 3) 在推理端做候选重排（仅 score/ranking）
python3 ml_pairing/rerank_pairs.py \
  --image ml_pairing/synth_data/test_000.png \
  --model ml_pairing/arc_pair_mlp.json \
  --out ml_pairing/reranked_pairs.json
```

## 跨数据集泛化建议

建议至少做两个 domain：
- 训练：干净边缘数据集（如工业件/合成图）
- 测试：自然场景或噪声更高数据集

并分别报告 `rule` 与 `mlp`：
- `avg_infer_ms`
- `recall`
- `false_positive_rate`
- `precision`

> 该模块是可插拔的排序器；若需要接入 C++ 主流程，可在候选 pair 生成处读取模型分数作为额外排序项。
