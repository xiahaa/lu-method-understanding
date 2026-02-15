# 论文与实验重排方案（LU 椭圆检测）

本方案用于将现有论文叙述与实验章节重排为“模块化贡献 + 可解释评测”的结构，突出 **grouping** 与 **validation** 的独立作用，并补齐复杂场景与效率剖析。

## 1. 单独报告 grouping quality 指标

在主实验中单列“Grouping Quality”小节，不与最终检测 F-measure 混排。建议至少包含：

- **Cluster Purity**：同一簇内边缘/弧段是否来自同一 GT 椭圆。
- **Over-merge rate**：一个簇混入多个 GT 椭圆实例的比例。
- **Over-split rate**：同一 GT 椭圆被拆分到多个簇的比例。

建议展示形式：

- 主表（平均值 ± 方差）按数据集给出 purity / over-merge / over-split。
- 附录给出按尺寸区间（小/中/大椭圆）与遮挡程度分桶统计。

## 2. 分离展示三类增益来源

在方法对比中明确拆分三种变体，避免“整体提升但来源不明”：

- **Only Grouping**：仅替换 grouping，validation 保持 baseline。
- **Only Validation**：仅替换 validation，grouping 保持 baseline。
- **Grouping + Validation**：两者同时替换。

推荐消融表结构如下：

| Variant | Grouping | Validation | Precision | Recall | F1 |
|---|---|---|---:|---:|---:|
| Baseline | Base | Base | - | - | - |
| Only Grouping | New | Base | - | - | - |
| Only Validation | Base | New | - | - | - |
| Grouping + Validation | New | New | - | - | - |

并在文字中报告：

- 绝对增益（ΔF1）
- 相对增益（%）
- 是否存在互补效应（联合改进是否大于单改进之和）

## 3. 增加复杂场景子集评测

在现有数据集中新增困难子集标注并单独汇报结果：

- **Occlusion**（遮挡）
- **Low contrast**（低对比）
- **Repeated texture**（重复纹理）

建议流程：

1. 依据可复现规则划分子集（阈值或人工复核）。
2. 在主文提供每个子集的样本数量与占比。
3. 单独报告 Precision / Recall / F1 与 grouping quality 指标。
4. 给出失败案例图，解释对应模块瓶颈（grouping 误聚类或 validation 误拒绝）。

## 4. 增加运行时间剖析

将总时延拆解为可解释分项，在相同硬件与线程设置下报告：

- **Edge Extraction**
- **Grouping**
- **Validation**

建议报告：

- 每张图平均耗时（ms）
- 分项占比（%）
- 随分辨率变化的曲线（如 640p / 1080p / 2K）

推荐表结构：

| Resolution | Edge Extraction (ms) | Grouping (ms) | Validation (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| 640p | - | - | - | - |
| 1080p | - | - | - | - |
| 2K | - | - | - | - |

## 5. 结论重写重点

结论部分应强调方法贡献的 **范式普适性**，而非对单一实现绑定：

- 本方法作用于“候选组织（grouping）+ 候选验证（validation）”两阶段接口。
- 可插拔于传统椭圆检测流水线（无论候选来自 ASLS、边缘片段或其他几何先验）。
- 实验显示其收益在多数据集与复杂场景子集上稳定，而非仅在某特定代码实现有效。

可直接使用的结论句式（示例）：

> Our improvements are not implementation-specific tweaks, but transferable upgrades to the classical proposal-grouping-validation paradigm for ellipse detection.

## 推荐章节重排（论文）

1. **Method**：先给接口层定义，再分别介绍 grouping / validation。
2. **Experimental Protocol**：数据集、子集划分、指标（含 grouping quality）、硬件配置。
3. **Main Results**：整体指标 + 三类消融（Only Grouping / Only Validation / Both）。
4. **Stress Tests**：遮挡、低对比、重复纹理子集。
5. **Efficiency Analysis**：分项时间剖析。
6. **Conclusion**：强调范式级可迁移性。

---

如需落地到本仓库，可在后续步骤补充：

- 脚本化统计 cluster purity / over-merge / over-split；
- 在测试输出中增加 edge/grouping/validation 分项计时；
- 在 `tests/` 中加入困难子集配置与批量评测入口。
