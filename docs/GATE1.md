# 关卡一页纸：取数前待确认件（确认后清单冻结）

## 1. 通道验证（已做，4 个校验值，SHA256）

| 文件 | 校验值 | 来源 |
|---|---|---|
| DMS_substitutions.csv（217 测定清单参考文件） | `a8f498…0308` | OATML-Markslab/ProteinGym @ commit `144fe22b07dfaeec2b366f2346203a9838a55b4c`（v1.3）reference_files |
| DMS_substitutions_Spearman_DMS_level.csv（A 段对照对象） | `f432423b…66437` | 同 commit benchmarks/DMS_zero_shot/substitutions/Spearman |
| ModelScope facebook/esm2_t33_650M_UR50D config.json | `539095c2…96b2b` | modelscope.cn master |
| 同上 vocab.txt | `0b82cc0a…598e03` | 同上 |

主数据通道（确认后下载）：`https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip`（v1.3，217 测定逐突变 CSV，含 mutated_sequence、DMS_score）。
权重通道：`modelscope.cn facebook/esm2_t33_650M_UR50D`（EsmForMaskedLM，651M 参数，已核 config/vocab）。两通道均已实测连通。

## 2. 测定清单及筛选规则（基于已下载参考文件实算）

- R1 仅替换基准（substitutions，217 测定全体）。
- R2 `seq_len ≤ 1022`（ESM2 单次前向上限；超长截窗会把窗口选择混入「模型内部信号」，弃 16 个）。
- R3 `DMS_total_number_mutants ≥ 500`（残基级统计与 bootstrap 需要样本量，弃 15 个）。
- 交集 **191 个测定**，合计 **2,411,093 个突变**，含多点突变测定 69 个。
- 清单全文：`data/assay_list_proposed.csv`（SHA256 `a21b607e…37b08`）。
- 墙钟对策：单卡 4090 / 6h；masked-marginals 前向上界 1,822,479 次（单点按唯一位置去重计，多点每个突变位点各一次，实际只会更多，取数后立即精确化；**该估计已作废**：官方口径为 WT 序列逐位置遮蔽，191 测定实际共 54,182 次前向，见 `docs/STATE.md` RSI-2 更正）。若超时，按 `DMS_index` 顺序处理，未跑完的测定如实列出，不缩小清单。

## 3. 「误差」的定义

- 逐突变残差 `r_m = z(DMS_score)_m − z(模型分)_m`，z 按测定内标准化（双侧）。
- 残基级误差 `e_i = mean(|r_m|)`，m 为落在残基 i 上的全部已测突变（测定内计算）；
  突变数 < 8 的残基不入样。
- 三分类判定规则：
  - **模型失效**：`e_i` 处于该测定前十分位，且命题信号在检验集上与 `e_i` 的关联成立（CI 不含 0）。
  - **测定本身不可预测**：该测定的残差分布与其官方基准参考上限（其他已发表模型在该测定的 Spearman，仅作解释背景、不作特征）不区分，或命题信号方向与残差关联不成立而信号与测定内部一致性指标关联成立。
  - **无法判断**：以上两条均不成立。

## 4. 一轮命题的格式（B 段每轮落盘件）

```
P<N>: <模型内部信号 X> → <残基级误差 e_i 的方向性预测>
动机：引用 P<N−1> 的具体数值（第 1 轮引用 A 段对照段的具体数值）
检验集：冻结清单中指定测定，P<N> 形成过程中未查看其残基级结果
推翻条件：检验集上 Spearman(X, e_i) 的 95% bootstrap CI（1000 次）覆盖 0
报告：效应量 + CI；三分类判定按第 3 节规则
```

命题文本在看到该轮检验结果之前写入 `results/P<N>_proposition.md`（时间戳留档）。

## 5. 已核实的公开基准口径（A 段对照）

masked-marginals（Meier et al. 2021）：每个突变位点在突变序列上下文中遮蔽该位点，
取 `log p(mut aa) − log p(wt aa)`，多突变求和。官方 ESM2-650M 测定级平均 Spearman ≈ 0.41–0.43
（后经实测为 0.4382，217 测定平均；见 `docs/STATE.md` RSI-2 更正）
（ProteinGym v1.3 口径）。A 段以此为对照，逐测定对比，|Δρ| > 0.05 须说明原因。

## 6. 待你确认的事项

1. 清单按 R2+R3 冻结为 191 个（含 69 个多点测定）——确认 / 修改筛选规则。
2. 误差与三分类定义如上——确认 / 修改。
3. 确认后：下载主数据 zip（1.1GB）与权重（2.6GB），A 段先行，B 段命题逐一落盘。
