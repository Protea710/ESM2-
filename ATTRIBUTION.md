# 第三方数据与模型归属

本仓库包含少量第三方文件（均在 `data/reference/`），用于离线校验与文档引用。

## ProteinGym（基准数据与参考文件）

- 来源：OATML-Markslab/ProteinGym，commit `144fe22b07dfaeec2b366f2346203a9838a55b4c`（v1.3）
- 许可：MIT
- 本仓库包含的存在副本：
  - `data/reference/DMS_substitutions.csv`（217 测定清单参考文件，含 GitHub API 标明的 MIT 许可）
  - `data/reference/DMS_substitutions_Spearman_DMS_level_esm2_650M.csv`（官方逐测定 Spearman，用作 A 段对照）
- 未包含：`DMS_ProteinGym_substitutions.zip`（43 MB 逐测定原始 CSV）。下载地址与 SHA256 见
  `docs/STATE.md`，`scripts/run_all.sh` 会校验后使用。

## ESM-2（模型配置与词表）

- 来源：`facebook/esm2_t33_650M_UR50D`，经 ModelScope 同步
- 许可：MIT（模型卡声明；fair-esm 仓库 LICENSE 同为 MIT）
- 本仓库包含的存在副本：`data/reference/ms_config.json`、`data/reference/ms_vocab.txt`
  （仅配置与 33-token 词表，**不含权重**）
- 权重不入库；`docs/STATE.md` 记录了权重文件的 SHA256（与外部发布值一致）

## 引用

若引用本仓库的结果，请同时引用上述数据与模型来源：

- Notin et al., *ProteinGym: Large-Scale Benchmarks for Protein Fitness Prediction and Design*, NeurIPS 2023 Datasets & Benchmarks
- Lin et al., *Evolutionary-scale prediction of atomic-level protein structure with a language model*, Science 379, 1123–1130 (2023)
