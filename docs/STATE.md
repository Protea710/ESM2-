# STATE

## 环境（lab 机实测 2026-11-16）
- GPU: RTX 4090 24G，驱动 595.91.07，CUDA 13.2
- torch 2.11.0+cu128，cuda available；Python 3.10.12；128 核；503G 内存；256G 盘（空）

## 通道校验值（SHA256，2026-11-16 实算）
- data/reference/DMS_substitutions.csv
  a8f498011532a74aa9fe556a50555a75e928c5837d19c06a87592ae04049b308
  来源：OATML-Markslab/ProteinGym @144fe22b (v1.3) reference_files
- data/reference/DMS_substitutions_Spearman_DMS_level_esm2_650M.csv
  f432423b87f79ac9778dfac86e3d95be041d246bc618dc0a406b35b0b7466437
  来源：同 commit benchmarks/DMS_zero_shot/substitutions/Spearman（A 段对照对象）
- data/reference/ms_config.json
  539095c22efc52a09d6147074ba4ca119f76a890df5901213b2b55f7d2f96b2b
- data/reference/ms_vocab.txt
  0b82cc0a7c7cf9e567b1e5892d793285b9fbae822c964ca48696f7db44598e03
  来源：modelscope.cn facebook/esm2_t33_650M_UR50D master（权重通道已通，2.6GB 权重体待确认后下载）

## 清单统计（scripts/build_assay_list.py 实算）
- 217 测定中：seq_len≤1022 → 201；总突变数≥500 → 202；交集 → 191
- 191 测定合计 2,411,093 突变；masked-marginals 前向上界 1,822,479 次
- 含多点突变测定 69 个（B 段非加性轮的检验集来源）
- 详细清单 data/assay_list_proposed.csv

## 更正（RSI-2）
- GATE1 第 2 节「masked-marginals 前向上界 1,822,479 次」作废：官方口径是**野生型序列逐位置**遮蔽，
  每测定 L 次前向，191 测定合计 **54,182 次**（build_assay_list 的 n_unique_positions_upper 不再用于预算）。
- GATE1 第 5 节「官方 ESM2-650M 测定级平均 Spearman ≈ 0.41–0.43」修正为实读值
  **0.4382**（217 测定平均，官方 CSV f432423b…，2026-11-16 lab 机实读）。

## A 段通道补充校验值（RSI-2 实算）
- /work/esm2_model/model.safetensors（ModelScope 下载，2.6GB）
  a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0
  = AtomGit 昇腾适配文档独立发布的同一 SHA256 → 权重通道交叉验证通过
- /work/DMS_ProteinGym_substitutions.zip（marks.hms.harvard.edu v1.3，43MB）
  3a83766254ac9ac9984ec25cb73c6e010ea4418f5e35f143933e6b6e6473b921
- 冒烟验证：NEIME 测定 our ρ=0.0297 vs 官方 0.030（逐位一致）；单测 3/3 通过

## RSI-1 自评
强：通道、清单、口径全部落在可校验文件上，本轮即可被 gate 证明。
弱：6h 预算按 1.82M 前向上界估，多点突变的 k 次前向未计，取数后需立即精确化；
    A 段若超时按 DMS_index 顺序报未完成项（已写入 GATE1）。
敌意攻击点：三分类判定规则在无重复测定的 assay 上可能退化为「无法判断」——
对策：测定不可预测的判据同时使用官方基准里其他模型的测定级 Spearman 上限作背景（只作解释，不作特征）。

## 当前阶段（RSI-12 末）
最终交付完成（gz 本体搬运通道待用户选定，其余全部落地）：
- 逐突变总表：singles 641,637 / doubles 826,245 / higher 943,211 = 2,411,093 行（逐位核对），
  SHA256 与列清单入 results/master_tables/MANIFEST.json；本体在 lab /work/master_tables/（132MB 实测）。
- 单命令重跑：bash /work/run_all.sh（lab）；requirements-lock.txt 落盘；重跑比对
  compare.csv 与三轮 per_assay 全部 byte_equal，θ/CI 逐位一致；唯一差异 seconds 列（±2.2s）。
- 最终报告 docs/FINAL_REPORT.md：三轮命题/冻结哈希/判定/效应量并列表、两处执行修正
  （CSV 误读、门禁误删）、P4 下一轮命题及理由。
- 本轮实现修正三处（均不影响口径，全部留档于 MANIFEST/report）：
  build_master_tables 的 DMS_id 缺列、空双点表 split KeyError、higher 表循环内 set_index（O(n²) 卡死）。

## P1 检验结果（RSI-4，冻结规则原样执行）
- 冻结完整性：三件冻结物 SHA256 与 P1_FREEZE.json（8a3a11c8…ed38，冻结于 17:43:34Z）逐位一致。
- 检验集 65 测定：61 合格，4 排除（合格位点<8：GCN4_Staller_2018=0、HIS7_Pokusaeva_2019=0、
  CASP7_Roychowdhury_2020=2、Q59976_Romero_2015=5）。
- θ = −0.1390（Δ_a 中位数，测定等权）；95% 簇 bootstrap CI = [−0.1767, −0.0839]（10,000 次，seed 20261117）。
- 判定：下限 −0.177 ≤ 0.10 → **REFUTED**；且点估计为负——margin 低四分位位点的
  残基级平均绝对残差反而**小于**高四分位位点（Δ_a<0 的测定 48/61；ERRATA E1）。
- 全表 results/p1/P1_per_assay.csv；结果文件 results/p1/P1_RESULTS.md（时间戳晚于冻结哈希）。

## 更正（RSI-4）
- pipeline/assay_scoring.py 中 e_i/n_single 曾按 0 基位点索引对 1 基 position 列合并，
  导致 lab 机 /work/results_a/positions/*.csv 的 e_i、n_single 两列错位 1 个位点。
  模型分、rho、逐突变残差不受影响。P1 检验未使用这两列（p1_test.py 从逐突变表重算）。
  已修复并加对齐单测；下次重跑 A 段时 positions 表将以正确对齐重建。
- RSI-3 verdict 为 blocked（buildconfig 模块缺 suite），本轮以 `make verify` suite 覆盖修复。

## A 段结果（RSI-2，lab 机实测）
- 官方口径确认：OATML-Markslab/ProteinGym @144fe22b `proteingym/baselines/esm/compute_fitness.py`
  的 masked-marginals 是在**野生型序列**上逐位置遮蔽取 log-softmax，LLR = log p(mt)−log p(wt)
  全部读自 WT 上下文行，多突变按 `:` 分隔求和；测定级 Spearman 含单点+多点全体突变。
- 运行：191/191 测定，总前向 54,182 次，fp32，batch 按长度自适应，墙钟 **899.4s**（预算 6h 的约 4%）。
- 对照（pipeline/compare_official.py）：
  - mean ρ (ours) = 0.45042 vs official same-assay mean = 0.45046（Δ = −4e-5）
  - max |Δρ| = 0.000506（191 测定全体）；|Δρ|>0.02 的测定数 = 0；>0.05 = 0
  - 逐测定明细 results/compare.csv（本地副本，与 lab 逐位一致）
- 结论：A 段完成，与公开基准一致，流程未跑偏。逐突变级总表与逐位特征表在
  lab 机 /work/results_a/{mutants,positions}（2.41M 行 + 54k 行；B 段特征已同步产出：
  logp_wt、entropy、margin、top1_prob、pos_norm、e_i、n_single、n_multi）。

## RSI-2 自评
强：口径先读官方源码再实现，冒烟（NEIME 0.0297 vs 0.030）先行，全量后逐测定对照全绿；
    单测 3/3、通道 4/4、清单校验通过。
弱：逐突变总表仍在 lab 机，尚未交付到用户可见处；B 段特征列尚未做残基级检验；
    fp32 HF 权重与官方 fair-esm 权重存在 5e-4 级数值差（Spearman 秩意义下无影响，已量化）。
敌意攻击点：检验集若与命题形成集来自同一蛋白家族，残基级相关可能被家族性伪影污染——
对策：B1 检验集指定时跨 coarse_selection_type 分层抽样，且命题形成只用测定级 ρ 与官方对比，不看残基级 e_i。
