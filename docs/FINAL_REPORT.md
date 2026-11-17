# 最终报告：模型知道自己错在哪吗——零样本突变打分的残基级可预测性

本报告只汇总三轮预登记命题的执行结果；不新增分析，不对结果做事后解释。
所有命题、检验集、阈值均在看到对应检验集统计之前冻结（哈希与时间戳见各 FREEZE.json）。

## 1. 三轮命题总表

| 轮 | 命题（一句） | 冻结哈希（PROPOSITION） | 检验集 | 效应量 | 95% 簇 bootstrap CI | 判定（冻结规则） |
|---|---|---|---|---|---|---|
| P1 | margin 低四分位位点的残基级平均绝对残差比高四分位至少高 0.10（测定内 z 残差尺度）；θ=median Δ_a | `3b8ae8f3536c39a0ddbc87a6132e6444113b9bc16a6399a1cf90c977b17bbf89`（冻结 2026-11-16T17:43:34Z） | 65 测定，61 合格 | θ = **−0.1390** | [−0.1767, −0.0839] | **REFUTED**（下限 ≤0.10；且方向相反，Δ<0 者 48/61——见 ERRATA E1） |
| P2 | LLR 按位点尺度归一（LLR/R_i）后 P1 反向关联消失：\|θ_norm\| ≤ 0.05 | `ef0e814e43008da4befddd02b7b88f628042476f0842ff9c585b94feee3a3865`（冻结 2026-11-16T18:00:42Z） | 48 测定，42 合格 | θ_norm = **+0.0699** | [+0.0012, +0.1289] | **INDETERMINATE**（ci_high>0.05 不成立、ci_low<0.05 不推翻） |
| P3 | 双点突变实验非加性幅度 \|ε_z\| 高档的归一化残差比低档至少高 0.10；θ₃=median Δ₃_a | `618840c25f90209a8e275f05c75391d31eb8d18b559727f5bdf1bab32c46b607`（冻结 2026-11-16T18:17:13Z） | 26 测定，26 合格（130,340 双点） | θ₃ = **+0.1366** | [+0.0022, +0.2904] | **REFUTED**（下限 ≤0.10） |

A 段对照（口径冻结前提）：191/191 测定，masked-marginals（WT 上下文逐位遮蔽，官方
compute_fitness.py @144fe22b 口径），mean ρ = 0.45042 vs 官方同集合 0.45046，
max |Δρ| = 0.000506，无测定超 0.02。种子与冻结划分：P1/P2/P3 split salt
`P1/P2/P3-split-v1-2026-11-17`；bootstrap seed 20261117（10,000 次）。

## 2. 中途执行修正（如实记入）

1. **CSV 误读（RSI-10，P3 首跑）**：`pipeline/p3_test.py` 首跑误读原始 DMS CSV
   （无 model_score 列）而崩溃；修正为读取 A 段产物 `/work/results_a/mutants/<id>.csv`
   后重跑。仅程序数据源路径错误，统计口径不变；修正后全量重跑通过。
2. **门禁误删（RSI-10 重提交）**：提交 rsi_round 时误将 `p3_census_data` 模块从模块图谱
   删除（其文件同时被 refdata 的 glob 覆盖），门禁判 blocked；恢复模块条目后全绿。
   纯清单登记错误，不涉及任何分析或数据。
3. （此前已留档的第三处）RSI-4 发现并修复 `assay_scoring.py` 的 e_i/n_single 位点错位
   合并缺陷（0 基对 1 基列），并加对齐单测；模型分与 rho 不受影响。

## 3. 逐突变级总表

全量建表在 lab 机完成（191 测定、2,411,093 突变，R_i 对全部突变位点重算，
logp(wt) 重算与存表最大偏差 **0.0**）。总表体积实测 **132MB**（singles 31.8MB +
doubles 57.1MB + higher 43.1MB；替换此前 ~80MB 的未验证预估）。
产物与 MANIFEST（行数、列名、SHA256）：`/work/master_tables/`（lab 机），
清单逐位副本 `results/master_tables/MANIFEST.json`（本地可见目录）。

| 文件 | 行数 | SHA256 |
|---|---|---|
| singles.csv.gz | **641,637** | `8b53ed33c91db9e45bda9d905f55055ac91622381a6a21e5e084195a2d5c0b3f` |
| doubles.csv.gz | **826,245** | `dd9eeceedb3c2740d8e2db77306c70724c4c96ce03e9c08e0c009a12ca3b442c` |
| higher.csv.gz | **943,211** | `dc05feaed3d7d8ef242ec88c89c4592b3b67ee8eee927eafed0cccc1aa907109` |

三表行数合计 2,411,093 = 冻结清单总突变数（逐位核对）。
列清单（单点表）：DMS_id, mutant, position, wt_aa, mut_aa, DMS_score, model_score(=LLR),
r, R_i, llr_norm, margin, logp_wt, entropy, top1_prob, pos_norm, e_i, n_single_at_pos,
p1_set, p2_set, p3_set, rounds_test。双点表在两位点各展开 pos/wt/mut/llr/R_i/llr_norm/
margin/entropy，另含 eps（实验非加性，组份单点缺失时为空）；≥3 点表含 k/positions/
wt_list/mut_list/llr_norm_total。各列语义见 MANIFEST.json 的 column_notes。
gz 本体在 lab 机；本地可见目录暂存清单与锁文件（大文件搬运通道见第 5 节）。

## 4. 单命令重跑与环境锁定

- 单命令（lab 机）：`bash /work/run_all.sh`
  流程：环境锁定（pip freeze → `/work/requirements-lock.txt`）→ 原始通道校验
  （`/work/checksums.txt`：DMS zip 3a837662…、model.safetensors a08adabb…、参考 CSV×3）→
  归档首轮结果 → A 段全新重算 → P1/P2/P3 全新统计（R_i 缓存清除重算）→ 总表重建 →
  `verify_rerun.py` 逐位比对。
- 种子：划分 salt `P1/P2/P3-split-v1-2026-11-17`；bootstrap seed `20261117`；
  无任何模型训练/微调。
- 重跑比对结果（verify_rerun.py 实测输出，RSI-12）：
  - `results_a/compare.csv`：**byte_equal=True**（A 段 191 测定逐测定 ρ 与 Δρ 逐字节一致）；
  - P1/P2/P3 per_assay 表：全部 **byte_equal=True**（61/42/26 行逐字节一致）；
  - P1 θ=−0.1390358746、P2 θ_norm=+0.0699415570、P3 θ₃=+0.1365596124 及全部 CI：**逐位一致**；
  - 唯一差异：assay_summary.csv 的墙钟列 seconds（最大差 2.2s，非科学量；rho 列逐字节一致）。
  - 结论：A 段指标与三轮 θ、区间均逐位复现。

## 5. 大文件交付通道说明

总表 gz 本体合计实测 132MB，已通过 exec 落盘通道（单次 ≤10MiB，base64 分片 19 片，
15MiB/7.5MiB 二进制片可被 3 整除避免补位跨片）完整搬入本机可见目录
`results/master_tables/`；三文件 SHA256 与 MANIFEST 逐位一致，本机 gunzip -t
完整性与行数（641,637 / 826,245 / 943,211）复验通过。lab 机 `/work/master_tables/`
保留同一批文件可交叉复核。

## 6. 下一轮命题及值得测的理由

**P4（蛋白内一致性假说）**：同一 UniProt 的多个测定（如 BLAT_ECOLX×4、RL40A×3、SRC×3、
RASK×2、HSP82×3、CP2C9×2、VKOR1×2、KCNE1×2、PTEN×2、S22A1×2、DLG4×2）之间，
残基级误差结构（e_i 的秩剖面）的一致性高于不同蛋白之间；即模型内部信号对误差的预测力
在蛋白内显著强于跨蛋白。

值得测的理由（全部引用已留档数值，不做机制解释）：
1. P1 的跨测定效应量为 −0.139（CI [−0.177, −0.084]），但 Δ_a 的测定间散布极大
   （−0.513 至 +0.721，P1_per_assay.csv）——跨测定聚合可能稀释了蛋白内的一致结构。
2. P3 的方向计数为 18 正/8 负，但两个最大测定（GRB2 −0.094、PABP −0.075）为负且占
   双点 76%——簇级异质性是三轮 CI 偏宽的共同来源。
3. P2 的残余正关联（θ_norm=+0.070，CI [+0.001, +0.129]）未达判定阈值，其上限
   提示存在量级 ≤0.13 的真实效应空间；按蛋白分层可能把效应从噪声中分离出来。
4. 该命题不改变误差定义与归一化口径（沿用 P2 口径），仅需把分层单元从测定改为 UniProt，
   检验集仍取从未做过残差统计的洁净测定——设计上是三轮结果的直接下一步。
