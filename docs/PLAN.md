# PLAN — 工作分解树（最终态）

- G1 关卡一页纸 —— done（用户确认后清单冻结 191 个）
  - G1a 通道验证（4 个校验值）—— done
  - G1b 测定清单 + 筛选规则（191/217）—— done
  - G1c 命题格式 + 误差定义 + 三分类判定规则 —— done
- A 段 对照段 —— done
  - A1 权重下载 + 权重指纹（a08adabb…4db0 与外部发布值一致）—— done
  - A2 masked-marginals 打分管线（fp32, 4090, 899s, 191/191）—— done
  - A3 测定级 Spearman vs 官方 CSV（max|Δρ|=0.0005, 0 超阈）—— done
- B 段 主体段 —— done
  - B1 P1 命题与检验 —— done（REFUTED）
  - B2 P2 命题与检验（引用 P1 数值）—— done（INDETERMINATE）
  - B3 P3 命题与检验（多点非加性）—— done（REFUTED）
- D 交付 —— done
  - D1 逐突变级总表（2,411,093 行，3 表 + MANIFEST）—— done
  - D2 重跑单命令 + 锁环境（run_all.sh + requirements-lock.txt）—— done
  - D3 下一轮命题（P4 蛋白内一致性）—— done（未检验）

检验入口：`make verify`（通道+关卡）、`scripts/check_p*.py`（冻结与结果一致性）、
`figures/audit.py`（图面数字逐项回表）、`bash /work/run_all.sh`（lab 机全量重跑）。
