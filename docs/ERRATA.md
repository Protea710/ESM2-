# ERRATA — 已发布文本中的数字更正

本文件记录已发现并更正的事实错误。冻结件（`*_PROPOSITION.md`、`*_FORMATION.csv`、
`*_TEST.csv`）哈希封印，一律不改；其内部错误在下方逐条指出，以本文件为准。

## E1（2026-11-17 发现）P1 方向计数：50/61 → 48/61

- **错在哪**：此前在 RSI-4 的回复、`docs/FINAL_REPORT.md` §1、`results/p1/P1_RESULTS.md`、
  `results/p2/P2_RESULTS.md`、`docs/STATE.md` 中写「P1 的 Δ_a < 0 的测定 50/61；Δ_a > 0 仅 11 个」。
- **正确值**：以 `results/p1/P1_per_assay.csv` 逐行实测——61 个合格测定中
  **Δ_a < 0 者 48 个，Δ_a > 0 者 13 个，无零值**（48 + 13 = 61）。
- **原因**：撰写时手抄计数错误，未回表核对。θ = −0.1390 与 CI [−0.1767, −0.0839] 不受影响
  （由同一张表计算，且重跑逐位复现），判定 REFUTED 不变。
- **处理**：`FINAL_REPORT.md`、`P1_RESULTS.md`、`P2_RESULTS.md`、`STATE.md` 已改为 48/61 并标注本条目；
  `figures/build_notes.py` 复核页与 `figures/NOTES_A4.pdf` 已重建，文字同步为「已按 ERRATA E1 更正」；
  `FIG_PLATE_A2.pdf` 自始使用表值 48/61，无需改动。
- **未改动的冻结件**：`results/p2/P2_PROPOSITION.md`（第 10 行）与 `results/p3/P3_PROPOSITION.md`
  （第 8 行）在动机部分引用了 50/61。二者哈希封印，**按要求保持原样**；其引用的正确值为 48/61，
  且该引用仅作动机陈述，不进入任何统计量或判定。

## E2（2026-11-17 发现）复核页文字与更正后的正文不同步

- `figures/NOTES_A4.pdf` 第 2 页原写「FINAL_REPORT §1 与 P1_RESULTS.md 写 50/61……两份文件的文字未改动」。
  E1 更正后该句失效，已改为「此前误写 50/61、11 正，已按 ERRATA E1 更正为 48/61」，PDF 已重建。
