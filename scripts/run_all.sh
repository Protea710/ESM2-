#!/bin/bash
# One command from raw data to all three-round conclusions.
# Usage (on the lab machine): bash /work/run_all.sh
set -euo pipefail
cd /work

echo "== [0/7] environment lock =="
python3 - <<'EOF'
import sys, torch, transformers, pandas, scipy, numpy
print("python", sys.version.split()[0], "torch", torch.__version__,
      "cuda", torch.cuda.is_available(), torch.cuda.get_device_name(0))
EOF
pip3 freeze > /work/requirements-lock.txt

echo "== [1/7] raw-channel checksums =="
sha256sum -c /work/checksums.txt

echo "== [2/7] archive first-run outputs =="
mkdir -p /work/run1_archive
for f in results_a/assay_summary.csv results_a/compare.csv \
         results_p1/P1_per_assay.csv results_p1/P1_summary.json \
         results_p2/P2_per_assay.csv results_p2/P2_summary.json \
         results_p3/P3_per_assay.csv results_p3/P3_summary.json; do
  [ -f "/work/$f" ] && cp "/work/$f" /work/run1_archive/
done

echo "== [3/7] A segment (fresh) =="
rm -rf /work/results_a
python3 -m pipeline.run_a_segment \
  --ref-csv /work/reference/DMS_substitutions.csv \
  --list-csv /work/reference/assay_list_proposed.csv \
  --data-dir /work/DMS_ProteinGym_substitutions/DMS_ProteinGym_substitutions \
  --out-dir /work/results_a
python3 -m pipeline.compare_official \
  --summary /work/results_a/assay_summary.csv \
  --official /work/reference/DMS_substitutions_Spearman_DMS_level_esm2_650M.csv \
  --out /work/results_a/compare.csv

echo "== [4/7] P1 (fresh stats over rerun scores) =="
rm -rf /work/results_p1
python3 -m pipeline.p1_test

echo "== [5/7] P2 (fresh R_i + stats) =="
rm -rf /work/results_p2
python3 -m pipeline.p2_test

echo "== [6/7] P3 (fresh R_i + stats) =="
rm -rf /work/results_p3
python3 -m pipeline.p3_test

echo "== [7/7] master tables + rerun comparison =="
rm -rf /work/master_tables /work/ri_all
python3 -m pipeline.build_master_tables
python3 /work/scripts_lab/verify_rerun.py
echo "RUN_ALL_DONE"
