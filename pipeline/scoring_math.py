"""Pure scoring math for the A segment (official ProteinGym masked-marginals protocol).

Protocol (verified against proteingym/baselines/esm/compute_fitness.py @144fe22b):
  - mask position i of the WILD-TYPE sequence, forward once, take log-softmax
    at the masked position
  - LLR(mut) = log p(mt aa) - log p(wt aa), both read from the same WT-context row
  - multi-mutants ("A12B:C34D") sum LLRs over their substitutions
  - Spearman is computed per assay over all scored mutants
"""

import numpy as np

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"


def parse_mutant(mutant: str, seq: str, offset: int = 1):
    """Parse 'A12B' or 'A12B:C34D' -> list of (wt, idx0, mt). 1-based idx per ProteinGym."""
    subs = []
    for part in mutant.split(":"):
        wt, idx1, mt = part[0], int(part[1:-1]), part[-1]
        idx0 = idx1 - offset
        if seq[idx0] != wt:
            raise ValueError(f"wt mismatch {part} seq[{idx0}]={seq[idx0]}")
        subs.append((wt, idx0, mt))
    return subs


def features_from_logprobs(row: np.ndarray, wt: str, aa_to_col: dict) -> dict:
    """row: log-probs over the 20 standard aa at one masked position; wt: wt aa there."""
    probs = np.exp(row)
    order = np.sort(probs)[::-1]
    return {
        "logp_wt": float(row[aa_to_col[wt]]),
        "entropy": float(-(probs * row).sum()),
        "margin": float(order[0] - order[1]),
        "top1_prob": float(order[0]),
    }


def llr_at(row: np.ndarray, aa_to_col: dict, wt: str, mt: str) -> float:
    return float(row[aa_to_col[mt]] - row[aa_to_col[wt]])


def zscore(v: np.ndarray) -> np.ndarray:
    sd = v.std()
    if sd == 0 or np.isnan(sd):
        return np.zeros_like(v, dtype=float)
    return (v - v.mean()) / sd


def residue_error_table(positions: list, singles_residual: dict, min_mut=8):
    """positions: iterable of idx0 with mutations; singles_residual: idx0 -> list of |r|."""
    rows = []
    for p in sorted(positions):
        rs = singles_residual.get(p, [])
        rows.append(
            {
                "position": p,
                "n_single": len(rs),
                "e_i": float(np.mean(rs)) if len(rs) else np.nan,
            }
        )
    return rows
