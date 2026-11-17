"""Torch-side masked-marginal engine (ESM2-650M, HF weights from ModelScope).

For one assay sequence of length L we run ceil(L/B) batched forwards, each batch
containing B copies of the tokenized WT sequence with a different position masked.
Output: log-softmax rows for the 20 standard aa at every masked residue position.
"""

import numpy as np
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

from .scoring_math import STANDARD_AA


class MaskedMarginalEngine:
    def __init__(self, model_path: str, device: str = "cuda", max_tokens: int = 32768):
        self.tok = AutoTokenizer.from_pretrained(model_path)
        self.model = AutoModelForMaskedLM.from_pretrained(
            model_path, torch_dtype=torch.float32
        )
        self.model.eval().to(device)
        self.device = device
        self.max_tokens = max_tokens
        vocab = self.tok.get_vocab()
        self.bos = vocab["<cls>"]
        self.eos = vocab["<eos>"]
        self.mask = vocab["<mask>"]
        self.aa_to_id = {aa: vocab[aa] for aa in STANDARD_AA}
        self.aa_to_col = {aa: j for j, aa in enumerate(STANDARD_AA)}
        missing = [aa for aa in "BXZU" if aa in vocab]
        self.other_id = {aa: vocab[aa] for aa in missing}
        # sanity: BOS/EOS/mask must exist and be distinct
        assert len({self.bos, self.eos, self.mask}) == 3

    def encode_wt(self, seq: str) -> torch.Tensor:
        ids = [self.bos] + [self.tok.convert_tokens_to_ids(c) for c in seq] + [self.eos]
        return torch.tensor([ids], dtype=torch.long)

    def position_logits(self, seq: str, positions: list) -> np.ndarray:
        """positions: 0-based indices into seq (token index = pos+1). Returns (len(positions), 20)."""
        wt = self.encode_wt(seq)[0]
        L = wt.shape[0]
        out = np.empty((len(positions), len(STANDARD_AA)), dtype=np.float64)
        B = max(1, self.max_tokens // L)
        with torch.no_grad():
            for s in range(0, len(positions), B):
                chunk = positions[s : s + B]
                batch = wt.repeat(len(chunk), 1)
                for r, p in enumerate(chunk):
                    batch[r, p + 1] = self.mask
                logits = self.model(batch.to(self.device)).logits  # (B, L, vocab)
                logprobs = torch.log_softmax(logits.float(), dim=-1)
                cols = [self.aa_to_id[a] for a in STANDARD_AA]
                rows = logprobs[
                    torch.arange(len(chunk)), torch.tensor([p + 1 for p in chunk])
                ][:, cols]
                out[s : s + len(chunk)] = rows.cpu().numpy().astype(np.float64)
        return out
