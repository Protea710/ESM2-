import hashlib
import json
import sys
from pathlib import Path

EXPECT = {
    "DMS_substitutions.csv": "a8f498011532a74aa9fe556a50555a75e928c5837d19c06a87592ae04049b308",
    "DMS_substitutions_Spearman_DMS_level_esm2_650M.csv": "f432423b87f79ac9778dfac86e3d95be041d246bc618dc0a406b35b0b7466437",
    "ms_config.json": "539095c22efc52a09d6147074ba4ca119f76a890df5901213b2b55f7d2f96b2b",
    "ms_vocab.txt": "0b82cc0a7c7cf9e567b1e5892d793285b9fbae822c964ca48696f7db44598e03",
}
REF = Path("data/reference")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


failures = []
for name, expected in EXPECT.items():
    p = REF / name
    if not p.exists():
        failures.append(f"missing {name}")
        continue
    got = sha256(p)
    if got != expected:
        failures.append(f"mismatch {name}: got {got} expected {expected}")
    else:
        print(f"OK {name} {got}")

if failures:
    print("\n".join(failures))
    sys.exit(1)
print("channel verification: all 4 checksums match")
