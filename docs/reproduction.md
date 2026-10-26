# Reproduction guide

## Level 1 — Offline verification

From the repository root, install the small base package and run:

```bash
python -m pip install -e .
pocketbench verify
python -m unittest discover -s tests -v
```

This needs Python and NumPy only. It checks stored predictions, input hashes, mappings and statistical aggregation. It does not fetch the external reference annotations.

## Level 2 — Independent rescoring

```bash
pocketbench prepare
pocketbench evaluate --published
```

The annotations are fetched from a pinned CryptoBench commit with SHA-256 validation. The labels themselves are not bundled or assigned this repository's MIT license. Keep their upstream source and terms when reusing them. The original-to-canonical and label-to-canonical maps are structural metadata, not a replacement reference annotation.

Expected recovery counts are P2Rank 44/61/61, learned 42/54/62 and learned + PLM 51/65/69 for Top-1/3/5. Rescoring also compares every stored per-target IoU and hit flag, not just these totals.

## Level 3 — Repeat prediction

Linux/WSL is the tested inference environment. Use a separate environment from unrelated research projects. The recorded reference stack is Python 3.11, Lacuna 1.2.0, P2Rank 2.5.1, NumPy 2.4.6, SciPy 1.17.1, Biopython 1.88, PyTorch 2.7.1+cu128, Transformers 4.57.1 and Safetensors 0.8.0. Recorded figure generation used Matplotlib 3.11.2. Java 25 was used for P2Rank.

The original machine had an Intel i9-13900KF and an NVIDIA RTX 5070 with approximately 12 GB VRAM. These are provenance, not minimum hardware requirements. GPU model loading, download, geometry generation and residue embedding timings must not be conflated.

```bash
python -m pip install -e '.[run,plm]'
python -m pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
pocketbench install-p2rank
pocketbench download-model
pocketbench run --p2rank tools/p2rank_2.5.1/prank \
  --model-dir models/esm2_t33_650M_UR50D --out runs/full-100
pocketbench evaluate --run-dir runs/full-100
```

`install-p2rank` verifies the release archive checksum and extracts it into `tools/`. It does not install Java. Python 3.11 with archive extraction filter support is required; use an up-to-date 3.11 patch release or newer.

`download-model` retrieves the five files recorded for the pinned ESM-2 revision and verifies every SHA-256. The approximately 2.6 GB encoder is not committed. An existing verified installation may be supplied instead.

The adapters check the archived Lacuna source/head hashes and P2Rank model hash. They intentionally fail on source drift. If a new upstream package differs, run it as a new experimental version and do not label its outputs an exact reproduction.

PDB inputs are marked as byte-preserved in `.gitattributes`. Do not convert their line endings: the frozen input SHA-256 checks intentionally detect any byte change.

Default PLM inference requires CUDA to match the original FP16 path. `--allow-cpu` permits CPU execution if CUDA is unavailable and records a changed numerical path; it does not force CPU when a GPU is available. Bitwise identical output across hardware and library builds is not promised.

Transformers may report newly initialized pooler weights when loading the encoder. Lacuna uses `last_hidden_state`, not that unused pooler. Its trained residue head is a separate verified upstream file; the workflow does not initialize a new pocket head or train one on the benchmark targets.

The three configurations reuse the same input. The two Lacuna rankers reuse the same complete NMA/alpha candidate lists; pool identity is checked. Neither a failed PLM call nor an absent GPU is silently relabeled as learned + PLM.

## Resume, subsets and failures

Use `--target B003`, repeated if needed, for a subset. A subset has its own selected denominator. Omitting it selects all 100. Output folders have a configuration fingerprint. Matching completed targets are reused; changed inputs, models or runtime metadata require a new output directory. Prediction failures are saved explicitly. Evaluation counts failed or missing selected outputs as unsuccessful and writes `failures.json`.

Do not change the frozen cohort, reference rule or primary threshold after inspecting outcomes. New target selection, new ranker training and a company-tuned P2Rank baseline are separate experiments.

## Figures

```bash
python -m pip install -e '.[plots]'
python scripts/make_figures.py
```

This regenerates English figures from the published results rather than new-run outputs. Figure 1 has Wilson 95% intervals; figures 2 and 9 have pointwise intervals. Figure 10 has paired target bootstrap intervals. Figure 3 is a deterministic target overlap matrix with no invented error bars and a piecewise blue-white-red scale centered at 0.25. Use numeric colorbar values when comparing color intensity across the center.
