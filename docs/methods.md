# Methods

## Task and frozen cohort

The experiment asks whether a known annotated cryptic site can be recovered from a single-chain, ligand-free input. The 100 target identities were fixed before their benchmark predictions were inspected. The cohort file preserves the original selection and execution hashes. The later PLM addition followed inspection of the initial two-configuration results, so the three-way study is an incremental comparison rather than a wholly blinded prospective study.

Inputs contain selected protein heavy atoms and one selected alternative conformation. Other chains, waters, ions, ligands, and biological assembly effects are not modeled. All residues were canonically numbered 1..N in chain A for all three configurations. This avoids negative-number and insertion-code problems in the Lacuna PLM path. Maps to author residue numbers and mmCIF label numbering are included. Retained coordinates, amino-acid identities, occupancies and B factors were checked against the original prepared chains.

## Predictors

P2Rank 2.5.1 uses the default model and four threads. Its native pocket residue assignments and ranks are evaluated.

Lacuna 1.2.0 uses the alpha detector and an ANM/NMA ensemble: the input plus 20 generated conformers. NMA samples linear combinations of low-frequency elastic-network modes. The installed backend uses an 8 Angstrom contact cutoff, ten modes and a maximum C-alpha RMS displacement parameter of 2 Angstrom. Those are model parameters, not a measured opening amplitude. The original source hashes are provided.

Both Lacuna rankers operate on the same complete candidate lists. Learned uses 23 geometry/ensemble features. Learned + PLM uses those features plus the mean, maximum, top-three mean and fraction >=0.5 of residue scores from the official ESM-2 linear head. Its geometric coefficients are separately fitted. We used the authors' weights without training on our 100 targets.

The encoder is `facebook/esm2_t33_650M_UR50D`, revision `08e4846e537177426273712802403f7ba8261b6c`. Original inference used CUDA FP16 with float32 residue-head calculation. The weights' SHA-256 and the supplied task head were independently checked. Changing the encoder requires a compatible refitted head; a different language model is not a drop-in replacement.

## Reference definition and scoring

The pinned CryptoBench `label_dataset.json` identifies reference residues with author chain ID and mmCIF label sequence number. We map them to the canonical input numbering. A reference may be the union of multiple annotated cryptic sites; we do not reconstruct a new holo-ligand contact set or choose a main pocket.

For a predicted pocket P and reference R:

```text
IoU(P, R) = |P intersect R| / |P union R|
Top-k hit = any individual pocket of rank <= k has IoU >= 0.25
Recovery = hit targets / selected targets
```

The union of multiple predicted pockets is not used to obtain a Top-k hit. Enlarged pockets are penalized through the union denominator. The complete candidate pool is retained to distinguish detection coverage from ranking conversion.

All 300 original outputs completed. New-run failures and missing selected outputs count as unsuccessful in the fixed selected-target denominator and are listed separately. No replacement targets are chosen. A subset run is reported with its own explicit denominator and must not be called the 100-target result.

## Statistical uncertainty

Marginal recovery proportions use asymmetric Wilson 95% confidence intervals. Curves use pointwise intervals; these are not simultaneous bands. The same target appears at every k and threshold, making the points correlated.

Top-5 PLM-minus-learned and PLM-minus-P2Rank differences use target-paired bootstrap intervals (10,000 resamples; seed 20261023). Exact two-sided McNemar tests use discordant target outcomes. Holm correction is applied to the two primary comparisons. Intervals are per-comparison, not simultaneous family-wise intervals. Top-1, Top-3, category breakdowns, and alternate IoU thresholds are exploratory. Overlap of marginal error bars is not a test of a paired difference.

## Numerical and causal interpretation

NMA samples are not an MD trajectory. Candidate persistence across these samples is not an equilibrium opening probability. The model does not establish biological activity, binding affinity, or clinical benefit.

The learned-to-PLM comparison isolates a changed ranking configuration over an unchanged candidate pool. It does not identify the individual causal contribution of every PLM feature or coefficient. Oracle coverage over all candidates is retrospective, and differing numbers of candidates prevent an equal-budget interpretation between tools.
