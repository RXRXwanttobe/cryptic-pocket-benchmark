# Limitations and appropriate use

This is an independent software comparison and learning project, not a newly trained pocket predictor or a peer-reviewed validation study.

The strongest supported result is improved Top-5 ranking by the PLM configuration relative to learned on the fixed 100 targets. Its observed lead over default P2Rank is uncertain. The sample selection is filtered and functionally imbalanced, and statistical intervals do not account for every source of selection bias or annotation uncertainty.

The reference residue union and residue-only criterion differ from the upstream benchmark's selected-site/distance criteria. Comparisons of absolute headline percentages across these protocols are not meaningful.

All configurations receive the same prepared single chain. That excludes cofactors, ions, partner chains, biological assembly changes, and ligand-driven effects. NMA samples a local elastic-network approximation and can miss large domain rearrangements. Side-chain coordinates are propagated rather than obtained from an atomistic MD simulation.

ESM-2 sequence familiarity is not equivalent to supervised pocket-label leakage, but its pretraining nonoverlap was not established. The documented training-list checks do not make the cohort an untouched external test.

No Pharmaco-net implementation, internal compound collection, or company-specific target information is included. Default P2Rank performance is not evidence of superiority over a separately tuned production model. A future integration should compare the actual tuned baseline and candidate fusion on a separate development set, with equal final candidate budgets and independent evaluation.

The figure uncertainty refers to target recovery estimates. It is not molecular binding uncertainty, an MD sampling error, or calibrated confidence in a ligand pose.
