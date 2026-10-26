# Target selection and leakage scope

The fixed cohort contains 100 distinct UniProt targets selected from the CryptoBench designated test fold. No target was replaced after its benchmark prediction was inspected.

## Eligibility and diversity

- One annotated author chain and one distinct UniProt per target.
- X-ray resolution <=3.0 Angstrom.
- Entity length and observed standard C-alpha residue count between 80 and 600.
- At least five annotated reference residues; all map to the observed input.
- Summed alternate C-alpha occupancy >=0.8 at every reference residue.
- No identical UniProt in the checked Lacuna training folds or P2Rank CHEN11/JOINED sets.
- No qualifying close-homology alignment to the fetched training/development polymer entities.
- No retained pair with sequence identity >=40% and aligned coverage >=80% on both targets.

The training-screen alignment used Smith-Waterman, BLOSUM62, gap-open 10 and extension 1. A match was excluded at >=40% identity over alignment columns and >=80% coverage of the candidate sequence. This screening covered 1,726 fetched training/development polymer-entity sequences. It reduces detectable close overlap under the stated rule; it does not eliminate all remote homology.

All eligible targets in the smaller functional categories were retained. Other enzymes were sampled by seeded round-robin over the first two EC fields, with seed 20261023. This is a diversity-oriented filtered cohort, not a representative random sample of protein families.

| Category | Targets |
|---|---:|
| Kinase | 2 |
| Protease | 5 |
| Other enzyme | 61 |
| Nucleotide binding / switch | 3 |
| Nonenzymatic binding / regulator / transport | 29 |

The two kinase-like cases are a histidine kinase and an MLKL pseudokinase. Isolated noncatalytic domains were assigned by the represented domain rather than the catalytic activity of the full parent protein. Class counts are not sufficient to establish class-wide model superiority.

## What was checked

The documented Lacuna supervised training lists contain 885 PDB entries across train folds 0–3. The checked P2Rank lists contain 241 unique CHEN11 PDB entries and 347 unique JOINED PDB entries. The retained 100 have no direct PDB overlap with those lists and passed the stated identity and homology screening. The counts refer to unique PDB identifiers, not all dataset lines or polymer entities.

The original cohort hash is `757b6ac4c967a54423a8844586dea09de71f4be105490b7bc243b13ce68e0af3`. A path-free cohort manifest, original training-set identifiers, and the final selection recheck are included. These files preserve the chosen cohort; this repository does not claim to automatically rederive every earlier eligibility decision from live external databases.

## What was not established

- Absence from ESM-2's large-scale sequence pretraining.
- Full historical absence from every calibration or development example of every tool.
- A completely untouched test distribution: Lacuna's authors report observing test performance during iterative development.
- Development-blind external validation or an exact reproduction of the authors' headline benchmark.

Appropriate wording is **a filtered designated-test cohort screened against documented supervised/development sets at specified cutoffs**. Do not describe it as 100 proteins unseen by every learning and development process.
