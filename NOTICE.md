# Attribution and component scope

The MIT license in this repository applies to the original workflow code, prose and figures contributed here. It does not replace upstream terms or transfer authorship of the predictors, dataset or encoder.

- **Lacuna** is developed by Clayton Moore. The fixed upstream source reviewed for this experiment carries the MIT license and copyright notice. This repository imports the installed package; it does not vendor the original detector, ranker implementation or learned weights.
- **P2Rank** is an upstream predictor distributed separately. The optional installer downloads the official pinned release and preserves its files. Its license remains the release's upstream license; consult the downloaded distribution.
- **ESM-2** is a Meta AI model with its own model card and license. The encoder and Lacuna's trained residue head are not authored here. The large encoder is downloaded separately and is not committed.
- **PDB coordinate data** are available under CC0 according to the wwPDB/RCSB usage policy. Bundled inputs are prepared derivatives of deposited coordinates, with source PDB IDs and reversible residue numbering maps. Retain the source IDs and cite structure authors when appropriate.
- **CryptoBench annotations** remain external. The CLI explicitly fetches a pinned author-provided annotation file into an ignored cache and maps it locally. We do not bundle that raw annotation file or assign it an MIT license. Consult the upstream dataset source and paper for its terms.

No association with or endorsement by the upstream authors, Recursion, or an employer is implied. No Pharmaco-net internal implementation or proprietary compound/target collection is included.
