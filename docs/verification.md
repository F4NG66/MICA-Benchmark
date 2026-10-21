# Verification

- All 284 data files passed SHA-256 checks against the server exports and again after packaging.
- All 30 AF3 jobs passed sequence and chain-count checks. All 25 MICA priors match the prepared input hashes.
- Input preparation ran successfully in both modes for all 16 targets, without reading reference coordinates.
- Rescoring reproduced all six metrics exactly for 95/96 rows. CryoAtom/9Y6N TM-score was 0.7054 versus archived 0.7057 across same-source USalign builds; its other five metrics were identical. The tables retain the archived values. [Full differences](../results/rescoring_verification.json).
- Summary tests check means, sample SD, duplicate targets and missing measurements. The archive passed a clean-unpack test.

The old USalign binary required newer GLIBC/GLIBCXX libraries than the second host. Rescoring used a compatible build from the same pinned source. Its 0.0003 TM difference was repeatable; the numerical cause remains unresolved.

Model inference was not rerun during packaging. Commands, source pins, patches and final structures are supplied for reproduction.
