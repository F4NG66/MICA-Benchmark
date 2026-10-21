# Six original MICA metrics on 16 common targets

Mean ± sample SD (ddof = 1), calculated independently across the same 16 targets for every method. Higher is better for all six metrics.

| Method | n | TM-score | Cα match (%) | Cα quality (%) | Aligned Cα length | Sequence identity (%) | Sequence match (%) |
|---|---|---|---|---|---|---|---|
| MICA + AF3 | 16 | 0.9502 ± 0.0459 | 92.98 ± 6.38 | 95.89 ± 2.96 | 588.5 ± 282.4 | 97.94 ± 2.19 | 95.68 ± 3.52 |
| ModelAngelo | 16 | 0.6714 ± 0.3075 | 96.73 ± 5.20 | 76.24 ± 24.80 | 393.1 ± 227.5 | 97.86 ± 3.12 | 97.26 ± 3.14 |
| EModelX + AF3 | 16 | 0.8948 ± 0.1045 | 90.72 ± 9.55 | 91.08 ± 8.09 | 560.4 ± 268.4 | 96.41 ± 4.35 | 92.01 ± 8.40 |
| CryoAtom | 16 | 0.7513 ± 0.2411 | 92.25 ± 8.75 | 89.63 ± 11.49 | 444.1 ± 224.4 | 98.01 ± 3.38 | 96.86 ± 3.79 |
| EMProt + AF3 | 16 | 0.9619 ± 0.0376 | 81.34 ± 14.09 | 95.55 ± 3.58 | 597.1 ± 291.7 | 97.74 ± 3.33 | 96.11 ± 4.11 |
| CryoAtom2 v2.0 | 16 | 0.8346 ± 0.2001 | 89.94 ± 13.16 | 96.19 ± 4.51 | 500.9 ± 233.6 | 97.67 ± 3.19 | 97.23 ± 2.80 |

CA match/quality and sequence match use Phenix chain_comparison at 3 Å in the deposited coordinate frame, with no fitting; these are not the added one-to-one registered precision/recall/F1. TM-score and aligned length/identity use sequence-independent USalign assembly alignment, normalized to all reference Cα residues for TM-score. See docs/metrics.md.

AF3 input audit: MICA used only entity monomers in 12 targets; four targets used entity homomers (9Y6N trimer; 9OKP two separate dimers; 31CK and 9HXD dimers). EMProt/EModelX used entity monomers throughout. This is a configuration difference; the table is not an all-monomer comparison.
