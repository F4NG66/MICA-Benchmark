# Metrics

Each target has equal weight. Summary = arithmetic mean ± sample standard deviation across 16 targets. Fractions are converted to percent only for display. No missing result is imputed as zero.

| Metric | Implementation | Interpretation |
|---|---|---|
| TM-score | USalign `-mm 1 -ter 1 -split 2 -mol prot -outfmt 2`, TM2 normalized by all reference protein Cα residues | Global assembly fold agreement after a single sequence-independent rigid alignment |
| Cα match | Phenix matched Cα / predicted Cα, maximum distance 3 Å | Fraction of predicted Cα near the reference in the deposited frame |
| Cα quality | Cα match × predicted Cα / reference Cα | Reference-normalized geometric coverage using the Phenix count |
| Aligned Cα length | USalign Lali | Number of structurally aligned residues; target-size dependent |
| Sequence identity | USalign IDali | Residue-type identity within the structural alignment |
| Sequence match | Phenix `Seq match (%)` / 100 | Sequence agreement reported among spatially matched residues |

Phenix uses full target/query, `unique_query_only=False`, `unique_part_of_target_only=False`, `test_unique_part_of_target_only=False`, `use_crystal_symmetry=False`, `max_dist=3.0`. Cα match and quality are the original-style Phenix metrics; they do not enforce the evaluator's added one-to-one residue assignment and must not be relabeled precision/recall/F1. Their values can differ from strict one-to-one scores, particularly for overlapping predictions.

The evaluator additionally reports one-to-one geometric, amino-acid-type and construct-position-registered precision/recall/F1 at 3 Å, without fitting. These complement the six metrics: good internal folds can be incorrectly placed in density; aligned-only identity can remain high for an incomplete model. Registered precision uses predictions mapped to reference-observed construct positions; unobserved/unregistered predictions are reported separately. This score does not measure peptide connectivity.

USalign excludes prediction chains shorter than three residues from its structural alignment; their number and the actual prediction denominator are recorded. Primary reference TM-score requires the complete reference denominator. Raw USalign and Phenix logs are saved by each fresh scoring run.
