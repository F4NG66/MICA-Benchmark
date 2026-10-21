# AF3 input audit

Original processed JSONs, selected CIFs and MICA adapter manifests were checked. All 30 AF3 jobs match their full input sequences and chain counts; all 25 entity priors match the SHA-256 of the CIF actually supplied to MICA.

| Target | MICA AF3 input |
|---|---|
| 9Y6N | One trimer |
| 9OKP | Two separate homodimers |
| 31CK | One dimer |
| 9HXD | One dimer |
| Other 12 targets | Separate entity monomers |

EMProt and EModelX used entity monomers for every target. MICA then extracted chains, segmented domains with Merizo and docked them into density.

The paper describes prediction for each unique chain. Its FASTA converter also supports repeated entity copies from header chain labels. This benchmark followed that copy-aware code path. An all-monomer comparison would require rerunning the four affected MICA cases.

Evidence: [chain counts and hashes](../manifest/af3_chain_audit.json), the original JSON/CIF files in the data archives, and [the pinned converter](https://github.com/jianlin-cheng/MICA/blob/726d1bb14be04887a43b7f5aefcff0d843761d01/utils/fasta_to_AF3_json.py).
