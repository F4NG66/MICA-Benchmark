# Modeling recipe

## Inputs

The archives retain both original EMDB maps and the exact canonical maps used for inference. Canonicalization permutes voxel axes and translates the origin without interpolation; `map_transform.json` records the inverse translation. Decompress `canonical.map.gz` to a fresh `input.map`. Use the supplied entity sequences, deposited copy counts, author contour level, and reported resolution from `manifest/targets.json`. Final scored predictions have already been restored to the original deposited frame.

Prepare assisted inputs in a fresh directory (requires Gemmi in the orchestration environment):

```sh
python scripts/prepare_inputs.py --data ../benchmark_data --pdb 9Y6N --mode mica --output ../run/9Y6N_mica
python scripts/prepare_inputs.py --data ../benchmark_data --pdb 9Y6N --mode monomers --output ../run/9Y6N_monomers
```

MICA output is native per-entity CIF layout with entity/copy FASTA headers. The monomer adapter writes one `ranked_0.pdb` per entity, plus a full-copy FASTA for EMProt. References are never read by this preparation.

## MICA

Clone the pinned revision from `manifest/methods.json` into `WORK/repos/MICA_merizo_fixed`. Acquire official checkpoints and build Pulchra according to upstream README. Apply the one exact patch `patches/merizo_terminal_edge_cases.patch` with `git apply`; verify the before/after hashes in `patches/merizo_provenance.json`. This fixes absent terminal/internal neighbor cases in Merizo without inventing domains or deleting short chains. Original source attribution/license is retained under `third_party/MICA`.

```sh
python scripts/setup_mica_workroot.py /absolute/WORK
export MICA_WORK_ROOT=/absolute/WORK
export MICA_PHENIX_PREFIX=/absolute/phenix-2.2.1-6174
```

Run the portable wrappers in the MICA environment with the repository on `PYTHONPATH`. The required Phenix prefix is explicit; the wrapper retains the native docking process rule and per-domain outcomes. CPU domain processing and GPU model inference are separate:

```sh
cd /absolute/WORK/repos/MICA_merizo_fixed
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
python /absolute/WORK/scripts/audited_mica_af3_processing_fixed_unbounded.py -f /absolute/run/9Y6N_mica/9Y6N.fasta -a /absolute/run/9Y6N_mica/AF3_results -d cpu
python /absolute/WORK/scripts/audited_mica_docking_unbounded.py --density_map /absolute/input.map --contour_level AUTHOR_CONTOUR --resolution RESOLUTION --fasta_file /absolute/run/9Y6N_mica/9Y6N.fasta --AF3_results_path /absolute/run/9Y6N_mica/AF3_results --phenix_act "$MICA_PHENIX_PREFIX/phenix_env.sh" --per_domain_timeout 0
python run.py -m /absolute/input.map -f /absolute/run/9Y6N_mica/9Y6N.fasta -a /absolute/run/9Y6N_mica/AF3_results -o /absolute/run/model --run_pulchra --pulchra_path modules/pulchra304/src/pulchra --resolution RESOLUTION --device cuda
```

Replace AUTHOR_CONTOUR and RESOLUTION with target metadata. In a fresh refinement directory activate Phenix and run `phenix.real_space_refine MODEL input.map historical/config/mica_phenix_2_2_compat.eff resolution=RESOLUTION dry_run=False output.overwrite=False`. Prefer the same-stem refined mmCIF deterministically after validating atom identity against the refined PDB (`select_refined_export.py`); never choose a stage using reference accuracy. Restore the original map frame before scoring with `restore_formal_prediction_frame.py`.

The docking wrapper uses Phenix's Python 3.11 layout from 2.2.1-6174. Its portable edits change location configuration only; environment-specific paths must exist. Historical sources show the original runtime settings. A 20 GiB host reserve guard was used on the old host; later wrappers honor effective cgroup headroom. These operational differences and CPU count affect runtime; this is not a controlled timing experiment.

## Other methods (separate environments)

Use pinned checkouts/checkpoints from `manifest/methods.json`, one job at a time on an assigned GPU. The benchmark commands were:

```sh
model_angelo build -v input.map -pf unique_sequences.fasta -o OUT --device 0
cryoatom build -v input.map -s unique_sequences.fasta -o OUT -d cuda:0 -k
# CryoAtom2 v2.0 environment, not the CryoAtom v1 environment:
cryoatom build -v input.map -ps unique_sequences.fasta -o OUT -d cuda:0 -k
emprot build --map input.map --seq full_copy_sequences.fasta --resolution RESOLUTION --output OUT --device 0 --nt 4 --chain ENTITY1_MONOMER.pdb ENTITY2_MONOMER.pdb
python EMODELX_REPO/run.py --protocol=temp_flex --EM_map=input.map --fasta=ENTITY_COPY_HEADER.fasta --template_dir=ENTITY_MONOMER_FOLDERS --output_dir=OUT --run_pulchra --pulchra_path=PULCHRA --mul_proc_num=4
```

For a one-entity target EMProt has one chain model argument; add one per entity otherwise. EModelX uses adapter FASTA copy labels and the folder structure produced by `prepare_inputs.py --mode monomers`. EModelX is followed by Phenix real-space refinement with the same compatibility parameters and deterministic export selection. Other methods retain their native final-stage exports. Restore any canonical-frame predictions to the deposited frame before comparison.

## AF3 inference itself

The exact processed JSONs allow native AF3 inference without rerunning MSA search. Use the pinned AF3 environment/revision and your separately licensed model parameters:

```sh
python AF3_REPO/run_alphafold.py --json_path=JOB_processed.json --output_dir=FRESH_OUTPUT --model_dir=YOUR_WEIGHTS --db_dir=YOUR_DATABASES --max_template_date=2025-12-31 --run_data_pipeline=false --run_inference=true --num_diffusion_samples=5 --compress_large_output_files=false
```

The original jobs used the intrinsic top-ranked native export, seed 1, and five samples; never rank candidates using the experimental reference. Hardware/library changes can alter numerical inference outputs even with the same seed. The archived selected CIFs reproduce the actual downstream priors exactly. Repeating inference is optional for rescoring or auditing this release.

## Verification status

Summary regeneration, file integrity and AF3 chain/sequence contracts are checked in this release. Fresh full modeling is not rerun; the source pins, adapters, compatibility changes and measured final outputs are supplied so an independent researcher can repeat it. `docs/verification.md` records rescoring checks separately.
