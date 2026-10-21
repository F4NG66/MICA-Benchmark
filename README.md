# Cryo-EM modeling benchmark

Six methods on 16 common EMDB/PDB targets released in 2026: MICA, ModelAngelo, EModelX, CryoAtom, EMProt and CryoAtom2 v2.0.

The full cohort has 25 targets. These 16 have valid outputs from all methods after excluding 10OY for exact training overlap. Results and sample SD are in [the summary table](results/common16_table.md); individual scores are in [per_target_metrics.json](results/per_target_metrics.json).

## Run

Python 3.11 or later:

```sh
python scripts/summarize.py
python -m unittest discover -s tests
```

Extract both data archives into the same directory, then verify:

```sh
python scripts/verify_data.py --data ../benchmark_data
```

To repeat scoring, install `requirements-evaluation.txt`, USalign and Phenix 2.2.1-6174:

```sh
python scripts/rescore.py --data ../benchmark_data --output ../rescored \
  --usalign /path/to/USalign --phenix /path/to/phenix_env.sh
```

See [modeling commands](docs/modeling.md), [metric definitions](docs/metrics.md) and [version pins](manifest/methods.json).

## AF3 inputs

MICA used entity monomers for 12 targets and homomers for four: 9Y6N, 9OKP, 31CK and 9HXD. EMProt and EModelX used entity monomers throughout. [Input audit](docs/af3_input_audit.md).

Data archives: `common16_inputs_and_models.tar.gz` and `common16_af3_inputs.tar.gz`. They contain maps, references, 96 final models and 30 AF3 input/output pairs. [Data terms](docs/data_licenses.md).
