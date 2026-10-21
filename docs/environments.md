# Environments and version pins

Use separate environments per method; do not modify a shared AF3 installation. Exact method revisions and source links are in `manifest/methods.json`. Checkpoints must match the method release, especially original CryoAtom2 v2.0, not later v2.1 weights. Large upstream weights are downloaded separately; this release does not redistribute them.

Historical MICA environment: Python 3.8 family, PyTorch 1.8.1+cu111, NumPy 1.19.1, SciPy 1.5.2, Biopython 1.78 and Open3D 0.13.0. Upstream environment files are supplied as source evidence; older CUDA wheels may need an appropriate container/driver. AF3 was run separately with native source revision ebfe70a27a6a1ad18c77664191c5f7fa486ebee9. Do not put AF3 dependencies in the legacy MICA environment.

Phenix 2.2.1-6174 was installed separately. Its compatibility parameter file replaces obsolete parameters from the original release; it is not an accuracy-tuned refinement configuration. Phenix docking and real_space_refine are CPU workloads in this protocol. Obtain Phenix under its own download terms.

The evaluation requirements are separate from all inference environments. USalign used source commit 1fa25a958fe3095900eba2ef9b48562bdc462251 (commit date 2026-09-21; program version 2026-09-20); the server-specific binary was rebuilt for its GLIBC. Install from https://zhanggroup.org/US-align/ and retain its version/source hash. Use the distributed final structures to rescore without requiring a GPU, AF3 databases or inference weights.

The benchmark used an A100 80 GB host (~108 GiB RAM, 28 logical CPUs) and later a MICA host (~1 TiB RAM, 128 logical CPUs, RTX 4090). MICA native docking selects min(int(0.6*visible_CPU),48) processes; changing visible CPUs may also change its spatial search subdivision. Do not treat wall times on these two hosts as controlled same-hardware comparisons.

The exact source archive and its upstream SHA-256 are under `third_party/USalign_source*`. Extract it and build `USalign` using its Makefile on the destination host. See `docs/verification.md` for the observed cross-build TM difference of 0.0003 on one row.
