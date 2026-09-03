# A Prescribed-Factor Sieve and an Unconditional Prime-Gap Bound of 216

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22261945.svg)](https://doi.org/10.5281/zenodo.22261945)

**Author:** Jingwei Liu  
**Archived release:** [Zenodo 10.5281/zenodo.22261945](https://doi.org/10.5281/zenodo.22261945)  
**License:** [CC BY 4.0](LICENSE)

This repository contains the LaTeX source, compiled paper, and reproducibility
materials for Jingwei Liu's prime-gap manuscript.

The archived Zenodo release is the fixed scholarly record. This GitHub
repository provides a browsable source tree and a convenient route for future
versioned corrections or extensions.

## Build the paper

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The resulting file is `main.pdf`.

## Reproduce the finite certificate

The `supplement` directory is self-contained for replaying the frozen finite
candidate apart from its Python dependencies.  From that directory:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
PYTHONPATH=scripts .venv/bin/python scripts/check_bounded_gap_226_layered_support.py \
  --k 46 \
  --candidate-json outputs/bounded_gap_216_address_state_rational_candidate.json \
  --output outputs/bounded_gap_216_address_state_support_exact.json
PYTHONPATH=scripts .venv/bin/python scripts/certify_bounded_gap_216_address_state_i_arb.py \
  --precision 1024
PYTHONPATH=scripts .venv/bin/python scripts/certify_bounded_gap_216_address_state_j_arb.py \
  --precision 1024
PYTHONPATH=scripts .venv/bin/python scripts/audit_bounded_gap_216_single_assumption.py
```

The generic support-checker filename retains an earlier `226` development
label.  The explicit `--k 46` and frozen candidate bind the calculation to
the diameter-216 certificate.
The floating-point search output that discovered the candidate is not included;
the package replays and certifies the frozen rational candidate rather than
reproducing its heuristic discovery process.  The provenance path stored in
the candidate JSON is nonportable metadata and is not part of its canonical
mathematical digest.

## Verification boundary

The exact-support and Arb certificates verify the finite `k=46` variational
module and the admissible diameter-216 tuple.  The regression scripts verify
finite algebra and artifact bindings in the prescribed-factor inheritance
chain; they do not replace its analytic proof.  The local descent, terminal
closure, outer dispersion closure, and the resulting prescribed-factor
Type-I theorem are proved in the body of the manuscript.  Combining that
theorem with the finite certificate and the sieve-transfer proposition gives
the manuscript's unconditional conclusion `H_1 <= 216`.

## Citation

Please cite the archived release:

```bibtex
@misc{liu2026prescribed,
  author    = {Jingwei Liu},
  title     = {A Prescribed-Factor Sieve and an Unconditional Prime-Gap Bound of 216},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22261945},
  url       = {https://doi.org/10.5281/zenodo.22261945}
}
```
