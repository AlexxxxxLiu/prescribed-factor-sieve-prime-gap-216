# Ancillary verification package

The `outputs` directory contains the frozen rational candidate, the exact
support certificate, the 1024-bit Arb certificates, and the finite regression
outputs.  The `scripts` directory contains every local Python module required
by the reproduction commands in Appendix B of the paper.

The package certifies a frozen rational candidate.  It does not include the
floating-point optimizer output that originally discovered that candidate, so
``reproduction'' means replay of the exact support and interval certificates,
not regeneration of the heuristic search history.  The absolute source path
in the candidate's provenance field is nonportable metadata and is excluded
from the canonical mathematical digest.

Run all commands from this directory so that the scripts resolve `outputs/`
relative to the package root.  Set `PYTHONPATH=scripts` for scripts that import
the shared exact-arithmetic modules.

The Z3 support checker records the complete rational input and partition
cover, but the present artifact does not include a separately replayable Z3
proof object.  The Arb computations use directed outward rounding.

The manuscript proves the analytic prescribed-factor Type-I theorem
separately.  These programs check finite algebra, support routing, exponent
bookkeeping, and numerical certificates; they are not substitutes for the
analytic proof.  The two components are combined in the paper to obtain its
unconditional diameter-216 conclusion.
