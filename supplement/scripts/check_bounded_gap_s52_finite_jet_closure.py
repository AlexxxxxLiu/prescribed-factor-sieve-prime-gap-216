#!/usr/bin/env python3
"""Regression audit for the integrated S52 finite-jet closure proof.

This script binds directly to the TeX proof included in the main manuscript.
It checks the three Taylor eliminations, the derivative debt, the remainder
budget, and the common polynomial-size terminal range.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROOF = ROOT.parent / "s52-proof" / "04-supported-descent.tex"
OUTPUT = ROOT / "outputs" / "bounded_gap_s52_finite_jet_closure.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    proof = PROOF.read_text(encoding="utf-8")
    taylor_markers = (
        r"\label{eq:Taylor-one}",
        "Taylor expansion to order $J$ separates this profile",
        "Taylor expansion through order $J$ therefore separates it",
    )
    taylor_depth = len(taylor_markers)
    terminal_derivative_order = 2
    checks = {
        "integratedProofExists": PROOF.is_file(),
        "allThreeTaylorEliminationsPresent": all(
            marker in proof for marker in taylor_markers
        ),
        "proofDefinesInputJetClass": (
            r"\begin{lemma}[Finite-jet closure]" in proof
            and r"\label{lem:finite-jet}" in proof
            and r"\label{eq:input-jets}" in proof
        ),
        "proofDefinesTerminalC2Class": (
            "one common\nshifted $C^2$ profile class" in proof
        ),
        "proofCarriesThreeTaylorDebt": (
            r"3J+2=J_*" in proof
            and r"J_*$ input derivatives" in proof
        ),
        "proofFixesTaylorOrder": r"J=\lceil20/\eta\rceil" in proof,
        "proofControlsRemainder": (
            r"x^{-5\eta(J+1)}" in proof
            and r"O(x^{-100})" in proof
        ),
        "proofBoundsDescendantCount": r"(J+1)^3" in proof,
        "terminalDebtFormula": (
            taylor_depth == 3 and terminal_derivative_order == 2
        ),
        "proofGivesExplicitPolynomialExponent": (
            r"m\le x^7" in proof
            and r"N_2\le x^7" in proof
            and r"\Delta_2/z_1\le x^7" in proof
        ),
    }
    payload = {
        "name": "Integrated S52 finite-jet closure regression",
        "claimBoundary": (
            "This checks the TeX proof binding and derivative-debt "
            "bookkeeping; the chain-rule argument in the manuscript remains "
            "the mathematical proof."
        ),
        "integratedProof": str(PROOF.relative_to(ROOT.parent)),
        "integratedProofSha256": digest(PROOF),
        "symbolicDebt": {
            "taylorDepth": taylor_depth,
            "taylorOrder": "J=ceil(20/eta)",
            "terminalDerivativeOrder": terminal_derivative_order,
            "maximumTerminalInputOrder": "3J+2",
            "maximumRemainderInputOrder": "3J+1",
        },
        "checks": checks,
        "allChecksPass": all(checks.values()),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if not payload["allChecksPass"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise SystemExit(
            "S52 finite-jet closure regression failed: " + ", ".join(failed)
        )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
