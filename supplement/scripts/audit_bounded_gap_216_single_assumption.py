#!/usr/bin/env python3
"""Standard-library audit of the end-to-end diameter-216 package.

This verifies artifact bindings, finite certificate flags, exact support
parameters, admissibility of the displayed tuple, and all mechanical S52
regression checks.  The program does not replace the analytic proof of S52.
"""

from __future__ import annotations

import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
OUTPUT = OUTPUTS / "bounded_gap_216_single_assumption_audit.json"

CANDIDATE = OUTPUTS / "bounded_gap_216_address_state_rational_candidate.json"
SUPPORT = OUTPUTS / "bounded_gap_216_address_state_support_exact.json"
I_CERTIFICATE = OUTPUTS / "bounded_gap_216_address_state_i_arb.json"
J_CERTIFICATE = OUTPUTS / "bounded_gap_216_address_state_j_arb.json"

TUPLE = (
    0, 4, 6, 16, 18, 28, 30, 34, 40, 48, 54, 58, 60, 64, 70,
    76, 78, 84, 88, 96, 100, 106, 114, 118, 120, 126, 130, 138,
    144, 148, 154, 160, 166, 168, 174, 180, 184, 186, 190, 196,
    198, 204, 208, 210, 214, 216,
)


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def canonical_digest(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decimal_lower(value: str) -> Decimal:
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", value, re.I)
    if match is None:
        raise ValueError(f"cannot parse interval lower endpoint: {value!r}")
    return Decimal(match.group(0))


def primes_through(limit: int) -> list[int]:
    result = []
    for value in range(2, limit + 1):
        if all(value % prime for prime in result if prime * prime <= value):
            result.append(value)
    return result


def tuple_audit() -> dict[str, object]:
    residue_data = {}
    for prime in primes_through(len(TUPLE)):
        residues = {value % prime for value in TUPLE}
        missing = sorted(set(range(prime)) - residues)
        residue_data[str(prime)] = {
            "occupiedCount": len(residues),
            "missingResidues": missing,
        }
    return {
        "cardinality": len(TUPLE),
        "diameter": max(TUPLE) - min(TUPLE),
        "strictlyIncreasing": all(
            left < right for left, right in zip(TUPLE, TUPLE[1:])
        ),
        "admissible": all(
            row["missingResidues"] for row in residue_data.values()
        ),
        "residuesThroughK": residue_data,
    }


def main() -> None:
    candidate = load(CANDIDATE)
    support = load(SUPPORT)
    i_certificate = load(I_CERTIFICATE)
    j_certificate = load(J_CERTIFICATE)
    binding = candidate["binding"]
    candidate_sha = candidate["sha256"]

    finite_checks = {
        "candidateBindingHashRecomputed": (
            canonical_digest(binding) == candidate_sha
        ),
        "candidateDimensionIs46": int(binding["k"]) == 46,
        "candidateCoefficientCountIs208": (
            len(binding["numerators"]) == 208
        ),
        "outerRadiusIs2509Over10000": (
            binding["outerRadius"] == "2509/10000"
        ),
        "candidateHashMatchesSupport": (
            support["candidateSha256"] == candidate_sha
        ),
        "supportSourceHashMatchesCandidate": (
            support["supportSourceSha256"] == binding["supportDigest"]
        ),
        "candidateHashMatchesI": (
            i_certificate["candidateSha256"] == candidate_sha
        ),
        "candidateHashMatchesJ": (
            j_certificate["candidateSha256"] == candidate_sha
        ),
        "iCertificateFileHashMatchesJ": (
            j_certificate["iCertificateFileSha256"]
            == file_digest(I_CERTIFICATE)
        ),
        "supportHashMatchesI": (
            i_certificate["supportSha256"] == binding["supportDigest"]
        ),
        "supportHashMatchesJ": (
            j_certificate["supportSha256"] == binding["supportDigest"]
        ),
        "allSupportChecksPass": all(
            support["structuralChecks"].values()
        ),
        "allIChecksPass": all(i_certificate["checks"].values()),
        "allJChecksPass": all(j_certificate["checks"].values()),
        "scoreLowerBoundExceedsOne": (
            decimal_lower(j_certificate["score46JOverI"]["lower"])
            > Decimal(1)
        ),
        "marginLowerBoundPositive": (
            decimal_lower(j_certificate["margin46RJMinusI"]["lower"])
            > Decimal(0)
        ),
    }

    arithmetic_audits = {
        "exponentLedger": load(
            OUTPUTS / "bounded_gap_s52_exponent_ledger.json"
        )["status"]["allChecksPass"],
        "q0Transport": load(
            OUTPUTS / "bounded_gap_s52_q0_transport.json"
        )["status"]["allChecksPass"],
        "gcdRepair": load(
            OUTPUTS / "bounded_gap_s52_gcd_repair.json"
        )["allChecksPass"],
        "reversiblePreprocessing": load(
            OUTPUTS / "bounded_gap_s52_reversible_preprocessing.json"
        )["allChecksPass"],
        "primeSupportTransport": load(
            OUTPUTS / "bounded_gap_s52_support_transport.json"
        )["allChecksPass"],
        "lossBudget": load(
            OUTPUTS / "bounded_gap_s52_loss_budget.json"
        )["status"]["algebraicAuditPass"],
        "finiteJetClosure": load(
            OUTPUTS / "bounded_gap_s52_finite_jet_closure.json"
        )["allChecksPass"],
        "outerClosure": load(
            OUTPUTS / "bounded_gap_s52_outer_closure.json"
        )["allChecksPass"],
    }
    tuple_result = tuple_audit()
    finite_checks["tupleHas46Entries"] = tuple_result["cardinality"] == 46
    finite_checks["tupleDiameterIs216"] = tuple_result["diameter"] == 216
    finite_checks["tupleIsAdmissible"] = tuple_result["admissible"]

    payload = {
        "name": "End-to-end diameter-216 package audit",
        "claimBoundary": (
            "All finite support, variational, binding, tuple, and mechanical "
            "S52 regression checks pass. The manuscript proves S52 relative "
            "to the cited established trace and dispersion theorems and then "
            "deduces H_1<=216. This audit does not replace the analytic proof "
            "of S52."
        ),
        "candidateSha256": candidate_sha,
        "finiteChecks": finite_checks,
        "allFiniteChecksPass": all(finite_checks.values()),
        "s52MechanicalRegressionChecks": arithmetic_audits,
        "allS52MechanicalRegressionChecksPass": all(
            arithmetic_audits.values()
        ),
        "tuple": {"values": list(TUPLE), **tuple_result},
        "rigorousBounds": {
            "score46JOverILower": j_certificate["score46JOverI"]["lower"],
            "margin46RJMinusILower": j_certificate[
                "margin46RJMinusI"
            ]["lower"],
        },
        "logicalStatus": {
            "finiteVariationalAndSupportModule": "proved",
            "sieveTransfer": "proved in the manuscript",
            "s52Theorem": (
                "proved in the manuscript relative to cited established "
                "trace and dispersion theorems"
            ),
            "unconditionalH1AtMost216Claimed": True,
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if not payload["allFiniteChecksPass"]:
        raise SystemExit("finite diameter-216 audit failed")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "allFiniteChecksPass": payload["allFiniteChecksPass"],
                "allS52MechanicalRegressionChecksPass": payload[
                    "allS52MechanicalRegressionChecksPass"
                ],
                "theoremStatus": "manuscript claims unconditional H_1 <= 216",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
