# Agent-readable research index

## Read first

This is an unrefereed candidate **computer-assisted proof**, not an unresolved
numerical search. Read [STATUS.md](STATUS.md), [the complete manuscript](paper.md),
and [the machine claim scope](CLAIMS.json) before reusing the conclusion.
The theorem covers all complex rank-at-most-two 27×27 matrices, and only the
specified three-copy qutrit endpoint. Corollary 7 gives the physical alpha threshold.

## Proof and exact object

The written semantic bridge is in [paper.pdf](paper.pdf) and [paper.tex](paper.tex).
The finite object is [certificate.json](certificate/certificate.json), with
[rational data](certificate/data/), [final multipliers](certificate/multipliers/),
and a [manifest](certificate/manifest.json). No numerical solver or historical
construction baseline is needed to check the theorem. See
[mathematical data identity](MATHEMATICAL_DATA_IDENTITY.json).

## Replay and ways to falsify acceptance

Follow [README.md](README.md), [the exact verifier instructions](certificate/README.md)
and [the acceptance record description](REPLAY_RECEIPT.md). Use new empty output
directories for normal and optimized runs. The [semantic controls](check_controls.py)
must reject negative Grams, an altered objective and a missing coupled cross factor.
The [review witness checker](check_review_congruences.py) uses exact congruences.
The [package checker](verify_package.py) checks all manifest entries.

## Reviews, sources and trust boundaries

Read [ASSURANCE.md](ASSURANCE.md), [the supplied referee report](reviews/referee-report-2026-10-04.md),
[the supplied audit](reviews/referee-audit-2026-10-04/README.md), and the
[point-by-point response](RESPONSE_TO_REFEREE.md). The alternative arithmetic
retains author foundation/coordinate definitions; human independence is not authenticated.
[Historical construction](HISTORICAL_CONSTRUCTION.md) is outside final replay coverage.
[CITATION_AUDIT.md](CITATION_AUDIT.md) and [SOURCES.md](SOURCES.md) distinguish
inherited methods and the scoped result. No exhaustive priority claim is made.

## Identity and reuse

[PROVENANCE.md](PROVENANCE.md), [CITATION.cff](CITATION.cff),
[licence map](LICENSES.md), and [environment](ENVIRONMENT.txt) specify identity,
source hashes, original-versus-supplied components and runtime assumptions.
The public receipt is a separate release asset so it can bind the final ZIP
without self-reference. Public availability never upgrades scientific assurance.

- [Contribution and limitation record](RESEARCH_GATES.json): retrospective intake search, exact source locators, closed extensions and remaining scope.
