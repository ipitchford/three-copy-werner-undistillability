# Three-copy qutrit Werner undistillability

Anonymous · Version 1.1.0-candidate · 4 October 2026 · DOI [10.5281/zenodo.23134454](https://doi.org/10.5281/zenodo.23134454)

This candidate contains a computer-assisted proof that every complex 27×27
matrix of rank at most two satisfies the three-copy Werner endpoint inequality.
There is no normality, Hermiticity or local-support restriction. A written
corollary establishes the three-copy threshold for the physical qutrit Werner
family: distillable exactly for −1 ≤ α < −1/2. No four-copy, all-copy or
unrestricted higher-dimensional theorem is claimed.

Start with the [paper PDF](paper.pdf), [accessible manuscript](paper.md), and
[agent index](AI_INDEX.md). The mathematical result consists of the written
reduction and its exact rational certificate, not the release badge or DOI.
Publication status remains **unrefereed candidate**; the supplied referee audit
has an explicitly bounded implementation-diversity scope and does not establish
authenticated human specialist or journal peer review.

## Reproduce

Use Python 3.11+ on Linux or macOS and install `certificate/requirements.txt`.
From outside the package, choose a new empty output directory:

```sh
python -B -s /absolute/package/certificate/verify.py --output /absolute/new-normal --phase all --max-seconds 1200
python -O -B -s /absolute/package/certificate/verify.py --output /absolute/new-optimized --phase all --max-seconds 1200
python -B -s /absolute/package/check_controls.py --replay /absolute/new-normal
python -O -B -s /absolute/package/check_controls.py --replay /absolute/new-optimized
python -B -s /absolute/package/check_review_congruences.py
python -O -B -s /absolute/package/check_review_congruences.py
```

Reserve one CPU thread and about 2.5 GB RAM. In the original shared environment,
admit the complete command tree using `research-run run --cpu 1 --memory-gb 2.5`.
Successful completion requires `COMPLETE_PORTABLE_EXACT_CERTIFICATE_PASS`, all
716 PSD checks and zero exact residual in all 2,374 coordinates. A small residual,
partial receipts or a resumed run is not a fresh independent replay. Runtime
depends on the host; the original endpoint replay took about 139 seconds.

See [certificate instructions](certificate/README.md) for bounded/resumed runs
and platform assumptions. [Release verification record](REPLAY_RECEIPT.md)
distinguishes historical evidence from the new frozen-package replay.

## Package map

- `paper.md`, `paper.tex`, `paper.pdf`: full scientific manuscript.
- `certificate/`: final rational data and exact verifier.
- `check_controls.py`: deliberate mathematical errors and resource-unit checks.
- `check_review_congruences.py`: exact recheck of supplied independent-arithmetic witnesses.
- `reviews/`: supplied review and audit, preserved without relicensing.
- `RESPONSE_TO_REFEREE.md`: disposition of every requested revision.
- `HISTORICAL_CONSTRUCTION.md`: optional construction history, outside final replay coverage.
- `CLAIMS.json`, `ASSURANCE.md`, `PROVENANCE.md`, `CITATION_AUDIT.md`: scope and evidence.
- `LICENSES.md`: CC0 original prose/data, MIT original code, retained third-party terms.

The public release assets include `verification-report.json` and `SHA256SUMS`;
these are produced after replay of the frozen archive, avoiding a circular
self-hash. Hosted CI checks the immutable payload commit. Provider availability
and website deployment are separate from scientific verification.
