# Referee audit: three-copy qutrit Werner endpoint

Review date: 4 October 2026.
Input: werner-problem6-proof.zip.

## Outcome

The alternative exact arithmetic replay passed. All 716 small multiplier
Grams were shown positive definite by exact integer congruence and diagonal
dominance. All 844 family functional vectors were recomputed from the final
N,Z data, local modules and source generators. All 2,374 coordinates of the
final rational identity had zero integer residual. The common denominator had
3,873 bits. The fresh replay took approximately 126.7 seconds on this host.
This is a run-specific measurement, not a general performance guarantee.

The native author verifier was NOT run: python-flint was unavailable and
installation/download attempts failed. The alternative replay retains the
submitted foundation and coordinate definitions; it is NOT a wholly separate
second implementation of the mathematical reduction. It replaces the modular
determinant, rational tensor arithmetic, PSD test and final summation backends.
Saved author PASS receipts are not trusted to establish any of these checks.
Author functional arrays are compared only after fresh recomputation.

An additional independently written test expands all 2,374 selected source
generators directly from the literal group-average definition, confirming
bosonic sandwiches, transpose and physical-site orbit weights. Small exact
Gaussian-integer tests check 24 complex nonnormal rank-at-most-two examples,
192 swap/partial-trace identities and 240 literal coupled cross-map entries.
Those finite tests are convention checks, not substitutes for the universal
certificate.

## Files

- summary.json: final alternative-replay result and explicit scope limits.
- hash_check.json: input bindings and complete manifest checks.
- foundation.json: recomputed physical embeddings, local metrics, source
  dimensions, modular basis minors, objective and selected equality rows.
- reviewer_checks.py: alternative exact arithmetic replay.
- replay.log: complete progress/output log from that replay.
- psd_congruence_certificates.json: all 716 explicit integer congruences R.
- verify_saved_congruences.py: fresh verification of those congruences using
  exact arithmetic only, with no author Python imports and no floating point.
- physical_checks.py / physical_checks.json: independent exact convention tests.
- projector_checks.py / projector_checks.json: independent literal expansion
  of every selected source generator.
- parameter_corollary.md: reviewer-derived full three-copy qutrit parameter
  threshold, conditional only on the submitted endpoint theorem.
- findings.md: editorial and claim-scope findings.
- literature_checked.json: bounded primary-source search record.
- environment.json: reviewer execution environment.

## Reproduction

Use Python 3.11 or later, NumPy, SymPy and SciPy. The successful review used
Python 3.13.5, NumPy 2.3.5, SymPy 1.14.0 and SciPy 1.17.0, without FLINT.
Do not use Python's -O option (some auxiliary checks use assertions).

The full alternative script expects the original ZIP two directories above
the extracted package root, in order to record its hash. A suitable layout is:

    work/werner-problem6-proof.zip
    work/extracted/werner-problem6-proof/
    work/referee-audit/  (this directory)

Run with a new output directory:

    python -B reviewer_checks.py /absolute/path/work/extracted/werner-problem6-proof /absolute/path/fresh-reviewer-output

The standalone PSD and projector checks accept the extracted package root:

    python -B verify_saved_congruences.py /absolute/path/work/extracted/werner-problem6-proof
    python -B projector_checks.py /absolute/path/work/extracted/werner-problem6-proof
    python -B physical_checks.py

The last two auxiliary scripts write their small JSON outputs next to the
scripts. The submitted package itself is never edited.

## What was not established

The discovery search was not rerun. Historical baseline Grams and correction
operations were not reconstructed. In particular, the optional relative
perturbation bound <1/8 and the correction-system minor in the submitted
construction narrative were not freshly reproduced: final corrected Grams
were verified directly instead. Literature priority was not established
exhaustively. There is no proof-assistant formalisation, independent external
human peer-review claim, or all-copy bound-entanglement conclusion here.
