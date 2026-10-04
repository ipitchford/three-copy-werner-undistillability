# Portable rational certificate for the three-copy Werner inequality

This directory contains the final rational certificate data and an exact-arithmetic
verifier. It is self-contained: historical solver files, floating-point models,
numerical changes of basis, and the zero-family calculations used to discover the
certificate are not inputs to verification. Historical paths in `certificate.json`
are provenance only.

**Acceptance rule:** a fresh replay succeeds only when it writes
`COMPLETE.json` with status `COMPLETE_PORTABLE_EXACT_CERTIFICATE_PASS` in a new
receipt directory outside this input directory. The original construction export's pending status
is a historical marker, not a substitute for this replay. A partial run, source
review, or numerical residual is not a complete verification.

## Run

Use Python 3.11 or later with the dependencies in `requirements.txt`. The listed
versions are the environment in which the verifier was developed; the exact
integer and rational operations do not use floating-point tolerances.

From a directory other than this input directory, with a fresh output path:

```sh
env -u PYTHONPATH python -B -s /absolute/path/certificate/verify.py \
  --output /absolute/path/new-replay --phase all --max-seconds 1200
```

The input directory is never modified. `-B` prevents Python bytecode files there;
`-s` excludes the user site-package directory. Reserve one CPU thread and about
2.5 GB RAM for a full replay. In the original shared environment, run this command
through `research-run` using the corresponding resource reservation.

For bounded runs, invoke `--phase foundation`, then repeat `--phase blocks`
with `--max-seconds 150`, and finish with `--phase finalize`. Successful block
receipts are atomic checkpoints bound to the input manifest. Disjoint
`--start N --stop M` block ranges can be processed separately in the same receipt
directory after the foundation phase. A timeout preserves committed blocks.
For an independent replay, always start from an empty receipt directory: resuming
intentionally trusts earlier receipts from that same replay.

## What is verified

1. Every input file is bound by `manifest.json` and its SHA-256 digest.
2. Every local rational module is checked against literal tensor-permutation
   matrices and their specified partial transpose, including the exact physical
   embedding and positive Gram metric. The joint off-diagonal permutation is
   also checked on literal tensor indices.
3. The rational source generators span spaces of dimensions 577, 1220, and 577.
   Exact averaging-projector traces give upper bounds, while actual projected
   generator minors over a finite field give matching lower bounds.
4. The objective and all 272 selected homogeneous equality rows are recomputed
   from the tensor diagrams. The saved equality matrix `B` is an integer matrix;
   the multiplier `eta` multiplies `B`, not `B / 1152`.
5. All 716 nonzero final rational Gram matrices `Z` are checked positive
   semidefinite. Symmetry and the alternating signs of the integer characteristic
   polynomial prove this exactly. The physical multiplier is `D = N Z N^T`.
   The remaining 741 multipliers are explicitly omitted as zero.
6. All multiplier functionals are recomputed on all 2,374 rational basis
   generators. The joint cross term is kept nonsymmetric and has its required
   factor of two. Exact integer arithmetic then verifies

   ```text
   c = sum(multiplier functionals) + B^T eta.
   ```

The source families are ordered `(3,1), (2,2), (1,3)`, with replica layouts
`A0 A1 A2 B3`, `A0 A1 B2 B3`, and `B0 B1 B2 A3`. The `A` factor is the normalized
isometry. Source coordinates represent `54^4` times the normalized moments, and
the objective includes that normalization. No additional scale factor is to be
applied to the exported multipliers.

## Contents

- `certificate.json`: coverage, source conventions, final multiplier references,
  explicit zero omissions, and provenance.
- `multipliers/`: 652 source and 64 joint nonzero final multipliers, including all
  exact affine corrections.
- `data/local_modules.json`: literal rational module embeddings and metrics.
- `data/basis_data.npz`: selected rational diagram generators, finite-field minor
  data, exact objective, and homogeneous equality rows.
- `data/equality_multiplier.json`: exact rational `eta`.
- `foundation.py`, `exact_coordinates.py`: physical and coordinate foundations.
- `exact_tables.py`, `basis_values.py`: exact tensor contractions and basis values.
- `verify.py`: checkpointed replay and final exact identity verification.

The external mathematical proof explains why these positive moment functionals
and valid homogeneous equalities imply the unrestricted rank-two inequality.
This bundle is the inspectable finite arithmetic part of that proof, rather than
a proof of any unrelated higher-copy or stronger inequality.

## Platform and revision notes

Version 1.1.0-candidate preserves every rational mathematical input from the
3 October certificate. Historical-status naming, resource reporting and explicit consistency checks
changed. The original archive and manifest are retained by hash in PROVENANCE.md.

The verifier targets Unix-like systems (tested on macOS; Linux is exercised in
CI). It uses `resource`, `SIGALRM`, and same-directory hard links for atomic
receipt commits. Native Windows is not supported; use a suitable Unix environment.
macOS reports `ru_maxrss` in bytes and Linux in KiB; the revised reporting helper
converts those platforms to bytes and retains raw units. Unknown platforms report
`peak_rss_bytes: null` instead of guessing. This measurement is process peak RSS,
not aggregate tree memory; the separate supervisor receipt reports sampled tree RSS.

Always start a separate replay with a new empty receipt directory. Normal and
optimized Python modes have the same explicit `require` acceptance predicates.
