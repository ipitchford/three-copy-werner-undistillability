# Historical construction account (not a portable replay claim)

This material was present in the 3 October 2026 manuscript. It describes the
producer's discovery and correction records. The final certificate verifier
does not reconstruct the baseline Grams, perturbations, or correction-system
minor discussed here. Those historical inputs are not shipped as a second
portable construction replay, and the supplied referee did not reproduce them.
The endpoint theorem uses only positivity of the final Grams and the complete
exact identity, which are separately replayed. The statements below are retained
as construction provenance, not as additional independently verified evidence.

For clarity, the final correction procedure is not an additional mathematical assumption. Starting with exact positive baseline Grams, it used 2,002 selected ordinary Gram-entry directions and 272 genuine equality rows, giving 2,274 selected columns. A square row minor had nonzero determinant \(664732\pmod{1000003}\). The exact solve was then checked on **all 2,374 rows**, with zero integer residual after denominators were cleared. The saved ordinary Gram-entry corrections were independently reapplied to the baseline blocks. All 64 changed Grams passed exact positive-definiteness checks. They are 64 changed blocks in total, not necessarily 64 coupled blocks. Untouched blocks retain their original exact positivity proofs.

There is also a uniform exact bound on the changed Grams. If \(Z_t^{(0)}\) is the positive definite baseline and \(\Delta_t\) its correction, rational comparisons verify

\[
\left\|(Z_t^{(0)})^{-1}\Delta_t\right\|_\infty<\frac18
\quad\text{for all 64 changed blocks}.
\tag{28a}
\]

The matrix in (28a) is similar to the symmetric relative perturbation \((Z_t^{(0)})^{-1/2}\Delta_t(Z_t^{(0)})^{-1/2}\). Its eigenvalues therefore lie strictly between \(-1/8\) and \(1/8\), and \(Z_t^{(0)}+\Delta_t\succ(7/8)Z_t^{(0)}\). The displayed bound is established by rational comparisons; the decimal maximum, approximately \(0.12407064145\), is only a readable summary.

An exact rank upper bound on the entire space of possible Gram corrections is unnecessary: the identity (28) is checked directly for the displayed certificate. A proof that a numerical candidate was close to feasible would not suffice; the exact zero residual and exact positivity checks are the load-bearing facts.

