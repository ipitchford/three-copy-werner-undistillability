# Referee report

**Submission:** *Three-copy undistillability of the boundary qutrit Werner state*  
**Review date:** 4 October 2026  
**Editorial recommendation:** **Major revisions—principally scholarly positioning and ancillary reproducibility claims, not replacement of the central proof.**  
**Prospective REF assessment:** **4\* candidate, with moderate confidence.**

## 1. Executive assessment

**This is a substantive computer-assisted mathematical result, not merely a numerical candidate.** I found no mathematical defect in the central reduction, and reproduced the decisive certificate using alternative exact arithmetic.

The fresh calculation verified all **716 multiplier Grams**, recomputed their functionals, and obtained **exactly zero residual in all 2,374 coordinates** of the final identity. Separate checks confirmed the physical permutation conventions and independently expanded every selected source generator. The unrestricted quantifier appears justified: the proof covers arbitrary complex, nonnormal, non-Hermitian rank-at-most-two matrices, rather than only special support configurations.

The distinction between verification and discovery is handled particularly well. Neither favourable optimiser output nor completeness of the zero families used during discovery is required for the final argument. The resulting theorem remains a finite-copy result: it does not establish undistillability for arbitrarily many copies or settle general negative-partial-transpose bound entanglement.

The principal shortcomings concern the article surrounding the proof. Its main text needs a proper account of the existing literature, and some optional claims about the historical correction procedure are not reproduced by the portable verifier. These are repairable without obtaining another major theorem.

**My assessment is that the supplied mathematics already supports a credible world-leading research candidate. Publication should follow revision and specialist external scrutiny, not another campaign of numerical searches.**

werner_referee_audit_2026-10-04.zip[Download the reviewer audit bundle](sandbox:/mnt/data/werner_referee_audit_2026-10-04.zip). It contains the checking scripts, exact positivity certificates, results, limitations and a useful parameter-threshold corollary.

## 2. Scope and limits of this review

I examined the complete main proof, the package and certificate documentation, all five verification modules, the certificate structure, relevant supporting audits, and the original problem statement and earlier partial-results note retained under `inputs/`. Manifest checks covered **2,181 packaged files**, including all **727 certificate inputs**.

Manuscript references below are to `proof.md`, with section and equation numbers.

There is an important computational qualification. **I did not run the unmodified FLINT-based verifier.** This environment lacked `python-flint`, and attempts to install or download it failed. Instead, I implemented alternative exact arithmetic for the modular determinants, rational tensor contractions, multiplier positivity and final functional summation.

That replay retained the submission’s foundation and coordinate definitions. It is therefore **not a wholly independent second implementation of the entire mathematical proof**. However, it did not trust the supplied pass receipts. The archived functional arrays were compared only after fresh recomputation. A further independently written check reconstructed all selected source generators directly from the literal group-average definition.

I did not rerun the discovery optimisation or reconstruct historical baseline Grams and correction steps. Those are unnecessary for verifying the final theorem, but this limits what I can endorse about the accompanying construction narrative.

The literature review was a bounded primary-source check. I inspected relevant theorem statements and scope restrictions, not independently reproved every external paper.

## 3. Contribution and positioning

The central assertion is

\[
q_3(C)=\|C\|_F^2
-\frac12\sum_i\|\operatorname{Tr}_i C\|_F^2
+\frac14\sum_{i<j}\|\operatorname{Tr}_{ij}C\|_F^2
-\frac18|\operatorname{Tr}C|^2\geq0
\]

for every complex \(27\times27\) matrix \(C\) of rank at most two. Through the stated vectorisation, this establishes three-copy undistillability of

\[
\rho=\frac{2I-F}{15}.
\]

The scientific contribution should be framed as **an explicit exact certificate for the unrestricted three-copy qutrit endpoint**.

It should not be framed as the invention of the partial-trace reformulation. Costa Rico’s 2025 paper already develops that connection between Werner-state distillability and rank-constrained partial-trace inequalities. Likewise, the maximally mixed auxiliary-qubit formulation has a close antecedent in Bharti, Gajjala and Haug’s 2026 treatment. [Springer](https://link.springer.com/article/10.1007/s11005-025-01935-y)

The nearest three-copy comparison is Wu and Zou. Their preprint proves the normal rank-two case in arbitrary local dimensions, together with specified nonnormal support sectors. It does not establish the unrestricted qutrit statement proved here. Removing those restrictions is a substantive advance, not just a larger numerical sample. [arXiv](https://arxiv.org/html/2608.02647v1)

The methodological contribution is more specific than “use semidefinite programming”. Entanglement-witness and robust-semidefinite approaches to Werner distillability substantially predate this submission. What matters here is the particular matched-moment construction, the coupled Gram constraint, and the successful conversion into a finite rational certificate. [arXiv](https://arxiv.org/abs/quant-ph/0608095)

An especially valuable feature is that the certificate is attached to a mathematical sufficiency argument. It is not merely an internally consistent collection of positive matrices: the paper explains why arbitrary admissible physical inputs satisfy the inequalities being combined.

## 4. Major comments

### 1. Bring the literature and novelty argument into the main article

**Issue and evidence.** The main proof has no conventional related-work section or bibliography. Relevant references occur in the historical partial-results note inside the nested input archive, but a reader should not have to unpack that archive to understand the contribution.

**Why this matters.** The distinction between inherited reductions and the genuinely new result is central to assessing originality. The recent two-copy papers, auxiliary-qubit formulations and restricted three-copy results are directly relevant. In particular, the July 2026 literature already reports a complete two-copy threshold; presenting this work against an older two-copy frontier would materially misstate its context. [arXiv](https://arxiv.org/abs/2607.21367)

**Required revision.** Add an introduction that identifies the exact remaining quantifier closed by this certificate. It should distinguish the general partial-trace reduction, the auxiliary-qubit lifting, previous two-copy results, previous three-copy support restrictions, and the explicit certificate supplied here.

Also separate **novelty of the theorem** from **novelty of the general techniques**. The theorem can be important even where several ingredients are established methods.

The manuscript already disclaims a priority claim. I found no evidence of deliberate misrepresentation; this is a substantial positioning omission.

### 2. Separate final-certificate verification from historical construction claims

**Issue and evidence.** Section 5 includes the relative correction bound

\[
\left\|(Z_t^{(0)})^{-1}\Delta_t\right\|_\infty<\frac18
\tag{28a}
\]

and describes a particular correction-system minor and reconstruction procedure. The supporting uniform-bound receipt refers to historical `work/...` inputs. The portable package supplies the final corrected Grams, but the portable verifier does not reconstruct the historical baselines and perturbations needed to recheck this relative bound.

**Why this matters.** Verifying the final corrected \(Z_t\) proves their positivity. It does **not**, by itself, verify a comparison between those matrices and absent historical baselines. These are different assertions.

The mathematical implication from the displayed relative bound to

\[
Z_t^{(0)}+\Delta_t\succ \frac78Z_t^{(0)}
\]

is sound. The issue is reproducibility of the numerical premise in exact arithmetic, not the implication.

**Required revision.** Either provide the minimal baseline data and an exact checker for these optional construction claims, or move them into a clearly labelled historical-construction account and narrow the verification language accordingly.

**This does not undermine the endpoint theorem.** My replay directly established positivity of every final Gram and the complete functional identity, without using the historical perturbation bound or correction-system minor.

## 5. Rigour, computations and inference

### 5.1 The physical reduction preserves the full problem

The normalisation and factorisation in §1 are correct. After setting \(\|C\|_F^2=1/2\), choosing an isometry \(A\) containing the column space of \(C\), and defining \(B=2C^*A\), one obtains

\[
A^*A=I_2,\qquad C=\frac12AB^*,\qquad \|B\|_F^2=2.
\]

Crucially, **\(B\) is not required to be an isometry**. Imposing that additional condition would exclude general singular-value configurations. The manuscript avoids that mistake. Rank-one inputs are also covered by extending their column space to a two-plane.

The conversion factor

\[
\langle\psi_C|(\rho^{T_B})^{\otimes3}|\psi_C\rangle
=\frac8{3375}q_3(C)
\]

is consistent with the unnormalised maximally entangled vector convention.

### 5.2 The real symmetry reduction does not restrict \(C\) to real matrices

The averaging in §4 acts on the matched moment tuple. It does not assert that a general complex \(C\) can be replaced by a real or normal matrix.

Each transformed physical pair remains admissible, the objective is preserved, and averaging with the complex-conjugate tuple preserves positivity and the homogeneous equalities. Therefore a hypothetical negative complex input would still produce a negative feasible averaged tuple.

This is the necessary direction of argument, and it is correctly stated.

### 5.3 The coupled Gram map is a genuine physical positivity constraint

Equations (17)–(20) are important. With

\[
u=a\otimes a\otimes b\otimes\overline a,\qquad
v=b\otimes b\otimes a\otimes\overline b,
\]

the block matrix is exactly

\[
|u\oplus v\rangle\langle u\oplus v|.
\]

The literal cross-index rule produces \(u_i\overline{v_j}\). It is not an assumed positivity property of an arbitrary realignment.

The implementation also retains the nonsymmetric cross block. The factor of two in the cross contribution comes from the two adjoint off-diagonal terms; it is not an unexplained normalisation adjustment.

I checked this algebraically and through separate exact Gaussian-integer index tests.

### 5.4 Basis completeness is established, not inferred from numerical rank

The source-space dimensions were reproduced as

\[
\frac{33\,264+8\,280}{72}=577,\qquad
\frac{33\,504+5\,536}{32}=1220.
\]

The corresponding modular determinants were nonzero at the stated prime \(1\,000\,003\).

The logic is sound: independence after a fixed linear compression implies independence before compression; a nonzero determinant modulo a prime establishes rational independence. Combined with the exact dimension upper bounds, this yields complete bases.

The compressed row selection does not introduce a probabilistic assumption into the final proof.

I additionally reconstructed **all 2,374 selected generators** using an independently written expansion of the bosonic sandwich, transpose and physical-site average. Their coefficients agreed exactly with the compiled evaluation rules.

### 5.5 Positivity is connected to the physical operators

The identities

\[
P^{\Gamma_T}J_t=J_t\pi_t(P^{\Gamma_T}),
\qquad G_t=J_t^TJ_t\succ0
\]

are a major strength. They identify the represented operators as literal physical compressions, rather than merely supplying abstract matrices with plausible multiplication rules.

Consequently, the use of \(G_t\pi_t(X)\), rather than treating \(\pi_t(X)\) as symmetric in the ordinary Euclidean metric, is justified.

Two possible objections are correctly avoided. The selected 272 equality rows need only be genuine identities on physical tuples; they need not span every raw equality row. Likewise, the zero families used to discover multiplier ranges need not be complete once the final positive multipliers satisfy the global identity.

### 5.6 Exact replay results

| Check | Result of this review |
|---|---|
| Package and certificate hashes | All checked hashes matched |
| Literal local embeddings and positive metrics | Recomputed; 44 local modules passed |
| Source dimensions and basis minors | Recomputed; \(577,1220,577\) confirmed |
| Literal projected-generator expansions | All 2,374 checked independently |
| Final multiplier Grams | All 716 proved positive definite |
| Multiplier family functionals | All 844 recomputed and matched |
| Global certificate identity | All 2,374 integer residuals exactly zero |
| Final common denominator | 3,873 bits |

The functional count is consistent with the block structure:

\[
652\times1+64\times3=844.
\]

For positivity, I used a different exact test from the submission’s characteristic-polynomial method. For each rational Gram \(Z\), the checker constructed an invertible integer lower-triangular matrix \(R\) such that, after clearing a positive denominator,

\[
H=RZ_{\rm int}R^T
\]

was symmetric with positive, strictly diagonally dominant rows. This proves positive definiteness by congruence.

Floating-point calculations were used only to **propose** \(R\). Acceptance depended entirely on exact integer inequalities. A separate checker then reverified every saved congruence with no floating-point arithmetic and no imports from the author’s Python modules.

The submission’s own characteristic-polynomial criterion is also mathematically valid for symmetric matrices. There is no reason to reject it simply because it is less familiar than an \(LDL^T\) test.

Positive definiteness of the small \(Z_t\) does not imply that each full \(D_t=N_tZ_tN_t^T\) is positive definite: the rectangular factors can leave kernels. The proof needs only positive semidefiniteness of the full multipliers.

### 5.7 Independent convention checks and boundaries

The additional exact tests covered 24 complex nonnormal rank-at-most-two examples, 192 swap/partial-trace identities, and 240 indexed coupled cross-map identities.

These are calibration tests, not the universal proof.

I also checked simple exact boundary examples. Writing \(E_{ij}=|i\rangle\langle j|\) and \(D=\operatorname{diag}(1,1,0)\),

\[
q_3(D\otimes E_{00}\otimes E_{00})=0,
\]

and the nonnormal example \(D\otimes E_{01}\otimes E_{00}\) also has value zero. In contrast,

\[
q_3(I_3\otimes E_{00}\otimes E_{00})=-\frac38.
\]

The latter has rank three. It provides a useful check that the rank restriction is essential and that the tested quadratic form has not inadvertently been replaced by an everywhere-positive expression.

## 6. External literature check

The search was conducted on **4 October 2026**, using live web discovery, arXiv records and full-text HTML, Springer’s article text, and official REF guidance. Queries included the exact manuscript title, “Werner three-copy 2026 undistillability”, the title of Costa Rico’s paper, and current REF panel-criteria terms. Citation-following was used because some broad searches returned poor matches.

**Costa Rico (2025):** establishes the relevant partial-trace approach and earlier restricted results. This is a direct antecedent of the mathematical formulation, not an earlier proof of the unrestricted result supplied here. [Springer](https://link.springer.com/article/10.1007/s11005-025-01935-y)

**Fu, Gao and Park (2026):** report the complete two-copy result in arbitrary dimension. I checked the current arXiv version record and main result. Their result is relevant context, but does not establish this three-copy theorem. [arXiv](https://arxiv.org/abs/2607.21367)

**Bharti, Gajjala and Haug (2026):** provide two-copy inequalities and finite-copy reformulations, including the auxiliary-qubit viewpoint. They explicitly explain why the two-copy argument does not automatically induct. [arXiv](https://arxiv.org/html/2607.24479v1)

**Wu and Zou (2026):** are the closest three-copy comparator. Their normal and specified nonnormal sectors leave a residual unrestricted problem; this submission addresses that remaining quantifier in the qutrit case. [arXiv](https://arxiv.org/abs/2608.02647)

**Vianna and Doherty (2006):** establish earlier use of entanglement witnesses and robust semidefinite programmes for Werner distillability. I inspected the abstract and bibliographic record, not the full technical argument, and make no assessment here of its numerical-certification claims. [arXiv](https://arxiv.org/abs/quant-ph/0608095)

**Novelty conclusion:** I did not locate an earlier proof of the same unrestricted three-copy qutrit endpoint in this bounded search. That supports a serious originality claim, but is not an exhaustive priority certificate. The 2026 comparisons above are preprints whose scope I checked; their full proofs were not independently certified in this review.

## 7. Minor comments

**Resource reporting needs a platform correction.** The verifier writes `ru_maxrss` into a field named `peak_rss_bytes` without conversion. On Linux the native value is in KiB, so that field is mislabelled there. This is separate from the supplied supervisor’s aggregate-memory receipt and has no mathematical consequence. [man7.org](https://man7.org/linux/man-pages/man2/getrusage.2.html)

**Clarify historical status metadata.** The frozen certificate contains `portable_replay_status: PENDING_FRESH_REPLAY`, while the completed replay records establish the later pass. The README explains the distinction, so this is not contradictory evidence. A construction-specific field name in the next release would nevertheless reduce confusion. Do not silently modify the already frozen archive.

**Improve the publication-facing package.** Add a clear article date/version, an abstract and bibliography, and explicit reuse terms. Document the Unix-oriented assumptions associated with `resource`, `SIGALRM` and hard links. Keep the earlier partial-results material clearly labelled as historical input.

The instruction to use a new empty output directory for an independent replay is good and should remain prominent. Resuming saved receipts is useful operationally, but is not the same exercise as recomputing them.

## 8. Prioritised revision plan

**Must:** Integrate the literature and contribution statement into the main article. Separate verification of the final certificate from optional historical correction claims. Retain the exact finite-copy claim boundary and the complete complex rank-two quantifier.

**Should:** Present the proof as a short sequence of named results: physical lifting, preservation under averaging, valid positive maps, complete source coordinates, and the finite certificate identity. Add the exact calibration examples and a clearly scoped independent verification record. Correct resource-unit reporting and release metadata.

**Could:** State the stronger parameter conclusion that already follows from the endpoint theorem. No new certificate is required:

\[
\boxed{\rho(3,\alpha)\text{ is three-copy distillable}
\iff \alpha<-\frac12.}
\]

Here is the derivation. The three-copy endpoint theorem implies the one- and two-copy statements by appending product vectors \(|0\rangle_A|1\rangle_B\), each with expectation one. For \(-1/2\leq\alpha\leq0\), put \(s=-2\alpha\). Then

\[
W_\alpha=I+\alpha|\Omega\rangle\langle\Omega|
=sW_{-1/2}+(1-s)I.
\]

Expanding the third tensor power gives nonnegative coefficients. Resolving the identity factors in product bases reduces each term to an endpoint expectation on a vector whose Schmidt rank remains at most two. For \(\alpha\geq0\), \(W_\alpha\) is positive semidefinite.

Conversely, use \(|00\rangle+|11\rangle\) on one copy and \(|01\rangle\) on each other copy. The unnormalised expectation is \(2+4\alpha<0\) when \(\alpha<-1/2\).

This is a **reviewer-derived corollary**, supplied with a complete argument in the audit bundle. It does not extend the result to four copies or arbitrary local dimensions.

## 9. Editorial recommendation

**Major revisions, with a positive assessment of the central mathematics.**

The required revisions concern the scientific presentation and the precise scope of reproducibility assertions. I am not recommending that the authors replace the theorem, search for more favourable numerical examples, or solve the all-copy problem before publication.

My confidence is strongest in the finite arithmetic and the physical meaning of the certificate. It is lower in exhaustive priority and the eventual breadth of methodological influence. A specialist referee should still scrutinise the physical lifting, coupled-map interpretation and representation-theoretic reduction before a major journal publication.

The next useful independent check would target those mathematical interfaces rather than simply repeat the same receipt-reading procedure.

**On the evidence reviewed, the work should be treated as a computer-assisted proof of the specified endpoint, subject to the ordinary possibility of an undiscovered error—not as an unresolved numerical conjecture.**

## 10. Prospective REF assessment

The most appropriate primary comparison is **UoA 10: Mathematical Sciences**, with **UoA 9: Physics** also plausible.

I checked the current REF 2029 Contributions to Knowledge and Understanding guidance. For star-level interpretation, I use the published Main Panel B definitions from REF 2021, rather than presenting them as newly finalised REF 2029 panel-specific criteria. The current framework recognises varied output forms, including software, datasets and multi-component outputs. [REF 2029](https://2029.ref.ac.uk/guidance/section-4-contributions-to-knowledge-and-understanding-cku-guidance/)

Main Panel B’s published four-star descriptors include work at the research frontier, substantial novelty and potential major influence; the three-star descriptors concern important internationally excellent contributions falling short of that highest level. [REF 2021](https://2021.ref.ac.uk/media/1450/ref-2019_02-panel-criteria-and-working-methods.pdf)

| Dimension | Provisional rating | Assessment |
|---|---:|---|
| **Originality** | **4\*** | The unrestricted three-copy qutrit certificate appears to establish a genuinely new endpoint result. This remains contingent on the bounded novelty assessment. |
| **Significance** | **4\*** | It closes a strategically important finite-copy case and supplies an explicit, inspectable certificate rather than another restricted-sector result. Influence beyond this case is promising but not yet demonstrated. |
| **Rigour** | **4\*** | The main argument is explicit, the complex quantifier is preserved, and the decisive finite certificate has passed alternative exact verification. Ancillary construction claims need narrower presentation or additional evidence. |
| **Overall** | **4\* candidate** | A holistic judgement of a potentially world-leading mathematical contribution, not an arithmetic average of component scores. |

**Confidence: moderate.** A conservative panel could place it at 3\* if it judged the contribution’s reach too narrow or found its novelty insufficiently distinguished. I would not impose that ceiling merely because the theorem concerns three copies and qutrits: a result can lead its research area without resolving every surrounding conjecture.

Nor would I downgrade it simply for being computer-assisted, lacking a Lean formalisation, or being supplied as a multi-component package.

This is a quality assessment, not an official REF outcome or a determination of institutional submission eligibility.

## 11. References

### Quantum-information literature

Bharti, K., Gajjala, R., & Haug, T. (2026). *Two-copy nondistillability of Werner states: Sharp partial-trace inequalities and finite-copy extensions* [Preprint]. arXiv:2607.24479v1. DOI: 10.48550/arXiv.2607.24479. [arXiv](https://arxiv.org/abs/2607.24479)

Costa Rico, P. (2025). New partial trace inequalities and distillability of Werner states. *Letters in Mathematical Physics, 115*, Article 47. DOI: 10.1007/s11005-025-01935-y. [Springer](https://link.springer.com/article/10.1007/s11005-025-01935-y)

Fu, J., Gao, L., & Park, S.-J. (2026). *A solution to 2-copy distillability of Werner states* [Preprint]. arXiv:2607.21367v2. DOI: 10.48550/arXiv.2607.21367. [arXiv](https://arxiv.org/abs/2607.21367)

Vianna, R. O., & Doherty, A. C. (2006). Study of the distillability of Werner states using entanglement witnesses and robust semidefinite programs. *Physical Review A, 74*, 052306. DOI: 10.1103/PhysRevA.74.052306. **Abstract and bibliographic record consulted.** [arXiv](https://arxiv.org/abs/quant-ph/0608095)

Wu, T., & Zou, Q. (2026). *Sharp Plücker geometry for three-copy Werner distillation* [Preprint]. arXiv:2608.02647v1. DOI: 10.48550/arXiv.2608.02647. [arXiv](https://arxiv.org/abs/2608.02647)

### Assessment and technical guidance

Research Excellence Framework. (2019). *Panel criteria and working methods* (REF 2019/02), Main Panel B, paragraph 202, pp. 36–37. [REF 2021](https://2021.ref.ac.uk/media/1450/ref-2019_02-panel-criteria-and-working-methods.pdf)

Research Excellence Framework. (2026, January 22). *Section 4—Contributions to Knowledge and Understanding (CKU) guidance*. [REF 2029](https://2029.ref.ac.uk/guidance/section-4-contributions-to-knowledge-and-understanding-cku-guidance/)

Linux man-pages project. (2026). *getrusage(2)*. Linux man-pages 6.19. [man7.org](https://man7.org/linux/man-pages/man2/getrusage.2.html)