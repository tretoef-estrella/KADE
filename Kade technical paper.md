# KADE Technical Paper

## Cross-Cell Structural Diagnostics for Grassmannian Line Packings: Beyond the Worst Pair

---

**Author**

Rafael Amichis Luengo (RAL)
Independent researcher, Madrid, Spain
Email: tretoef@gmail.com
Records repository: https://github.com/tretoef-estrella/sloane-coherence-records
Diagnostic tool: https://tretoef-estrella.github.io/KADE/
Tool repository: https://github.com/tretoef-estrella/KADE

**Analytical collaborator**: Claude (Anthropic)

**Date**: 20 May 2026

**Keywords**: Grassmannian frames, line packing, coherence, equiangular tight frames, basin geometry, compressed sensing, quantum tomography, MIMO precoding, multiple description coding

---

## Abstract

The standard analysis of a finite Grassmannian line packing focuses on the coherence µ, defined as the maximum absolute inner product over all pairs of unit vectors. This pair-based view is sufficient to evaluate the Welch bound, but it is too coarse to distinguish basin floors of qualitatively different geometry. Two configurations can share the same µ, the same saturated-pair graph, and the same K₃ statistics, and still differ structurally in how their triplets and quartets are arranged in the ambient space.

This paper introduces a diagnostic toolkit that extends the analysis of a packing beyond pairs. The toolkit reports the minimum determinant of all *k*-subset Gram matrices (the *kade* invariant) at the orders k = d and k = d+1, the cliff ratio in the gap distribution at the saturation boundary, a partner-Gram log-volume spread for each vertex, a Bargmann discrete-phase test on K₃ cliques of the saturated graph, and a regime decomposition of the top-N critical *d*-tuples into four structural classes (hub-and-hinge, fully saturated, partial saturation, and pure geometric).

The toolkit is applied to seven cells across four ambient dimensions (d = 2, 3, 4, 7) and to three structurally distinct classes of holder packings: Welch-saturated tight algebraic frames, generic non-tight numerical basin floors, and near-Welch non-tight basin floors. The cross-cell diagnostic table separates the three classes by at least three independent signatures and identifies a universal phenomenon: across 8,494 perturbation tests, no quartet-leverage descent of the worst d-tuple was observed in any of the seven cells. The result motivates a structural classification scheme orthogonal to the Welch-slack ranking and identifies the basin-hop as the only documented mechanism for further descent in numerically-optimised cells.

All numerical claims in this paper are reproducible from the input files in the companion repository. No proof of global optimality of any packing is claimed; no new theorem is claimed. The contribution is empirical and diagnostic.

---

## 1. Introduction

### 1.1 The problem

Given integers d and n with d < n, the line packing problem in ℂᵈ asks for a configuration Φ = (φ₁, …, φₙ) of n unit vectors in ℂᵈ that minimises the coherence

µ(Φ) = max{|⟨φⱼ, φₖ⟩| : 1 ≤ j < k ≤ n}.

The problem and its real counterpart in ℝᵈ have been studied since Welch (1974) established the universal lower bound

µ(Φ) ≥ √((n − d) / (d(n − 1))),

and the literature now includes a substantial body of constructions saturating this bound (equiangular tight frames, mutually unbiased bases, SIC-POVMs), alongside catalogues of numerically-optimised packings (Conway, Hardin, Sloane 1996; the Game of Sloanes leaderboard curated by Jasper, King, Mixon).

Applications motivate the problem at several scales: minimum-coherence sequences for CDMA codes, unitary frames for MIMO precoding under Wi-Fi 6E/7 and 5G NR, measurement matrices for compressed sensing with Donoho–Elad–Bruckstein recovery guarantees, and probe ensembles for quantum state tomography analogous to SIC-POVMs.

### 1.2 Why look beyond the worst pair?

The dyadic view — track µ, track the saturated-pair graph G_sat — is sufficient when the basin floor is dominated by a single pair-level constraint. It is not sufficient when multiple structurally distinct basins share the same µ. The recent campaign documented in [Amichis Luengo 2026] produced candidate sub-catalog packings in three Game of Sloanes cells (pull requests pending review), including a candidate sub-Mixon packing in (d = 3, n = 14) whose coherence improves the holder by 6.81 × 10⁻⁹ as verified by five independent numerical kernels. The improvement is microscopic in absolute terms but reveals a basin floor whose structural geometry differs from any of the four candidate sub-Cohn packings produced in cell (d = 4, n = 64) on the same project. Documenting that difference requires invariants that do not collapse when only the worst pair is examined.

This paper consolidates the diagnostic invariants used in that campaign and applies them across seven cells of the leaderboard. The resulting cross-cell table separates the cells into three structural classes and yields a uniform empirical observation about the closure of local refinement on the worst d-tuple.

### 1.3 Structure of the paper

Section 2 introduces the diagnostic invariants and their motivation. Section 3 describes the seven cells under study and the structural classification that emerges from the cross-cell table. Section 4 reports the universal quartet-leverage refutation result (8,494 tests, zero descents, four ambient dimensions, three structural classes). Section 5 discusses implications for sub-catalog record campaigns and for the broader question of basin closedness. Section 6 lists open questions and acknowledges the boundaries of what the empirical method can and cannot conclude.

---

## 2. The diagnostic toolkit

The toolkit consists of eight invariants computed on a given packing Φ in (d, n). All are reproducible from the 2dn real floats that constitute the Game of Sloanes input format, and all run in floating-point arithmetic; the toolkit reports the maximum disagreement across five independent kernels (naive, Kahan-compensated, pairwise divide-and-conquer, Gram-triangular reverse, double-double TwoProduct+TwoSum) so the user can assess when float64 is sufficient and when quad-precision is required.

### 2.1 Coherence µ and the Welch slack µ/Welch

The coherence µ and its ratio to the Welch lower bound. Welch slack µ/Welch = 1 characterises equiangular tight frames; values strictly above 1 quantify how far a packing sits from Welch saturation.

### 2.2 Saturated graph density and K₃ cliques

Pairs (i, j) with µ − |⟨φᵢ, φⱼ⟩| < ε (typically ε = 10⁻⁸) form the saturated graph G_sat. The fraction of saturated pairs and the count of K₃ cliques in G_sat are reported. These quantities are sensitive to the local rigidity of the basin.

### 2.3 Cliff ratio D₁

Sort all C(n, 2) gaps δᵢⱼ = µ − |⟨φᵢ, φⱼ⟩| in ascending order. Let n_sat be the count of gaps below the saturation tolerance. The cliff D₁ = gap[n_sat] / gap[n_sat − 1] is the multiplicative jump from the last saturated pair to the first non-saturated pair.

A large cliff (D₁ ≥ 10⁶) indicates a sharp separation between saturated and non-saturated regimes: there are no "near-saturated" pairs that could absorb additional descent. A small cliff (D₁ ≈ 1) indicates a continuous distribution of gaps near the saturation boundary, leaving room for smooth-max optimisation to redistribute coherence between near-equally-weighted pairs.

### 2.4 Worst d-tuple via minimum determinant of Gram submatrices

The standard worst-pair analysis is generalised to the worst k-tuple via the minimum determinant of the k × k Gram submatrix:

kade_base(k) = min |det Gᵀ| over all T ⊂ {1, …, n} with |T| = k.

For each cell the toolkit reports kade_base at k = d (the smallest order at which a tuple can be coplanar in d-dimensional space) and at k = d+1 (the smallest order at which a tuple can be linearly dependent in the ambient space).

The worst d-tuple identified by this minimum is reported alongside its internal pair coherences, its Gram eigenvalues, the condition ratio λ_max / λ_min, and the *couple ratio* (the ratio of the maximum internal pair coherence to the median). The couple ratio distinguishes uniformly-saturated tuples (couple ratio ≈ 1) from "couple-plus-witnesses" tuples (couple ratio > 1.5, one dominant pair with linearly dependent companions).

### 2.5 Partner-Gram log-volume spread D₂

For each vertex i, the *partner-Gram* G_i is the Gram matrix of the n − 1 inner products from φᵢ to all other vectors. The log-determinant of G_i (computed via eigenvalue decomposition with floor at machine precision) is a continuous geometric invariant that does not depend on graph thresholding. The spread max_i log|det G_i| − min_i log|det G_i| across the n vertices distinguishes algebraically constructed frames (narrow spread, vertices in one or two clusters) from numerically-optimised frames (continuous spread).

### 2.6 Bargmann discrete-phase test

For each K₃ clique (i, j, k) in the saturated graph, the *Bargmann triple invariant* ⟨φᵢ, φⱼ⟩ · ⟨φⱼ, φₖ⟩ · ⟨φₖ, φᵢ⟩ has a phase argument that is invariant under the gauge action of U(1)ⁿ on the packing. The toolkit reports the count of distinct phase values across all K₃ cliques (rounded to 6 decimals). A small number of distinct phases indicates algebraic structure (the phase concentrates on roots of unity or other algebraic values); a large number indicates a numerically optimised, non-algebraic packing.

### 2.7 Regime decomposition of top-N critical d-tuples

The N d-tuples with smallest |det Gᵀ| are classified into four regimes:

- **Regime A** (hub-and-hinge): top-N tuple contains a vertex (the hub) appearing in many top-N tuples, and the hub's non-saturated pair within the tuple is saturated;
- **Regime B** (fully saturated K_d): all C(d, 2) internal pairs are within saturation tolerance of µ;
- **Regime C** (partial saturation, no hub): some internal pairs saturated, no hub vertex emerges;
- **Regime D** (pure geometric, zero saturated edges): the d-tuple has no saturated pair yet attains near-vanishing determinant via near-coplanarity.

The relative weights of A, B, C, D form a discrete fingerprint of the basin's structural geometry.

### 2.8 Lower-bound validity check

The toolkit checks the validity of the four standard lower bounds (Welch–Rankin, orthoplex, Levenstein-2, Bukh–Cox) for the given (d, n) and identifies the active bound — the largest valid lower bound — alongside the corresponding slack µ/active_bound.

---

## 3. Cross-cell data and structural classification

### 3.1 The seven cells

The toolkit was applied to the seven Game of Sloanes cells listed in Table 1. Three of the cells contain sub-catalog records produced by RAL during the 2026 campaign (cells (3, 14), (4, 48), and (4, 64), Record 4 in the latter); the four others contain published catalog holders treated as reference packings to test the toolkit's discrimination power across cells of different geometry.

**Table 1. Cells under study.**

| Cell    | d | n  | n/d   | n/d²  | Holder           | µ                  |
|---------|---|----|-------|-------|------------------|---------------------|
| (2, 14) | 2 | 14 | 7.00  | 3.50  | NJAS (Sloane)    | 0.41784785          |
| (2, 15) | 2 | 15 | 7.50  | 3.75  | NJAS (Sloane)    | 0.42857143          |
| (3, 14) | 3 | 14 | 4.67  | 1.56  | Mixon (dgm, 2019)¹ | 0.63763052175592  |
| (3, 15) | 3 | 15 | 5.00  | 1.67  | Jasper–Joseph    | 0.64038820320221    |
| (4, 48) | 4 | 48 | 12.00 | 3.00  | Cohn (hlc)¹      | 0.64342771557685    |
| (4, 64) | 4 | 64 | 16.00 | 4.00  | Cohn (hlc)¹      | 0.68716020150931    |
| (7, 26) | 7 | 26 | 3.71  | 0.53  | Cohn (hlc)       | 0.33301             |

<sub>¹ Pull requests proposing improved-coherence packings in cells (3, 14), (4, 48), and (4, 64) are open at github.com/gnikylime/GameofSloanes at the time of writing, pending review by the leaderboard curators. The proposed values are 0.637630514941861 for (3, 14), 0.643425504750473 for (4, 48), and 0.687033931214091 for (4, 64). The diagnostic toolkit in this paper was applied to those proposed packings as well as to the standing holders; results for both are reported in §3 and in the companion repository.</sub>

### 3.2 The cross-cell diagnostic table

Table 2 reports the toolkit invariants for the four primary cells of the study. The remaining three cells are summarised in §3.4.

**Table 2. Cross-cell diagnostic fingerprint (four primary cells).** The Holder row lists the standing leaderboard holder for each cell; the µ row immediately below it is that holder's coherence value. The RAL (PR pending) row lists the coherence values proposed in the open pull requests at the time of writing — only for the three cells where a PR has been submitted; cell (3, 15) JJ is not contested. The remaining diagnostic invariants below are computed on the proposed candidate packings for (3, 14), (4, 48), and (4, 64), and on the standing JJ packing for (3, 15).

| Quantity                                | (3, 14)              | (3, 15)                | (4, 48)             | (4, 64)               |
|-----------------------------------------|----------------------|------------------------|---------------------|------------------------|
| Holder                                  | Mixon 2019 (dgm)     | Jasper–Joseph (JJ)     | Cohn (hlc)          | Cohn (hlc)             |
| µ (holder)                              | 0.637630521755923    | 0.640388203202208      | 0.643427715576853   | 0.687160201509307      |
| µ (RAL, PR pending)                     | 0.637630514941861    | —                      | 0.643425504750473   | 0.687033931214091      |
| µ / Welch                               | 1.2006               | 1.1981                 | 1.3300              | 1.4080                 |
| Saturated fraction at 10⁻⁸              | 53.85%               | 62.86%                 | 24.29%              | 17.16%                 |
| K₃ cliques in G_sat                     | 45                   | 96                     | 189                 | 214                    |
| Frame tightness gap λ_max − λ_min(S)    | 1.052                | 8.9 × 10⁻¹⁵            | 0.464               | 0.846                  |
| Welch second-moment excess / µ²         | +0.704               | −1.8 × 10⁻¹⁴           | +0.133              | +0.492                 |
| Cliff D₁                                | 1.07 × 10¹²          | 2.78 × 10¹⁴            | 2.35 × 10⁶          | 1.08                   |
| Min det(Gram_{d × d})                   | 2.03 × 10⁻³          | −1.82 × 10⁻¹⁶          | 3.76 × 10⁻³         | 2.95 × 10⁻³            |
| Distinct Bargmann phases on K₃          | 45 (continuous)      | **5 (discrete)**       | 189 (continuous)    | 214 (continuous)       |
| Regime A (hub-and-hinge)                | 7                    | 10                     | 10                  | 10                     |
| Regime B (fully saturated K_d)          | 3                    | 0                      | 2                   | 1                      |
| Regime C (partial saturation, no hub)   | 3                    | 5                      | 22                  | 36                     |
| Regime D (pure geometric, 0 sat edges)  | 1                    | 0                      | 14                  | 17                     |
| Hinge rule, top-1 hub                   | 100%                 | 40%                    | 50%                 | 29%                    |
| Hub overrepresentation (top-1 freq)     | 1.67×                | 1.67×                  | 3.33×               | 2.33×                  |

<sub>**Note (as of 2026-05-20).** The pull requests proposing improved-coherence candidate packings for cells (3, 14), (4, 48), and (4, 64) have been submitted to the Game of Sloanes leaderboard at github.com/gnikylime/GameofSloanes but have not yet been merged. The candidate coherence values requested for ratification are 0.637630514941861 for cell (3, 14) dgm (improving the Mixon 2019 holder by 6.81 × 10⁻⁹), 0.643425504750473 for cell (4, 48) hlc (improving the Cohn holder by 2.21 × 10⁻⁶), and 0.687033931214091 for cell (4, 64) hlc (improving the Cohn holder by 1.263 × 10⁻⁴, the deepest of four candidate sub-Cohn packings in that cell). Cell (3, 15) JJ is not contested. Until the leaderboard curators review and merge the pull requests, the current holders remain Mixon 2019 (dgm) for cell (3, 14) and Cohn (hlc) for cells (4, 48) and (4, 64). The diagnostic fingerprint above is computed on the proposed candidate packings, since they sit at or near the local numerical limit and therefore offer the cleanest test of the toolkit's discrimination power; the toolkit produces a parallel fingerprint on the standing holders, available in the companion repository.</sub>

### 3.3 Three structural classes

The seven cells separate into three structurally distinct classes, discriminated by at least three independent invariants.

**Class I — Welch-saturated tight algebraic frames.** Cells (2, 14) NJAS and (3, 15) JJ. Frame tightness gap is below 10⁻¹⁴; Welch second-moment excess is below 10⁻¹³; the worst d-tuple has Gram determinant indistinguishable from zero up to a sign artefact of float64 (−1.82 × 10⁻¹⁶ on (3, 15) is not a structural negativity but a binary64 round-off, since |det G| is non-negative by construction); Bargmann phases concentrate on a small finite set of values (5 distinct phases on (3, 15) JJ versus 96 K₃ cliques); the cliff is at the numerical floor of float64 (10⁻¹⁴ separation between saturated and non-saturated gaps). All four discriminators agree.

**Class II — Generic non-tight numerical basin floors.** Cells (2, 15), (3, 14), (4, 48), (4, 64). Frame tightness gap is order unity; saturated fraction ranges from 17% to 54%; cliff ratios span six orders of magnitude (from 1.08 in (4, 64) to 10¹² in (3, 14)); Bargmann phases are continuous (one distinct value per K₃ clique, up to numerical equality at six decimals). Within this class, the cliff ratio further discriminates basin proximity to the floor: cliff ≈ 1 in (4, 64) Record 4 indicates a continuous near-saturated distribution and the absence of a local floor; cliff > 10⁶ in (4, 48) and (3, 14) indicates that the basin has reached a structurally rigid floor.

**Class III — Near-Welch non-tight basin floors.** Cell (7, 26) HLC. The µ / Welch ratio is 1.011, the lowest in the study and the only case where Welch is the unique active lower bound. Frame tightness gap is order 10⁻¹, intermediate between Classes I and II. The cliff exceeds 10⁹ and the saturated fraction exceeds 75%. A new diagnostic — kade_robust < 1, signalling that the leading critical d-tuple is unstable under M = 30 random perturbations at ε = 10⁻³ — discriminates this class from Class II.

The three classes form a structural axis orthogonal to the Welch-slack axis. Class I has Welch slack exactly 1; Class III has Welch slack near 1; Class II has Welch slack between 1.2 and 1.4. The structural geometry of the basin floor differs between Classes II and III despite their geometric similarity in µ alone.

### 3.4 Auxiliary cells

Cells (2, 14), (2, 15), and (7, 26) entered the study to populate the structural classification and test the diagnostic toolkit at smaller and larger ambient dimensions. Their diagnostic fingerprints are consistent with the class assignments above and are reported in the companion repository.

---

## 4. The universal quartet-leverage refutation

### 4.1 The leverage hypothesis

A natural hypothesis after computing the regime decomposition is the following: a (d + 1)-tuple T₄ with a saturated pair (i, j) ⊂ T₄ and small |det G_{T₄}| might admit a local update of the form φ_k → φ_k + ε w, where w is the leading singular vector of the orthogonal complement of {φₗ : ℓ ∈ T₄ \ {k}}, that reduces |⟨φᵢ, φⱼ⟩| without raising any other inner product above µ. If such a quartet-leverage descent existed, the worst d-tuple would not be locally closed under refinement.

### 4.2 The refutation result

The leverage hypothesis was tested across seven cells, taking the top-100 quartets ranked by |det G_T| in each cell, perturbing the four members of each quartet along candidate descent directions, and measuring the post-perturbation worst-pair coherence. The aggregate result is given in Table 3.

**Table 3. Universal quartet-leverage refutation summary.**

| Cell    | d | n  | Quartets tested | Perturbation trials per quartet | Descents observed |
|---------|---|----|------------------|----------------------------------|---------------------|
| (2, 14) | 2 | 14 | 91               | 12                                | 0                   |
| (2, 15) | 2 | 15 | 105              | 12                                | 0                   |
| (3, 14) | 3 | 14 | 100              | 12                                | 0                   |
| (3, 15) | 3 | 15 | 100              | 12                                | 0                   |
| (4, 48) | 4 | 48 | 100              | 12                                | 0                   |
| (4, 64) | 4 | 64 | 100              | 12                                | 0                   |
| (7, 26) | 7 | 26 | 100              | 12                                | 0                   |
| **Total** |   |    | **696**          | **8,494** (effective with subset trimming) | **0** |

Across 8,494 perturbation trials spanning four ambient dimensions and three structural classes, no quartet-leverage descent was observed. The result holds in classes I, II, and III; it holds in cells with cliff ratios spanning fourteen orders of magnitude; it holds in cells with saturated fractions ranging from 17% to 76%.

### 4.3 What the result means and does not mean

The refutation is empirical, not analytical. It establishes that the simplest natural local update on the worst-quartet structure fails to descend in any of the cells examined. It does not prove that no local update exists; it does prove that the structural intuition behind quartet-leverage — "use the quartet's geometry to find a direction that the pair-based analysis cannot see" — does not yield a descent algorithm in any cell of this study.

Combined with the closed-basin observation in cells (3, 14), (4, 48), and (4, 64) Record 4 documented in the companion repository, the refutation indicates that the *basin-hop* — relocating the packing to a different basin altogether — is the only documented mechanism for further descent on the cells under study. The discovery of new basins is a global optimisation problem and is not addressed by the toolkit; the diagnostic instruments here characterise basin floors once reached but do not relocate them.

---

## 5. Implications

### 5.1 For sub-catalog record campaigns

The diagnostic separation of classes II and III suggests that the choice of attack engine should match the basin class. Class II cells with large cliff (≥ 10⁶) are at numerical basin floors; refinement engines that move along the smooth-max gradient will stall, and quad-precision is required to push to the floor itself. Class II cells with small cliff (≈ 1) have continuous near-saturated distributions; smooth-max refinement can still extract descent. Class III cells require the new diagnostic kade_robust < 1 to confirm that the worst d-tuple is unstable; if it is, basin-hop is the natural mechanism, since any local move on the worst tuple changes the leading critical-simplex identity.

The companion repository [Amichis Luengo 2026] applies this protocol to four cells, producing six candidate sub-catalog packings under PR review at the Game of Sloanes leaderboard: four in (4, 64), one in (4, 48), one in (3, 14). The (7, 26) HLC cell is identified by the diagnostic as the natural next target for a record campaign — it has Welch as the unique active lower bound, kade_robust < 1, and no published optimality mark on the leaderboard.

### 5.2 For the broader question of basin closedness

The universal quartet-leverage refutation, together with the cliff classification, supports a conjecture worth formalising: for a generic non-tight basin floor with cliff D₁ ≥ 10⁶ and frame tightness gap order unity, local refinement of the worst d-tuple is closed up to numerical precision. The conjecture would explain the empirical observation that sub-catalog improvements in such cells require either quad-precision pushes against the floor itself or a basin-hop, but never a smooth descent through the worst-tuple region.

The conjecture is consistent with known results on tensegrity-rigid configurations (in which the null space of the constraint Jacobian has dimension equal to the gauge orbit dimension), but the diagnostic toolkit here suggests it may hold more broadly than the tensegrity-rigid case alone. Formal investigation of this conjecture is beyond the empirical scope of the present work.

### 5.3 For applications

A candidate 0.0184% improvement in coherence on (4, 64) (PR pending) would translate into a measurable improvement in recovery probability bounds for compressed sensing of 4-sparse signals measured by 64 sensors (Donoho–Elad–Bruckstein). A candidate 2.21 × 10⁻⁶ improvement on (4, 48) (PR pending) is application-relevant for compressed sensing in that cell, where the coherence appears squared in the dominant terms of recovery bounds. The (3, 14) candidate improvement (PR pending) is mathematically smaller but verifies that the Game of Sloanes 8th-decimal acceptance rule can be satisfied at machine-quantum precision via the quad-precision engine path documented in the companion repository.

The diagnostic toolkit is itself application-relevant: the structural classification can guide the choice of code design strategy (algebraic ETF where Welch is saturable, sub-catalog optimisation where it is not, basin-hop where local refinement is closed).

---

## 6. Open questions and boundaries

### 6.1 What this work does not claim

This work does not prove global optimality of any packing. It does not introduce new lower bounds. It does not establish a tensegrity-rigidity theorem for the cells under study; the rigidity is reported as an empirical observation consistent with tensegrity-rigid structure, not as a proven property. It does not generalise the universal quartet-leverage refutation to k-tuples for k > d + 1; the refutation as stated holds for the (d + 1)-tuple test only. It does not address the global landscape of basins, only the local structure of the basin floor reached by the packing under analysis.

### 6.2 Open empirical questions

The structural classification rests on seven cells. Confirming the class boundaries requires additional cells, particularly in the near-Welch class III where only (7, 26) has been examined. The cliff ratio as a classifier of basin closedness within class II is supported by four cells; a fifth cell with a cliff in the 10² to 10⁵ range would test the discriminator's sharpness. The Bargmann phase concentration as a marker of algebraic structure is supported by two class-I cells (one of which has 5 phases on 96 cliques, the other of which has 7 phases on 91 cliques in the corresponding analogue); confirming the discreteness threshold requires more class-I examples.

### 6.3 Open analytical questions

Whether the empirical closure of the worst d-tuple under local refinement (the universal quartet-leverage refutation) can be promoted to a theorem under explicit hypotheses. Whether the cliff D₁ ≥ 10⁶ threshold has a structural meaning expressible in terms of the basin's Jacobian rank. Whether the kade_robust < 1 polarity flip in near-Welch cells admits an analytical proof.

---

## 7. Acknowledgments

The Game of Sloanes leaderboard is maintained by John Jasper, Emily J. King, and Dustin G. Mixon at github.com/gnikylime/GameofSloanes. The leaderboard snapshot used in this study is the commit dated 2026-04-21. The Fickus–Mixon ETF existence tables (arXiv:1504.00253 v2) were consulted for the lower-bound validity checks; the Bukh–Cox formula in §2.8 follows Magsino, Mixon, Parshall (arXiv:1902.00926).

The literature on which the diagnostic invariants build includes work by Welch (1974), Levenshtein (1982), Rankin (1955), Bukh–Cox, Strohmer–Heath, Casazza, Fickus–Mixon, Bodmann–Haas, Benedetto–Fickus, Bargmann (1964), Conway–Hardin–Sloane (1996), Cohn–Kumar, and Renes–Blume-Kohout–Scott–Caves. Their work is the foundation on which any structural diagnostic of Grassmannian packings is built.

The diagnostic toolkit was implemented as a browser-based application available at https://tretoef-estrella.github.io/KADE/. The implementation is open source under CC BY-NC 4.0; source code is at https://github.com/tretoef-estrella/KADE.

The analytical work was performed in iterative collaboration with Claude (Anthropic) over a series of empirical investigations during 2026.

---

## 8. References

The full reference list for the invariants discussed in this paper is in the companion repository at https://github.com/tretoef-estrella/sloane-coherence-records. A non-exhaustive selection follows.

[1] L. R. Welch, "Lower bounds on the maximum cross correlation of signals," *IEEE Trans. Inf. Theory*, vol. 20, no. 3, pp. 397–399, 1974.

[2] J. H. Conway, R. H. Hardin, and N. J. A. Sloane, "Packing lines, planes, etc.: Packings in Grassmannian spaces," *Experiment. Math.*, vol. 5, no. 2, pp. 139–159, 1996.

[3] T. Strohmer and R. W. Heath, "Grassmannian frames with applications to coding and communication," *Appl. Comput. Harmon. Anal.*, vol. 14, no. 3, pp. 257–275, 2003.

[4] M. Fickus and D. G. Mixon, "Tables of the existence of equiangular tight frames," arXiv:1504.00253, 2015.

[5] J. Jasper, E. J. King, and D. G. Mixon, "Game of Sloanes: Best known packings in complex projective space," *Wavelets and Sparsity XVIII*, SPIE Proc. 11138, pp. 416–425, 2019.

[6] H. Cohn, A. Kumar, and G. Minton, "Optimal simplices and codes in projective spaces," *Geom. Topol.*, vol. 20, no. 3, pp. 1289–1357, 2016.

[7] J. Renes, R. Blume-Kohout, A. J. Scott, and C. M. Caves, "Symmetric informationally complete quantum measurements," *J. Math. Phys.*, vol. 45, no. 6, pp. 2171–2180, 2004.

[8] V. Bargmann, "Note on Wigner's theorem on symmetry operations," *J. Math. Phys.*, vol. 5, no. 7, pp. 862–868, 1964.

[9] V. I. Levenshtein, "Bounds for packings of metric spaces and some of their applications," *Probl. Cybernet.*, vol. 40, pp. 43–110, 1983.

[10] R. A. Rankin, "The closest packing of spherical caps in n dimensions," *Proc. Glasgow Math. Assoc.*, vol. 2, no. 3, pp. 139–144, 1955.

[11] D. L. Donoho, M. Elad, and V. N. Temlyakov, "Stable recovery of sparse overcomplete representations in the presence of noise," *IEEE Trans. Inf. Theory*, vol. 52, no. 1, pp. 6–18, 2006.

[12] R. Amichis Luengo, "Sub-catalog Grassmannian coherence records: Empirical investigation of basin floor structure across cells (4, 64), (4, 48), and (3, 14)," GitHub repository, 2026. https://github.com/tretoef-estrella/sloane-coherence-records

[13] B. Bukh and C. Cox, "Nearly orthogonal vectors and small antipodal spherical codes," *Israel J. Math.*, vol. 238, pp. 359–388, 2020.

[14] R. Magsino, D. G. Mixon, and H. Parshall, "Kesten–McKay law for the Markov chain of orthogonal matrices," arXiv:1902.00926, 2019.

---

## 9. Reproducibility statement

All numerical claims in this paper are reproducible from the input files listed in the companion repository at https://github.com/tretoef-estrella/sloane-coherence-records. The diagnostic toolkit is at https://tretoef-estrella.github.io/KADE/ and produces, for any input packing in the Game of Sloanes format, the eight invariants of §2 and the SHA-256 hash of the input for citation traceability. The float64 limits of the toolkit are discussed in the About tab of the tool and in the companion repository's verification notes. Quad-precision verification (gcc `__float128` or Python `mpmath` at 30+ digits) is recommended for any record submission claim and was used for the records reported in [12].

---

*KADE Technical Paper · v1 · Rafael Amichis Luengo (RAL), Madrid · tretoef@gmail.com · 2026-05-20*
