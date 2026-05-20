# KADE Packing Diagnostic

**A browser-based structural diagnostic for complex Grassmannian line packings.**

**Live tool**: https://tretoef-estrella.github.io/KADE/

Drag a `.txt` packing in. Read its structural fingerprint in under three seconds. No install. No upload. Nothing leaves your browser.

---

## Why this tool exists

Every Grassmannian packing has a coherence value μ — the worst inner-product magnitude between any two of its unit vectors. The Game of Sloanes leaderboard ranks packings on this single number. **μ alone, however, is not the same packing.** Two configurations with identical μ can sit in basins of completely different geometry, behave differently under refinement, and respond differently to the lower bounds that constrain them.

KADE computes the rest. Given a packing in standard Game of Sloanes format (a flat list of 2·*d*·*n* real floats, real block then imaginary block), it returns in a single screen:

- **μ across five independent float64 kernels** with the max disagreement reported, so you know exactly when you are looking at structural truth and when you are looking at floating-point noise.
- **The Welch slack μ/Welch** — how far above the universal lower bound the packing actually sits.
- **The saturated-pair graph G_sat** and its cliff D₁ — the ratio between the smallest saturated gap and the largest non-saturated gap. A large cliff means the basin floor is dominated by a clean structure; a cliff near 1 means the basin is gravelly.
- **The worst d-tuple** — the d-element subset whose Gram matrix is closest to singular. Whether its near-singularity is hub-driven, fully saturated, or geometric, classified into four regimes (A / B / C / D).
- **Frame-operator tightness gap** λ_max − λ_min(S) — whether the packing is a tight frame or how far from one.
- **Bargmann triple-product phases** on K₃ cliques of G_sat — whether the basin sits on a small finite set of discrete phases (algebraic ETF signature) or on a continuous distribution.
- **Validity of four standard lower bounds** — Welch–Rankin, orthoplex, Levenstein-2, and notes on Bukh–Cox — telling you which bound is the active constraint for the cell.
- **Side-by-side comparison against the current Game of Sloanes holder** for the same cell, when the cell is on the leaderboard.

Everything runs locally. The only outbound traffic is an optional fetch to `raw.githubusercontent.com` for the holder-comparison feature.

---

## What you see in three seconds

Drop a packing in. You get something like this — these are the actual numbers for the Mixon (2019) holder of cell (3, 14), one of the most-studied cells on the leaderboard:

```
─── KADE fingerprint ────────────────────────────────────────────
Cell:                      (d=3, n=14) — C^3, 14 unit vectors
Coherence μ:               0.637630521755923
   5-kernel max disagreement: 8.3 × 10⁻¹⁶   (float64 is comfortable)

Lower bounds:
   Welch–Rankin (active):  0.531085...   ratio μ/Welch = 1.2006
   Orthoplex:              valid (μ above 1/√d)
   Levenstein-2:           0.644503...   (NOT saturated by this packing)

Saturated-pair graph G_sat:
   Saturated fraction:     53.85% of all C(14,2) = 91 pairs
   K₃ cliques in G_sat:    45 distinct triangles
   Cliff D₁:               1.07 × 10¹²   ← clean basin floor

Worst d-tuple (d=3):
   Min det(Gram_{3×3}):    2.03 × 10⁻³
   Regime classification:  Regime A (hub-and-hinge, 7 of 14 triples)
                           Regime B (fully saturated, 3)
                           Regime C (partial, 3) — Regime D (geometric, 1)
   Top-1 hub overrepresentation: 1.67× expected frequency

Bargmann phases on K₃:     45 distinct (continuous distribution)
Frame tightness gap:       1.052  (NOT a tight frame)
─────────────────────────────────────────────────────────────────
```

That is the fingerprint. The number on top, μ, is what the leaderboard ranks. Everything below it is structural information about *why* this packing sits where it sits. The same screen also offers a JSON export with SHA-256 hash of the input, for citation in downstream work.

---

## How to read the fingerprint

Five numbers do most of the work. Once you can read these, you can tell at a glance what kind of basin a packing lives in:

1. **5-kernel max disagreement.** Below ≈ 10⁻¹⁴ means float64 is reproducing μ stably and the value you see is structural. Above ≈ 10⁻¹⁰ means the cell has reached the limit of double-precision arithmetic and any record submission for it needs quad-precision verification. KADE is honest about this; it does not hide kernel disagreement under a single rounded number.

2. **μ / Welch.** The ratio between the packing's coherence and the Welch lower bound. Ratio = 1 means the packing saturates Welch and is therefore an equiangular tight frame (ETF) if it exists. Ratio between 1.2 and 1.5 is the typical territory of hard cells where Welch is the best known active bound. Ratio close to 1 with a non-tight frame operator means something interesting is happening — possibly a Levenstein-2 or orthoplex regime.

3. **Cliff D₁.** Large (10⁶ and above) means the basin floor is *structurally clean* — the saturated and non-saturated pair populations are clearly separated, the worst-pair indices are stable across verification kernels, and the packing is sitting on a well-defined local minimum. Small (near 1) means the floor is gravelly — many pairs are nearly tied and the geometry is essentially Welch-saturated.

4. **Distinct Bargmann phases on K₃ cliques.** A small discrete set (3, 5, 7 distinct phases) is the signature of an algebraic ETF or a structurally constructed frame, where the Bargmann invariant lives on a finite Galois-like orbit. A continuous spread of phases means the basin was reached by numerical optimisation and the configuration is empirical, not algebraic.

5. **Regime distribution of the worst d-tuples.** Whether the d-element subsets closest to singularity are hub-and-hinge (one vector shared by many near-singular triples), fully saturated (every pair maxed out), partial, or purely geometric. This is what makes two packings with identical μ live in different basins.

The full theory of these quantities — citations, derivations, edge cases, the cliff D₁ definition, the Bargmann triple-product computation — lives in the **About & Theory** tab inside KADE itself, and in the companion technical paper [`KADE_Technical_Paper.pdf`](KADE_Technical_Paper.pdf), available in this repository.

---

## Three tabs

**1. Calculator.** Paste or upload a packing in Game of Sloanes format. Auto-detection of (d, n) from the standard `dxn_xxx.txt` filename pattern. Full diagnostic fingerprint computed in under three seconds for cells with C(n, d) ≤ 700 000 d-tuples. Export as JSON with input SHA-256 for citation. Cells beyond the d-tuple budget receive the rest of the diagnostic and skip the worst-d-tuple section explicitly — the tool tells you, never silently.

**2. Screener.** A static snapshot of the Game of Sloanes leaderboard (commit dated 2026-04-21 upstream). Each row links to the exact holder file in the upstream repository. A **Check upstream** button performs a live fetch of the upstream README and reports whether the leaderboard date has changed since the snapshot — so you know when KADE's embedded view is stale and a manual refresh is worth doing.

**3. About & Theory.** Every quantity explained — what it measures, why it matters, the literature it draws on, and an honest list of what the tool does *not* claim. Designed to be read first by anyone who wants to understand the numbers before trusting them. Includes the precedents in Welch, Levenshtein, Rankin, Bukh–Cox, Strohmer–Heath, Casazza, Fickus–Mixon, Bodmann–Haas, Benedetto–Fickus, Bargmann, Conway–Hardin–Sloane, Cohn–Kumar, and Renes–Blume-Kohout–Scott–Caves.

---

## Why five kernels for μ

The Game of Sloanes 8th-decimal acceptance rule means that two packings can differ in coherence by ~10⁻⁹ and one of them is a record. At that scale, the way you compute μ matters: a naive `max(abs(G[i,j]))` over a float64 Gram matrix can drift in the 12th–14th decimal depending on accumulation order. For most cells this drift is invisible. For some cells — especially those whose basins are at the local numerical limit — it is decisive.

KADE computes μ five independent ways:

1. **Naive** — straightforward double-precision max-abs over the Gram matrix.
2. **Kahan-compensated** — Kahan summation in the inner product, drift-corrected per pair.
3. **Pairwise divide-and-conquer** — recursive halving of the sum, classical numerical-analysis stable.
4. **Gram-triangular with reverse accumulation** — Gram computed bottom-up, indices walked in reverse to expose accumulation asymmetry.
5. **Double-double TwoProduct + TwoSum** — roughly 31 extra bits of precision over float64, used as the in-browser quad-precision approximation.

The maximum disagreement among the five is reported alongside μ. **You see the noise, not the headline.** When the disagreement is comfortable (< 10⁻¹⁴) you can use the value with confidence; when it is uncomfortable (> 10⁻¹⁰) the cell has saturated float64 and you should trust only the double-double kernel — or, for a serious record submission, fall back on an external quad-precision verifier (gcc `__float128` or Python `mpmath` at ≥ 30 digits).

This discipline is what makes KADE usable for pre-submission sanity checks. The float64 disclaimer is the second sentence of the About tab, not buried at the bottom.

---

## Who this is useful for

- **Researchers in frame theory** who want a structural fingerprint of a packing in seconds without a Jupyter session.
- **Students** learning Grassmannian geometry who want to see what lives beyond μ.
- **Game of Sloanes contributors** preparing a submission and wanting a sanity check before opening the PR.
- **Practitioners** in compressed sensing, MIMO precoding, CDMA, quantum tomography, or multiple-description coding who use Grassmannian frames and want a fast structural reading of the configurations they work with.
- **Anyone curious about the geometry** of finite-dimensional packings, who wants a tool that does not hide what it is doing.

KADE does not replace mathematical reasoning, quad-precision verification, or peer review. It is one diagnostic instrument among others. It does not claim global optimality of any packing it analyses, and it does not claim that float64 is sufficient for an 8th-decimal record submission. What it offers is the convenience of computing a serious structural fingerprint in the browser, on any packing on your disk, in seconds, with full reporting of the numerical regime you are in.

---

## How KADE came to exist

KADE was not designed in the abstract. It emerged from an empirical campaign on the Game of Sloanes leaderboard during 2026, in which six candidate sub-catalog packings were produced across three cells — (4, 64) hlc, (4, 48) hlc, and (3, 14) dgm — and submitted as pull requests to the upstream leaderboard. Every quantity now reported in KADE earned its place by being the thing that distinguished one basin from another during that work, or by being the thing that revealed when a refinement strategy was working and when it was hitting a structural wall.

The records, the verification protocols, the engine codenames (Aquiles, Boagrius, Sanjuanbautista, Rasputin), and the basin-floor technical notes are all documented in the companion repository [sloane-coherence-records](https://github.com/tretoef-estrella/sloane-coherence-records). Readers who want to see KADE applied to non-trivial cases will find the full fingerprint of every record packing there, along with the equivalent fingerprint of the standing holder it improves upon — six worked examples, all reproducible byte-exact from the raw `.txt` files in that repository.

The companion technical paper [`KADE_Technical_Paper.pdf`](KADE_Technical_Paper.pdf) (in this repository) documents the toolkit at the level expected of a research preprint: definitions, cross-cell tables, the structural-class taxonomy, the universal quartet-leverage refutation that emerged from running 8 494 worst-quartet tests across seven cells with no observed descents, and the open questions the toolkit makes precise.

---

## Technical notes

- Single HTML file. No build step, no server, no external JavaScript libraries. Google Fonts CSS is loaded for typography; everything else is self-contained.
- Hermitian eigenvalues computed via Jacobi rotation on the 2d × 2d real-block expansion of the d × d complex Gram. Adequate for the small *d* range typical of the Game of Sloanes leaderboard (d ≤ 8).
- Worst d-tuple enumeration capped at C(n, d) ≤ 700 000 to keep the browser responsive. Cells beyond that budget show the rest of the fingerprint and skip the worst-d-tuple section, with an explicit message rather than a silent omission.
- Tested no-regression against several published holder packings from the Game of Sloanes leaderboard. The Mixon (3, 14) dgm, Jasper–Joseph (3, 15), and Cohn (4, 48) and (4, 64) hlc holders are part of the regression suite.

---

## Citation

If you use KADE in published work, please cite this repository:

```
Amichis Luengo, R. (2026). KADE Packing Diagnostic.
https://github.com/tretoef-estrella/KADE
```

The companion technical paper, when cited as a methodological reference for the diagnostic invariants themselves, is:

```
Amichis Luengo, R. (2026). KADE: A Structural Diagnostic for
Complex Grassmannian Line Packings. Technical paper, available at
https://github.com/tretoef-estrella/KADE/blob/main/KADE_Technical_Paper.pdf
```

---

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

You are free to share and adapt this work for **non-commercial** purposes (research, education, personal use), provided you give appropriate credit and link to the license. Commercial use requires explicit permission — contact `tretoef@gmail.com`.

---

## Author

Rafael Amichis Luengo (RAL), Madrid — `tretoef@gmail.com`

Independent investigator. The diagnostic framework was developed in iterative collaboration with Claude (Anthropic) during empirical cross-cell investigations of Grassmannian frame packings. The records those investigations produced are documented in the companion repository [sloane-coherence-records](https://github.com/tretoef-estrella/sloane-coherence-records).

---

## Acknowledgements

The Game of Sloanes leaderboard is curated by John Jasper, Emily J. King, and Dustin G. Mixon at [github.com/gnikylime/GameofSloanes](https://github.com/gnikylime/GameofSloanes) — the reference repository against which any tool of this kind is measured.

The literature cited inside the **About & Theory** tab and in the companion technical paper includes foundational work by Welch, Levenshtein, Rankin, Bukh–Cox, Strohmer–Heath, Casazza, Fickus–Mixon, Bodmann–Haas, Benedetto–Fickus, Bargmann, Conway–Hardin–Sloane, Cohn–Kumar, and Renes–Blume-Kohout–Scott–Caves. Their work is the foundation on which any tool of this kind is built.

Anthropic — for building Claude, the AI assistant that collaborated on the engineering, the structural-classification taxonomy, and the in-browser verification kernels.

---

*Last updated: 2026-05-20.*
