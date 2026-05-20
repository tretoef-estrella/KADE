# KADE Packing Diagnostic

A browser-based structural diagnostic tool for finite-dimensional complex
Grassmannian line packings.

**Live tool**: [https://tretoef-estrella.github.io/KADE/](https://tretoef-estrella.github.io/KADE/)

---

## What this is

Given a packing Φ of *n* unit vectors in ℂ<sup>d</sup> in the standard
[Game of Sloanes](https://github.com/gnikylime/GameofSloanes) format (a list
of 2·*d*·*n* real floats, real block then imaginary block), the tool computes
a compact structural fingerprint and presents it in a single page.

The fingerprint includes the coherence μ (verified across five independent
float64 kernels and cross-checked for internal consistency), the Welch slack
μ/Welch, the saturated-pair graph and its cliff D₁, the minimum determinant
across all *d*-subset Gram matrices (the *worst d-tuple*), its internal pair
morphology, the frame-operator tightness gap, a Bargmann discrete-phase test
on K₃ cliques in the saturated graph, the validity of the four standard
lower bounds (Welch–Rankin, orthoplex, Levenstein-2, and notes for
Bukh-Cox), and — when the cell (d, n) is present in the embedded
leaderboard snapshot — an automatic side-by-side comparison against the
current Game of Sloanes holder packing for that cell.

Everything runs locally in the browser. No data is uploaded anywhere. The
only outbound requests are optional fetches to
`raw.githubusercontent.com` to retrieve public packing files from the
Game of Sloanes repository for the holder-comparison feature.

---

## What this is not

This tool does not introduce new theorems, new bounds, or new optimality
proofs. The individual quantities it computes — μ, the Welch bound, the
frame operator, Gram subdeterminants viewed as RIP-style invariants,
Bargmann triple products — all have well-established precedents in the
cited literature. What the tool offers is the convenience of computing
them all together, in the browser, on any packing the user has on disk,
without installing anything.

The tool does not claim global optimality of any packing it analyses, and
it does not claim that float64 precision is sufficient for a Game of
Sloanes record submission at the 8th-decimal byte-exact rule.
Quad-precision verification (gcc `__float128` or Python `mpmath` at ≥ 30
digits) is required for record submissions. The Calculator tab reports
the maximum disagreement across its five kernels so the user can judge
when float64 is comfortable and when it is not.

---

## Who this is useful for

- Researchers in frame theory who want a structural fingerprint of a
  packing in a few seconds without firing up a Jupyter session.
- Students learning the geometry of Grassmannian configurations who want
  to look at what lies beyond the worst pair.
- Anyone preparing a Game of Sloanes submission who wants a
  pre-submission sanity check (with the float64 disclaimer above).
- Practitioners in compressed sensing, MIMO precoding, quantum
  tomography, or multiple description coding who use Grassmannian frames
  and want a quick structural reading of the configurations they work
  with.

The tool does not replace mathematical reasoning, quad-precision
verification, or peer review. It is one diagnostic instrument among
others. It emerged from empirical work on sub-catalog packings documented
in the companion repository
[sloane-coherence-records](https://github.com/tretoef-estrella/sloane-coherence-records),
where the same structural fingerprint conventions are applied to record
packings in cells (4, 64), (4, 48), and (3, 14) — readers interested in
worked examples will find the fingerprint of every record there.

---

## Three tabs

1. **Calculator.** Paste or upload a packing (.txt in the Game of Sloanes
   format) and get the full diagnostic. Auto-detection of (d, n) from the
   standard `dxn_xxx.txt` filename pattern. Export the result as a JSON
   report with SHA-256 hash of the input for downstream citation.

2. **Screener.** A static snapshot of a subset of the Game of Sloanes
   leaderboard (commit dated 2026-04-21 upstream). Each row links to the
   exact holder file in the upstream repository. A "Check upstream"
   button performs a live fetch of the upstream README and reports
   whether the leaderboard date has changed since the snapshot.

3. **About & Theory.** An explanation of every quantity the tool
   computes, why each one matters, the literature it draws on, and an
   honest list of what the tool does *not* claim. Designed to be read
   first by anyone who wants to understand what the numbers mean before
   trusting them.

---

## Technical notes

- Single HTML file, no build step, no server, no external JavaScript
  libraries. Google Fonts CSS is loaded for typography.
- Five independent kernels for μ: naive, Kahan-compensated, pairwise
  divide-and-conquer, Gram-triangular with reverse accumulation, and
  double-double (TwoProduct + TwoSum, roughly 31 extra bits of
  precision). The maximum disagreement across kernels is reported in the
  output so the user knows when float64 has reached its limit for a
  given packing.
- Hermitian eigenvalues computed via Jacobi rotation on the 2*d* × 2*d*
  real-block expansion of the *d* × *d* complex Gram. Adequate for the
  small *d* range typical of the Game of Sloanes leaderboard.
- Worst d-tuple enumeration limited to C(n, d) ≤ 700 000 to keep the
  browser responsive. Cells beyond that limit display the rest of the
  diagnostic and skip the worst-d-tuple section.
- Tested no-regression against several published holder packings from
  the leaderboard.

---

## License

This work is licensed under the
[Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

You are free to share and adapt this work for **non-commercial** purposes
(research, education, personal use), provided you give appropriate credit
and link to the license. Commercial use requires explicit permission —
contact `tretoef@gmail.com`.

Please cite this repository if you use the tool in published work:

```
Amichis Luengo, R. (2026). KADE Packing Diagnostic.
https://github.com/tretoef-estrella/KADE
```

---

## Author

Rafael Amichis Luengo (RAL), Madrid — `tretoef@gmail.com`

Independent investigator. The diagnostic framework was developed in
iterative collaboration with Claude (Anthropic) during empirical
cross-cell investigations of Grassmannian frame packings. The records
those investigations produced are documented in the companion repository
[sloane-coherence-records](https://github.com/tretoef-estrella/sloane-coherence-records).

---

## Acknowledgements

The Game of Sloanes leaderboard is curated by John Jasper, Emily J.
King, and Dustin G. Mixon at
[github.com/gnikylime/GameofSloanes](https://github.com/gnikylime/GameofSloanes).
The literature cited inside the *About & Theory* tab includes work by
Welch, Levenshtein, Rankin, Bukh–Cox, Strohmer–Heath, Casazza, Fickus–Mixon,
Bodmann–Haas, Benedetto–Fickus, Bargmann, Conway–Hardin–Sloane,
Cohn–Kumar, and Renes–Blume-Kohout–Scott–Caves. Their work is the
foundation on which any tool of this kind is built.
