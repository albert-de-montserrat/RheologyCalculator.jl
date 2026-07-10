# Plan: JOSS paper for RheologyCalculator.jl

## Context

`RheologyCalculator.jl` builds and solves local (point-wise) rheological
constitutive models assembled from viscous, elastic, and plastic building
blocks, composed in series / parallel / nested networks and solved with a
Newton–Raphson scheme. It targets geodynamics / Stokes solvers that must solve
the local rheological system at every quadrature point each time-step. It is
MIT-licensed, has Documenter-based docs, a test suite (unit tests + Aqua), CI,
and a large set of validated examples. This makes it a strong candidate for a
**JOSS** (Journal of Open Source Software) submission.

JOSS reviews the *software*; the paper is a short (~250–1000 word) advertisement
plus a "Statement of Need". The bulk of the work is ensuring the repository meets
JOSS submission requirements, then writing `paper.md` + `paper.bib`.

**Decisions locked in with the user:**
- **Authors:** Albert de Montserrat, Boris Kaus, Thibault Duretz (the three main
  git contributors; ORCIDs/affiliations to be confirmed during execution).
- **Scope:** produce a *full draft* of `paper.md` and `paper.bib`, not just a scaffold.

## Deliverables

1. This plan, stored in **two** locations (per the request):
   - `/Users/albert/.claude/plans/make-a-plan-to-sparkling-kazoo.md` (auto-loaded here).
   - `JOSS_PAPER_PLAN.md` in the repo root (for other agents). *Written on approval —
     plan mode forbids editing repo files now.*
2. `paper/paper.md` — JOSS paper with YAML front-matter + prose.
3. `paper/paper.bib` — BibTeX references.
4. (Optional) `paper/` figure(s) reusing existing `docs/assets/*.png`.
5. A short pre-submission checklist confirming repo meets JOSS criteria.

## JOSS submission requirements (verify during execution)

- [x] Open source, OSI license — MIT `LICENSE` present. ✅
- [x] Public VCS repo with substantial scholarly effort. ✅
- [x] Documentation — Documenter site under `docs/`. ✅ (confirm it covers
      install, API/usage, and a statement of functionality)
- [x] Automated tests — `test/runtests.jl`, Aqua, per-model tests. ✅
- [ ] **Community guidelines** — JOSS expects `CONTRIBUTING`/contribution +
      issue-reporting + support guidance. Repo has none → add a short
      `CONTRIBUTING.md` (or a README section) during execution.
- [ ] Version tagged + archived (Zenodo DOI) — needed at submission time, not draft time.
- [ ] Paper length: target **around 1000 words** (JOSS accepts 250–1000; aim for ~950–1000 to use the full budget).

## `paper.md` structure

YAML front-matter:
```yaml
---
title: 'RheologyCalculator.jl: Composable local rheological models for geodynamics'
tags:
  - Julia
  - geodynamics
  - rheology
  - constitutive models
  - viscoelastoplasticity
authors:
  - name: Albert de Montserrat
    orcid: 0000-0000-0000-0000   # confirm
    affiliation: 1
  - name: Boris J. P. Kaus
    affiliation: 2
  - name: Thibault Duretz
    affiliation: 3
affiliations:
  - name: ETH Zürich, Switzerland   # confirm
    index: 1
  - name: Johannes Gutenberg University Mainz, Germany
    index: 2
  - name: Goethe University Frankfurt, Germany   # confirm
    index: 3
date: <submission date>
bibliography: paper.bib
---
```

Prose sections:
1. **Summary** — what the package does, in accessible terms: composable
   viscous/elastic/plastic elements → series/parallel/nested networks →
   nonlinear residual solved by Newton iteration; allocation-free via Julia
   `@generated` / type dispatch; returns stress, strain-rate, pressure at a point.
2. **Statement of Need** — geodynamic Stokes solvers evaluate rheology at every
   quadrature point every time-step; hand-coding each viscoelastoplastic
   combination (Maxwell, Kelvin–Voigt, Burgers, VEP, VEVP, Drucker–Prager, cap
   plasticity, Cam-Clay, rate-and-state...) is error-prone and rigid. This
   package generates the local nonlinear system automatically from composed
   elements and solves it, decoupling the material catalogue from the solver.
   Position relative to the Julia geodynamics ecosystem (GeoParams.jl,
   JustRelax.jl, LaMEM/JustPIC) and to hand-written constitutive updates.
3. **Design / functionality** (brief) — state-function interface, `SeriesModel`/
   `ParallelModel`, `generate_equations`, `solve`, ForwardDiff Jacobian, elastic
   strain-rate correction. Keep concise — JOSS is not a manual.
4. **Example / validation** — one code snippet (Maxwell from README) + reference
   to examples validated against analytical solutions
   (`examples/Maxwell_KV_Maxwell.jl`, `docs/src/strain_rate_correction.md`).
   Optionally embed one figure from `docs/assets/`.
5. **Acknowledgements** — funding/ecosystem (confirm with user).
6. **References** — via `paper.bib`.

## `paper.bib` — references to gather

- Julia language (Bezanson et al. 2017).
- ForwardDiff.jl (Revels et al. 2016), StaticArrays.jl.
- Ecosystem: GeoParams.jl, JustRelax.jl, LaMEM (Kaus et al.), JustPIC.jl.
- Rheology/plasticity sources already cited in README: Popov et al. 2025
  (VEP+Cap), Abbo & Sloan 1995 (Hyperbolic), Golchin et al. 2021, de Souza Neto
  (Cam-Clay), Herrendörfer et al. 2018 (rate-and-state).
- General geodynamics rheology background (e.g. Gerya textbook, Moresi/May
  Stokes-solver refs) as needed for the Statement of Need.

## Key repo files to reuse (no core code changes needed)

- `README.md`, `docs/src/index.md` — source prose for Summary.
- `CLAUDE.md` — authoritative architecture description for the Design section.
- `examples/Maxwell_KV_Maxwell.jl`, `examples/Maxwell_VEP.jl`, etc. — validation.
- `docs/assets/*.png` — ready-made figures.
- `Project.toml` — deps for the software-paper acknowledgements/citations.

## Execution steps (after approval)

1. Write `JOSS_PAPER_PLAN.md` to repo root (copy of this plan for other agents).
2. Create `paper/` directory; write full `paper.md` (front-matter + all sections,
   within word limit) and `paper.bib`.
3. Optionally copy/reference one figure into `paper/`.
4. Add a minimal `CONTRIBUTING.md` (or README "Contributing" section) to satisfy
   JOSS community-guidelines criterion.
5. Leave ORCIDs, affiliations, funding, and submission date as clearly-marked
   TODOs for the user to confirm.

## Verification

- Word count of `paper.md` body is **around 1000** (target ~950–1000, JOSS hard cap 1000; `wc -w` excluding front-matter/bib).
- YAML front-matter parses (valid keys per JOSS schema; every `affiliation`
  index has a matching `affiliations` entry).
- Every in-text `[@key]` citation resolves to an entry in `paper.bib` (grep keys).
- Optionally build the paper PDF locally with the JOSS/Open Journals
  `inara`/`openjournals` Docker action to confirm it compiles, or rely on the
  JOSS "Compile PDF" bot after opening the submission.
- Repo still passes existing checks (no source changes, so `test/runtests.jl`
  unaffected).
```
