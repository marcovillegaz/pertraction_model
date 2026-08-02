# UNIFAC multi-model refactor — 2026-08-01

## Context

`docs/next_steps.md` asked for three selectable UNIFAC variants (legacy,
UNIFAC-FV, UNIFAC-vdW-FV), a way to pick a variant on the model class, and
a test workflow. Investigation found the repo already had two near-duplicate
~200-line implementations of one variant (`UNIFAC_test.m` / `UNIFAC_vdW_FV.m`,
differing only in a swapped argument order and library API), true legacy
UNIFAC existing only in `obsolete/UNIFACtest.m`, and live code
(`thermodynamicsFactors.m`) depending on a broken function that had been
moved to `obsolete/`.

## What changed

- **Shared terms extracted** into `source/calculations/thermodynamic-models/unifac/`:
  `unifacGroupParams.m`, `unifacCombinatorial.m`, `unifacResidual.m`,
  `unifacFreeVolume.m`. The residual term (~80 lines, previously
  byte-identical in both old files) now exists once.
- **Three variant functions** compose those terms: `unifacLegacy.m` (C+R),
  `unifacFV.m` (C+R+FV, replaces `UNIFAC_test.m`), `unifacVdwFV.m`
  (C+R+FV, replaces `UNIFAC_vdW_FV.m` — see Outstanding below).
- **Class hierarchy**: `ThermoModel` (abstract base, `methods (Abstract)`
  pattern matching `DiffusionModel`) → `UNIFACModel` (abstract, holds
  polymerization degree) → `UNIFACLegacyModel` / `UNIFACFVModel` /
  `UNIFACvdWFVModel`. Each declares `RequiredCompoundProps` (Constant) —
  legacy needs nothing extra, FV and vdW-FV need `FVP` and `density`.
- **Factory**: `ThermoModel.create(variantName, compLib, unifacLib, opts)`,
  used in `main.m` instead of constructing a class directly.
- **Validation**: `ThermoModel.validateCompoundData()` runs inside
  `create()`, checking every compound against the variant's
  `RequiredCompoundProps` — including catching *present but empty* FVP
  rows, not just missing fields. Verified against two real compound files:
  `omimTf2N.xlsx` (no `FVP` group at all) and `Water.xlsx` (`FVP` present,
  rows 3/4/5 empty).
- **Fixed the live→obsolete dependency**: `thermodynamicsFactors.m` now
  takes a `ThermoModel` object directly (`thermodynamicsFactors(thermoModel, T, x)`)
  instead of routing through `activityCoefficients.m`. That wrapper had
  moved to `obsolete/` and had a live argument-order bug — it passed
  `(x, T)` into `UNIFAC_test`'s `(T, x)` slots, so anything reaching it
  computed garbage. Likely related to the "physically suspicious" Lngamma
  magnitudes noted in `docs/project-state-2026-06-22.md`.
- **`data/` reorganized** (separate from but preceding this work) into
  `input/` / `intermediate/` / `output/`, with a new `source/utils/projectPath.m`
  helper so paths resolve from the repo root instead of the CWD.
- **Tests**: `tests/` using `matlab.unittest` — `tUnifacGolden.m` (regression
  against a pre-refactor golden fixture), `tUnifacVariants.m` (shape/finite
  checks + factory dispatch for all three variants), `tCompoundValidation.m`
  (the two real-world validation-failure cases above),
  `tUnifacNayakAkhouri.m` (placeholder — loads `test_nayakAkhouri.xlsx`,
  no reference values yet). Run via `cd tests; runAllTests`.
- **Deleted** `UNIFAC_test.m` and `UNIFAC_vdW_FV.m` (fully superseded).
- **Docs**: `docs/theory/thermodynamic-models/unifac.md` rewritten with one
  section per variant, each citing its source paper; citations added to
  `docs/theory/references.md`.

## Decisions

- **Factory + subclasses, not a router switch.** The repo had already tried
  a string-dispatch router (`obsolete/activityCoefficients.m`) and
  abandoned it with a live bug in place. Subclasses let each variant
  declare its own compound-data requirements; a router can't express that
  without a growing central `switch`.
- **Shared terms extracted before writing new variants**, verified against
  a golden fixture (`tests/fixtures/unifac_fv_golden.mat`,
  `[65.9712; 71.2320; -14.4071]` at T=300K, x=[0.1 0.6 0.3]) at each step,
  so the refactor is provably behavior-preserving up to the point where new
  behavior was intentionally added.
- **Did not invent the UNIFAC-vdW-FV physics.** Both pre-existing
  implementations computed the same free-volume formula despite being
  named differently, and there was no reference distinguishing them.
  Ported the code as-is and marked it `NEEDS REFERENCE`/`NEEDS FORMULA`
  rather than guess.

## Outstanding / next

- **UNIFAC-vdW-FV free-volume term is still a placeholder.** It reuses
  `unifacFreeVolume.m` (UNIFAC-FV's Bondi hard-core-volume term) unchanged,
  so `UNIFACvdWFVModel` currently produces output numerically identical to
  `UNIFACFVModel`. *Update since this entry's source content in
  `docs/next_steps.md`*: the citation has since been supplied — Kannan,
  Duda & Danner (2005), "A free-volume term based on the van der Waals
  partition function for the UNIFAC model" — and is now in
  `docs/theory/references.md` and cited in `unifac.md`. The actual
  equations still need to be extracted into `unifac.md` and implemented as
  `unifacFreeVolumeVdw.m`.
- `tests/tUnifacNayakAkhouri.m` needs real expected values from the
  Nayak & Akhouri literature case once available.
- The "physically suspicious" Lngamma magnitudes (large positive values for
  benzene/methylAcetate) were not investigated — out of scope for this
  refactor, which preserved them exactly via the golden-fixture check.

## Related

- [`docs/theory/thermodynamic-models/unifac.md`](../theory/thermodynamic-models/unifac.md)
- [`docs/theory/diffusion-models/maxwell-stefan.md`](../theory/diffusion-models/maxwell-stefan.md) (consumer of `thermodynamicsFactors`)
- [`docs/migration-plan.md`](../migration-plan.md) (Priority 2c, marked resolved)
- [`docs/next_steps.md`](../next_steps.md) (original ask; source of this entry)
