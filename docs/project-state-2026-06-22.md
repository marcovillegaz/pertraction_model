# Project State Report — Pertraction Model
**Date:** 2026-06-22 | **Branch:** `feature_ms`

---

## 1. What actually runs

`main.m` is the only active entry point. It executes this chain:

```
main.m
  addpath(genpath('source'))          ← makes all source/ visible globally
  CompoundsLibrary(list, folder)
    └─ loadCompoundData(path)
         ├─ propertyFit()             ← fits density & viscosity on every load
         ├─ GCVOL60()                 ← group-contribution density estimate
         └─ sastriRao()               ← group-contribution viscosity estimate
  UNIFACLibrary(folder, file)
  UNIFACModel(compoundsLib, unifacLib)
  thermoModel.computeActivityCoefficient(T, x)
    └─ UNIFAC_test(compLib, thermoLib, T, x, 1000)
```

This chain is complete and tested. `Lngamma` output values at T=300K, x=[0.1, 0.6, 0.3] are very large (~66, ~71 for solvents) — physically suspicious and needs investigation.

There is a commented-out block in `main.m` referencing `computeFickDiffusivity()`. That function exists but is not yet callable (see §3).

---

## 2. The diffusion layer — exists, not wired

The full Maxwell-Stefan → Fick pipeline is coded in `source/calculations/diffusion-models/`. Intended call graph:

```
computeFickDiffusivity(compLib, unifacLib, T, x)     ← fick/
  mutualDiffusion(compLib, T, x)                      ← maxwellStefan/
    ├─ siddiqiLucas(T, compLib)                       ← infinite-dilution D, solvent pairs
    ├─ kooijmanTaylor(x, D0_inf)                      ← mutual D, non-polymer pairs
    ├─ vrentasVrentas(x, T, compLib)                  ← self-diffusion via FVT
    └─ kubaczka(x, D_self, D_mutual)                  ← mutual D, polymer pairs
  Bmatrix(D_ms, x)
  thermodynamicsFactors(compLib, unifacLib, T, x)
    └─ activityCoefficients('unifac-test', ...)
         └─ UNIFAC_test(...)
fickDiffusivity = B \ Gamma
```

`kooijmanTaylor`, `kubaczka`, and `Bmatrix` are self-contained (operate only on arrays). `siddiqiLucas` and `vrentasVrentas` are blocked by API mismatch (see §3.1).

---

## 3. Concrete bugs blocking the diffusion layer

### 3.1 — API mismatch in `siddiqiLucas.m` and `vrentasVrentas.m`

Both use the old procedural struct API:
- `fieldnames(compoundLibrary)` — assumes a plain struct, not a `CompoundsLibrary` object
- `extractPropertyAsArray(compoundLibrary, ...)` — calls the superseded global function in `source/processing/`

OOP equivalents already exist: `compoundLibrary.list()` and `compoundLibrary.extractPropertyAsArray(path)`. These two functions need updating before the diffusion layer can run.

### 3.2 — `MaxwellStefanModel.computeDiffusivity()` has two independent bugs

```matlab
A = 1e-9; n = 1.5;
D = A * temperature^n;  % Placeholder — physically meaningless
```

Additionally:
- References `obj.compoundLibrary` — but the property is named `compLib` in `DiffusionModel`
- Calls `compoundLibrary.getProperties(name)` — method does not exist; it is `compoundLibrary.get(name)`

This class would crash immediately on first call for both reasons.

### 3.3 — `thermodynamicsFactors.m` calls old `activityCoefficients()` wrapper

The wrapper correctly routes to `UNIFAC_test`, but passes `compoundLibrary` as a struct. This likely works only because `UNIFAC_test` internally uses `fieldnames()` on the struct stored inside `CompoundsLibrary`. Needs verification when wiring.

### 3.4 — Duplicate `mutualDiffusion.m`

Two files, same function name:
- `source/calculations/diffusion-models/mutualDiffusion.m` — older (2023), uses `fprintf` directly
- `source/calculations/diffusion-models/maxwellStefan/mutualDiffusion.m` — newer, uses `debugMsg`

`addpath(genpath('source'))` adds both directories. MATLAB resolves to whichever appears first — non-deterministic. The root-level copy must be deleted.

### 3.5 — `PerstractionModel.m` calls undefined functions

`applyFit()` calls `updateUnifac()` and `updateCompounds()` — neither exists. `computeMixtureProperties()` calls `arrayfun(@(c) c.getDensity(T), obj.compoundsLibrary)` — `CompoundsLibrary` is not an array of objects; `getDensity()` does not exist (it stores a function handle at `compound.density`).

### 3.6 — `Membrane.m` calls undefined functions

`initialize()` calls `computeXpoly()` and `computeMolarDensity()` — neither exists anywhere in the project.

---

## 4. Root-level clutter

The project root mixes six categories of files. Items that need relocation:

**Move to `obsolete/`:**

| File | What it is |
|---|---|
| `semiTransient_fit.m` | Old fitting script, 2023 |
| `semiTransient_fit_both.m` | Old fitting script, 2025 |
| `semiTransient_plot.m` | Old plot script, 2023 |
| `semiTransient_plot_both.m` | Old plot script, 2023 |
| `UNIFACtest.m` | Old standalone UNIFAC test, 2023 |
| `config.m` | Skeleton config function, not called from anywhere |
| `systemConstants.m` | PCB77 constants script, not wired to anything |

**Move to `data/` (appropriate subfolder):**

`.mat` files: `aq_fit`, `compoundData`, `compoundDataFVP`, `eqData`, `expCharge`, `expData_water_POMS_omim`, `expPerstract`, `ext_fit`, `flux_aq`, `flux_both`, `flux_ext`

`.xlsx` files: `expDensity`, `expViscosity`, `GCVOL`, `PCB77_perstract`, `SastriRao`, `variablesMatrix`

**Delete (MATLAB autosave — already gitignored):**
`Kubaczka.asv`, `main.asv`, `VrentasVrentas.asv`

**Images at root (move to `images/` or delete — gitignored):**
`fit_c.jpg`, `flux_fit.jpg`, `pcb77_FVP.jpg`, `pcb77FVP.jpg`

---

## 5. Documentation — scattered and incomplete

| File | Status |
|---|---|
| `source/README.md` | Broken links (`clasess/` typo; wrong relative paths) |
| `source/calculations/README.md` | Three bullet points, no content |
| `source/calculations/NOTES.md` | Useful — describes thermodynamic factor Γ matrix |
| `source/calculations/diffusion-models/README.md` | Partial, missing formulas |
| `source/calculations/diffusion-models/maxwellStefan/README.md` | References `inverseDrag.m` — file doesn't exist (it is `Bmatrix.m`) |
| `source/calculations/diffusion-models/fick/README.md` | Empty |
| `source/io/README.md` | Describes superseded procedural functions; not updated for OOP |
| `docs/NOTES.md` | Good FVT theory reference with equations and citations — keep |

---

## 6. `scripts/` and `main/` folders

**`scripts/processExperiments.m`** — processes PCB77 perstraction data, computes partition constants, saves `.mat` files. References `PCB77_MW` from `systemConstants.m` which is not on the path when called. Needs to be refactored as a function (currently a script).

**`main/FVP_fit.m`** — fits Free Volume Parameters via genetic algorithm. Loads `compoundData.mat` and accesses properties via `S.(component).FVP{row,col}` — entirely based on the old struct API. Incompatible with the current OOP layer.

---

## 7. Priority order for next work

| Priority | Action | Reason |
|---|---|---|
| 1 | Delete `source/calculations/diffusion-models/mutualDiffusion.m` | Duplicate causes non-deterministic function resolution |
| 1 | Delete `source/classes/PertractionModel.m` | Typo duplicate |
| 1 | Remove `disp(fh_vector)` at `CompoundsLibrary.m:147` | Debug statement in production code |
| 2 | Fix `siddiqiLucas.m` and `vrentasVrentas.m` | API mismatch blocks entire diffusion layer |
| 2 | Replace `MaxwellStefanModel.computeDiffusivity()` | Fix property/method errors + wire to real MS chain |
| 3 | Reorganize root-level files | Move docs → `docs/`, legacy → `obsolete/`, data → `data/` |
| 4 | Design `Membrane` and `PerstractionModel` | Undefined dependencies need architecture decision |
| 5 | Update READMEs across `source/` | Broken links, wrong function names, empty files |
