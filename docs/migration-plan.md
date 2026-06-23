# Migration Plan — Pertraction Model OOP Refactor
**Date:** 2026-06-22 | **Branch:** `feature_ms`

---

## 1. Where the project came from

The codebase contains three generations of code that coexist in the same repository:

| Generation | Period | What it looks like |
|---|---|---|
| **G1 — Scripts** | 2023 | `semiTransient_*.m` at root, monolithic `compoundData.mat` struct, everything hardcoded |
| **G2 — Procedural** | 2024–2025 | Functions in `source/calculations/`, compound data in Excel files per compound, fitting in `main/FVP_fit.m` |
| **G3 — OOP** | 2025–present | Classes in `source/classes/`, library objects, model objects — thermodynamics complete, diffusion stubbed |

The migration was going in the right direction. The G2 calculation functions are clean and mathematically correct. The G3 OOP layer is a good design. The problem is that the migration stalled halfway: G3 classes were created but not wired to G2 functions, and G1/G2 artifacts were never removed from the root.

---

## 2. What belongs to which system

Two separate physical systems are mixed in the repository. Understanding this is essential before reorganizing files.

**System A — Benzene + methyl acetate + polystyrene**
This is the current active system. `main.m` runs this. Test compound Excel files in `data/test-compounds/`. The OOP migration targets this system.

**System B — PCB77 + water/extractant + POMS membrane**
This was the original research campaign. Most root-level files belong here:
- `PCB77_perstract.xlsx`, `systemConstants.m` (PCB77 MW)
- `expCharge.mat`, `expPerstract.mat`, `expData_water_POMS_omim.mat`
- `flux_aq.mat`, `flux_ext.mat`, `flux_both.mat`, `aq_fit.mat`, `ext_fit.mat`
- `semiTransient_fit.m`, `semiTransient_fit_both.m`, `semiTransient_plot.m`, `semiTransient_plot_both.m`
- `scripts/processExperiments.m`, `scripts/fit_plot.m`, `scripts/flux_plot.m`
- `FVP_results.txt`, `FVPResultsWENO.txt`, `waterPOMS_fit_c.txt`
- Compound fit images in `images/` (acetonitrile, omimtf2n, water, POMS)

System B code is not being migrated to OOP — it ran, produced results, and is now legacy. It should be archived, not deleted, in case experimental data or parameters need to be recovered.

---

## 3. Current state of each layer

### 3.1 OOP layer (`source/classes/`)

| Class | State | Notes |
|---|---|---|
| `CompoundsLibrary` | **Working** | Loads Excel files, stores density/viscosity as `@(T)` handles. Minor bug: `disp(fh_vector)` at line 147. |
| `ThermoLibrary` / `UNIFACLibrary` | **Working** | Loads UNIFAC groups and interaction matrix. |
| `ThermoModel` / `UNIFACModel` | **Working** | Wraps `UNIFAC_test`. Tested in `main.m`. |
| `DiffusionModel` | **Stub** | Abstract base class. Correct structure. |
| `MaxwellStefanModel` | **Broken stub** | Placeholder formula. Two property/method name bugs. Cannot call. |
| `Membrane` | **Empty skeleton** | Calls two undefined functions. No design yet. |
| `PerstractionModel` | **Broken skeleton** | Has more code but calls three undefined functions and uses wrong API for `CompoundsLibrary`. |
| `PertractionModel.m` | **Delete** | Typo duplicate of `PerstractionModel.m`. |

### 3.2 Calculation layer (`source/calculations/`)

| Function | State | Notes |
|---|---|---|
| `UNIFAC_test.m` | **Active** | Called by `UNIFACModel`. Working. |
| `UNIFAC_vdW_FV.m` | **Exists, not wired** | Older UNIFAC variant. |
| `activityCoefficients.m` | **Exists, not wired** | Router function; called by `thermodynamicsFactors.m` only. |
| `computeFickDiffusivity.m` | **Exists, not wired** | Top-level entry for Fick matrix. Calls the full MS chain. Commented out in `main.m`. |
| `mutualDiffusion.m` (maxwellStefan/) | **Exists, not wired** | Correct version. Calls siddiqiLucas, kooijmanTaylor, vrentasVrentas, kubaczka. |
| `mutualDiffusion.m` (diffusion-models/) | **Delete** | Duplicate of above. On MATLAB path simultaneously — non-deterministic. |
| `Bmatrix.m` | **Exists, not wired** | Correct. Self-contained. |
| `thermodynamicsFactors.m` | **Exists, not wired** | Calls old `activityCoefficients()`. Needs API check. |
| `siddiqiLucas.m` | **Blocked** | Uses old struct API. Must be updated before diffusion chain can run. |
| `kooijmanTaylor.m` | **Ready** | Self-contained array math. No changes needed. |
| `vrentasVrentas.m` | **Blocked** | Uses old struct API. Must be updated before diffusion chain can run. |
| `kubaczka.m` | **Ready** | Self-contained array math. No changes needed. |
| `propertyFit.m`, `GCVOL60.m`, `sastriRao.m` | **Active** | Called by `loadCompoundData`. Working. |
| `TynCalus.m` | **Exists, unknown** | Not traced in active chain. |

### 3.3 IO layer (`source/io/`)

| File | State |
|---|---|
| `loadCompoundData.m` | **Active** — called by `CompoundsLibrary` |
| `initCompoundLibrary.m` | **Superseded** — replaced by `CompoundsLibrary` constructor |
| `loadUnifacData.m` | **Superseded** — replaced by `UNIFACLibrary` |

### 3.4 Processing layer (`source/processing/`)

| File | State |
|---|---|
| `extractPropertyAsArray.m` | **Superseded** — replaced by `CompoundsLibrary.extractPropertyAsArray()` |
| `extractPropertyAsFunction.m` | **Superseded** — replaced by `CompoundsLibrary.extractPropertyAsFunction()` |

---

## 4. The missing wiring — what needs to happen

The G2 diffusion chain is complete and correct, but `siddiqiLucas` and `vrentasVrentas` still use the old G1/G2 struct-based API. The fix is localized: update those two functions to call the `CompoundsLibrary` OOP methods, then wire `computeFickDiffusivity` into `MaxwellStefanModel`.

The full path that needs to work:

```
MaxwellStefanModel.computeDiffusivity(T, x)
  └─ computeFickDiffusivity(obj.compLib, unifacLib, T, x)
       ├─ mutualDiffusion(compLib, T, x)
       │    ├─ siddiqiLucas(T, compLib)         ← needs OOP API fix
       │    ├─ kooijmanTaylor(x, D0)             ← ready
       │    ├─ vrentasVrentas(x, T, compLib)     ← needs OOP API fix
       │    └─ kubaczka(x, D_self, D_mutual)     ← ready
       ├─ Bmatrix(D_ms, x)                       ← ready
       └─ thermodynamicsFactors(compLib, unifacLib, T, x)  ← verify API
```

---

## 5. Recommendations

### Priority 1 — Immediate cleanup (no design decisions required)

These are purely mechanical. Do before anything else to eliminate confusion and the non-deterministic path bug.

1. **Delete** `source/calculations/diffusion-models/mutualDiffusion.m` — duplicate causes non-deterministic function resolution
2. **Delete** `source/classes/PertractionModel.m` — typo duplicate
3. **Remove** `disp(fh_vector)` at `CompoundsLibrary.m:147`
4. **Delete** `main.asv`, `Kubaczka.asv`, `VrentasVrentas.asv` — MATLAB autosave files (gitignored, physically present)

### Priority 2 — Wire the diffusion layer

This is the core of the current sprint. Two files to fix, then one class to complete.

**Step 2a** — Fix `siddiqiLucas.m`

Replace:
```matlab
compound_names = fieldnames(compoundLibrary);
viscosityFuncs = extractPropertyAsArray(compoundLibrary, {"Viscosity"})
boilingVolumes = extractPropertyAsArray(compoundLibrary, {"boilingPoint",1})
```
With:
```matlab
compound_names = compoundLibrary.list();
viscosityFuncs = compoundLibrary.extractPropertyAsArray({"Viscosity"});
boilingVolumes = compoundLibrary.extractPropertyAsArray({"boilingPoint",1});
```

**Step 2b** — Fix `vrentasVrentas.m` (same pattern)

Replace all `fieldnames(compoundLibrary)` and `extractPropertyAsArray(compoundLibrary, ...)` calls with `compoundLibrary.list()` and `compoundLibrary.extractPropertyAsArray(...)`.

**Step 2c** — Verify `thermodynamicsFactors.m`

It calls `activityCoefficients('unifac-test', compoundLibrary, unifacLibrary, x, T, n)` which routes to `UNIFAC_test`. Confirm the signature is compatible with the current `UNIFACLibrary` object, then decide whether to keep the wrapper or call `UNIFACModel` directly.

**Step 2d** — Rewrite `MaxwellStefanModel.computeDiffusivity()`

Fix property name (`compLib`, not `compoundLibrary`), remove placeholder, wire to `computeFickDiffusivity`. The class also needs access to `unifacLib` — decide whether to pass it at construction or at call time.

**Step 2e** — Test the full chain in `main.m`

Uncomment the `computeFickDiffusivity` block. Run at T=300K, x=[0.1, 0.6, 0.3]. Validate output dimensions and order of magnitude.

### Priority 3 — Reorganize root-level files

Once the diffusion chain runs, the old files become clearly identifiable as legacy. Group by destination:

**Archive in `data/experimental-data/pcb77/`** (System B experiment data):
`PCB77_perstract.xlsx`, `expCharge.mat`, `expPerstract.mat`, `expData_water_POMS_omim.mat`, `flux_aq.mat`, `flux_ext.mat`, `flux_both.mat`, `aq_fit.mat`, `ext_fit.mat`

**Archive in `data/fit-results/`** (fitting outputs):
`compoundData.mat`, `compoundDataFVP.mat`, `eqData.mat`, `FVP_results.txt`, `FVPResultsWENO.txt`, `waterPOMS_fit_c.txt`

**Archive in `data/reference/`** (property reference tables used by old scripts):
`GCVOL.xlsx`, `SastriRao.xlsx`, `expDensity.xlsx`, `expViscosity.xlsx`, `variablesMatrix.xlsx`

**Move to `obsolete/`** (superseded scripts):
`semiTransient_fit.m`, `semiTransient_fit_both.m`, `semiTransient_plot.m`, `semiTransient_plot_both.m`, `UNIFACtest.m`, `config.m`, `systemConstants.m`

**Delete** (gitignored generated outputs physically present):
`fit_c.jpg`, `flux_fit.jpg`, `pcb77_FVP.jpg`, `pcb77FVP.jpg`

### Priority 4 — Design `Membrane` and `PerstractionModel`

These require architectural decisions, not just wiring. The key questions:

- What does `Membrane` own? Likely: geometry (area, thickness), composition (molar fractions, total concentration), and a reference to `MaxwellStefanModel`. The `computeDiffusivity()` call should come through `MaxwellStefanModel`, not be reimplemented here.
- What does `PerstractionModel` own? Likely: the two bulk phases, the membrane, equilibrium constants (partition coefficients), and the ODE system for the semi-transient mass balance. The old `semiTransient_fit.m` is the reference implementation.
- How are partition coefficients handled? They appear in `PerstractionModel` as `equilibriumConstants` — currently assumed to be fitted externally and passed in. This needs to be explicit.

Create an ADR before implementing these classes.

### Priority 5 — Migrate `FVP_fit.m` and `processExperiments.m` to OOP

Both workflows (`main/FVP_fit.m` and `scripts/processExperiments.m`) use the old G1 struct API and cannot be called from the new layer. `FVP_fit.m` in particular is critical — it produces the Free Volume Parameters that `vrentasVrentas.m` depends on.

`FVP_fit.m` should become a method or script that takes a `CompoundsLibrary` object, fits FVP per compound, and writes the parameters back into the compound Excel files (or a separate FVP file). This removes the dependency on the old `.mat` struct format.

---

## 6. What to ignore

The following files exist and can be left alone indefinitely. They are not blocking anything and are not worth cleaning up until the active model is working end-to-end:

- `source/io/initCompoundLibrary.m` and `loadUnifacData.m` — superseded but harmless
- `source/processing/extractPropertyAsArray.m` and `extractPropertyAsFunction.m` — superseded; kept so old scripts in `obsolete/` remain readable
- `scripts/finalFit_c.m`, `fit_plot.m`, `flux_plot.m` — System B analysis scripts, archive when System B files are moved
- `config/` folder — investigate before deciding; may contain something relevant or be entirely empty

---

## 7. Target state

When complete, `main.m` should run end-to-end:

```matlab
compLib  = CompoundsLibrary(COMPOUNDS_LIST, COMPOUNDS_FOLDER);
unifacLib = UNIFACLibrary("data/unifac-data", "unifac-test.xlsx");

thermoModel  = UNIFACModel(compLib, unifacLib);
diffusModel  = MaxwellStefanModel(compLib, unifacLib);
membrane     = Membrane(area, thickness, diffusModel, thermoModel);
model        = PerstractionModel(membrane, systemConfig);

Lngamma      = thermoModel.computeActivityCoefficient(T, x);
D_fick       = diffusModel.computeDiffusivity(T, x);
[flux, conc] = model.simulate(experimentalData);
```

Everything above that line is already working. `D_fick` is the next milestone.
