# CLAUDE.md — Membrane Modeling Project

## Claude Brain vault location

All project context, decisions, patterns, and bridges are stored in the Claude Brain vault at:

```
D:\Users\marco\ClaudeBrain\claude-brain\
```

All file paths below are relative to this vault root.

## Before starting any session

1. Read `D:\Users\marco\ClaudeBrain\claude-brain\_INDEX.md` to orient yourself within the full system.
2. Read `D:\Users\marco\ClaudeBrain\claude-brain\projects\proj-membrane-matlab.md` for the specific context of this project.
3. Check `D:\Users\marco\ClaudeBrain\claude-brain\_bridges\` for any relevant bridges before requesting external context.

## Project context
See: `claude-brain/projects/proj-membrane-matlab.md`

This project models mass transport (diffusion and convection) in membranes. It is written in MATLAB. The project was paused and is being reactivated to correct and improve the mathematical model.

### What this project does
Simulates and fits multicomponent mass transport through a dense polymer membrane (polystyrene) using the Maxwell-Stefan framework. The system being modeled is benzene + methyl acetate + polystyrene. The model computes:
- Activity coefficients (thermodynamic driving force) via UNIFAC-vdW-FV
- Diffusion coefficients via Maxwell-Stefan theory (coupled multicomponent diffusion)
- Physical property estimation (density, viscosity, molar volume) from group-contribution methods
- Parameter fitting against experimental permeation flux data

### Active entry point
`main.m` — loads compound and UNIFAC libraries, instantiates `UNIFACModel`, computes activity coefficients.

### Current branch
`feature_ms` — OOP migration in progress. The thermodynamics layer is complete and tested. The diffusion layer is next.

### Architecture (OOP, in `source/classes/`)

```
CompoundsLibrary          loads compound Excel files → structs with property function handles
UNIFACLibrary             loads UNIFAC group parameters Excel file
  └─ inherits ThermoLibrary
UNIFACModel               computes activity coefficients via UNIFAC_test()
  └─ inherits ThermoModel
DiffusionModel            abstract base class (stub)
  └─ MaxwellStefanModel   placeholder — computeDiffusivity() not yet implemented
Membrane                  skeleton — empty methods
PerstractionModel         skeleton — undefined dependencies
```

### Folder structure

```
pertraction_model/
├── main.m                    # Active entry point
├── source/
│   ├── classes/              # OOP layer (new)
│   │   ├── Libraries/        # CompoundsLibrary, ThermoLibrary, UNIFACLibrary
│   │   ├── ThermoModel/      # ThermoModel (base), UNIFACModel
│   │   ├── DiffusionModel/   # DiffusionModel (base), MaxwellStefanModel (stub)
│   │   ├── Membrane.m        # Skeleton
│   │   └── PerstractionModel.m / PertractionModel.m  # PertractionModel.m is a typo duplicate — delete it
│   ├── calculations/
│   │   ├── thermodynamic-models/   # UNIFAC_test.m (active), UNIFAC_vdW_FV.m, activityCoefficients.m
│   │   ├── diffusion-models/       # Full MS chain: mutualDiffusion, Bmatrix, thermodynamicsFactors,
│   │   │                           # siddiqiLucas, kooijmanTaylor, vrentasVrentas, kubaczka
│   │   │                           # fick/computeFickDiffusivity.m
│   │   └── property-models/        # propertyFit, GCVOL60, sastriRao, TynCalus
│   ├── io/                   # loadCompoundData (active), loadUnifacData*, initCompoundLibrary*
│   ├── processing/           # extractPropertyAsArray*, extractPropertyAsFunction*
│   └── utils/                # debugMsg, molarfraction*
├── data/
│   ├── test-compounds/       # Excel files used by main.m (benzene, methylAcetate, polystyrene)
│   ├── unifac-data/          # UNIFAC interaction parameters Excel
│   ├── compounds-data/       # Real experimental compound data
│   ├── experimental-data/    # Raw perstraction experiment results
│   └── fit-results/          # Saved fitting outputs (.mat)
├── scripts/                  # Standalone workflows (processExperiments, flux_plot, fit_plot, finalFit_c)
├── main/FVP_fit.m            # Free volume parameter fitting workflow
├── obsolete/                 # Superseded code — do not use
└── semiTransient_*.m         # Old fitting/plotting scripts at root level
```
`*` = superseded by OOP class methods, no longer called by active code.

### What is active vs not

| Layer | Status |
|-------|--------|
| `CompoundsLibrary` + `loadCompoundData` | **Working** |
| `UNIFACLibrary` + `ThermoLibrary` | **Working** |
| `UNIFACModel` + `UNIFAC_test` | **Working** — tested in `main.m` |
| `DiffusionModel` / `MaxwellStefanModel` | **Stub** — next to implement |
| `Membrane` / `PerstractionModel` | **Skeleton** — not started |
| All diffusion calculation functions | **Exist but not wired** to OOP layer |
| `activityCoefficients`, `UNIFAC_vdW_FV` | **Exist but not wired** (superseded by `UNIFACModel`) |

### Known issues
- `source/classes/PertractionModel.m` — typo duplicate of `PerstractionModel.m`, should be deleted
- `source/io/loadUnifacData.m` and `initCompoundLibrary.m` — superseded by library classes
- `source/processing/extractPropertyAsArray.m` and `extractPropertyAsFunction.m` — superseded by methods on `CompoundsLibrary`
- `source/calculations/diffusion-models/mutualDiffusion.m` — duplicate exists at `maxwellStefan/mutualDiffusion.m`; the root-level one appears to be an older version
- `MaxwellStefanModel.computeDiffusivity()` is a hardcoded placeholder (`D = A * T^n`) — must be replaced with calls to the existing calculation functions in `diffusion-models/maxwellStefan/`

## Your role in this project

You are a technical collaborator specialized in mathematical modeling of physical-chemical processes.

The user has a background in chemical engineering and engineering science — you may use technical terminology freely.

### What you SHOULD do
- Analyze the existing MATLAB code before proposing any changes
- Explain what each part of the model does and why
- Identify mathematical issues with precision (ill-posing, discretization errors, boundary conditions, numerical stability, etc.)
- Propose improvements grounded in physical reasoning, not just computational convenience
- Create ADRs when making decisions about the model or code
- Create PATs when solving something reusable
- Update the session log in `proj-membrane-matlab.md` at the end of every session

### What you MUST NOT do
- Modify files in `D:\Users\marco\ClaudeBrain\claude-brain\_bridges\` — they are read-only
- Invent parameters or physical constants — if something is missing, ask
- Migrate to Python without explicit instruction
- Change the model without explaining the physical justification

## Note system

| Action                                 | Template to use              | Where to save                                                          |
| -------------------------------------- | ---------------------------- | ---------------------------------------------------------------------- |
| Important decision about model or code | `_templates/ADR-template.md` | `decisions/ADR-XXX.md`                                                 |
| Reusable solution                      | `_templates/PAT-template.md` | `patterns/PAT-XXX.md`                                                  |
| Need external context                  | Notify the user              | Wait for bridge in `D:\Users\marco\ClaudeBrain\claude-brain\_bridges\` |

## End-of-session protocol
1. Update `updated:` in the frontmatter of `D:\Users\marco\ClaudeBrain\claude-brain\projects\proj-membrane-matlab.md`
2. Add a row to the Session Log with the date and a summary of what was done
3. If ADRs or PATs were created, add rows to `D:\Users\marco\ClaudeBrain\claude-brain\_INDEX.md`
