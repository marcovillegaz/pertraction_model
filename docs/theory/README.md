# Theory

Mathematical foundations for the pertraction model, organized by physical
model rather than by source file. Each doc here is the authority for its
equations — implementation status, bugs, and TODOs belong in
[`docs/migration-plan.md`](../migration-plan.md) and
[`docs/project-state-2026-06-22.md`](../project-state-2026-06-22.md), not here.

- [`nomenclature.md`](nomenclature.md) — shared symbol table. Every doc below uses these symbols instead of redefining them.
- [`references.md`](references.md) — consolidated bibliography, cited by author-year from the docs below.

## Diffusion models

The multicomponent mass-transfer chain follows Maxwell-Stefan theory
(Taylor & Krishna, 1993), assembled into an effective Fick diffusivity.
Call graph (see [`fick.md`](diffusion-models/fick.md) for the assembly step):

```
computeFickDiffusivity(compLib, unifacLib, T, x)
  ├─ mutualDiffusion(compLib, T, x)                     → dilute-correlations.md + maxwell-stefan.md
  │    ├─ siddiqiLucas      (D at infinite dilution, non-polymer pairs)
  │    ├─ kooijmanTaylor    (mutual D, non-polymer pairs)
  │    ├─ vrentasVrentas    (self-diffusion via Free Volume Theory)   → free-volume-theory.md
  │    └─ kubaczka          (mutual D, polymer pairs)                 → maxwell-stefan.md
  ├─ Bmatrix(D_ms, x)                                    → maxwell-stefan.md
  └─ thermodynamicsFactors(compLib, unifacLib, T, x)     → maxwell-stefan.md
       └─ activityCoefficients → UNIFAC                  → thermodynamic-models/unifac.md
fickDiffusivity = B \ Gamma
```

- [`diffusion-models/dilute-correlations.md`](diffusion-models/dilute-correlations.md) — infinite-dilution and mutual diffusion for non-polymer pairs (Siddiqi & Lucas; Kooijman & Taylor)
- [`diffusion-models/free-volume-theory.md`](diffusion-models/free-volume-theory.md) — self-diffusion in the polymer-solvent system (Vrentas & Vrentas FVT)
- [`diffusion-models/maxwell-stefan.md`](diffusion-models/maxwell-stefan.md) — polymer mutual diffusion (Kubaczka), the inverse-drag $[B]$ matrix, and the thermodynamic factor $[\Gamma]$ matrix
- [`diffusion-models/fick.md`](diffusion-models/fick.md) — assembly into the effective Fick diffusivity matrix

## Thermodynamic models

- [`thermodynamic-models/unifac.md`](thermodynamic-models/unifac.md) — UNIFAC-FV activity coefficient model

## Property models

- [`property-models/density-viscosity.md`](property-models/density-viscosity.md) — group-contribution density (GCVOL), viscosity (Sastri & Rao), and molar volume at boiling point (Tyn & Calus)

## How this maps to `source/`

Each doc above documents a *model*, which may span several `.m` files. The
top of each doc lists exactly which files in `source/calculations/` it
governs. When you add a new correlation to an existing model (e.g. another
dilute-diffusion correlation), add a section to the existing doc — don't
create a new file per function. Create a new doc only when a genuinely new
physical model enters the pipeline.

## Two physical systems

The math in this folder is model-level and applies to both systems
currently in the repo (see [`docs/migration-plan.md`](../migration-plan.md) §2):

- **System A** (benzene + methyl acetate + polystyrene) — the active OOP
  migration target. Uses the full chain above.
- **System B** (PCB77 + water/extractant + POMS membrane) — legacy,
  archived in `obsolete/` and `data/experimental-data/pcb77/`. Not
  reorganized here; its math, where it diverges, stays with the legacy
  scripts.
