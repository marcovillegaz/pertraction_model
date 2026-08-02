# Perstraction model using Maxwell-stefan approuch. 
The perstracion model of my Msc. Thesis.

This code is a revisited version of the code that i made in Matlab for my thesis in 2023. The project structure was improved and new an enehece experimental data was added.

## Initualization procedure
The properties of each component are stored in individual Excel files, with each property organized into separate sheets within the file.

If experimental data for viscosity or density as a function of temperature is available, it is included as a table in the corresponding sheet. When such data is not available, the code estimates viscosity using the Sastri-Rao method and density using the GCVOL60 method.

The function loadCompoundData() imports and organizes the data for each compound into a structured format. Beyond that, the function initCompoundLibrary() builds a master struct that aggregates all the information for the system being modeled.

The function loadUnifacData() imports and organize unifac groups that describe the activity model for the specified system and store tha data as struct for activity coefficients computation. 


## Project structure

- [Source Overview](source/README.md) — MATLAB code: io, calculations, classes, processing, utils
- [Data Overview](data/README.md) — input / intermediate / output data lifecycle
- [Scripts Overview](scripts/README.md) — legacy System B (PCB77) workflow scripts
- [Theory](docs/theory/README.md) — mathematical reference (diffusion models, thermodynamics, property models), organized by physical model rather than by source file
- [Project log](docs/log/README.md) — chronological record of what happened in this project and why
- [Next steps](docs/next_steps/next_steps.md) — current, forward-looking ask/instructions
- `tests/` — MATLAB `matlab.unittest` suite, run via `cd tests; runAllTests`

```
.
├── main.m                        # Active entry point (System A: benzene + methyl acetate + polystyrene)
│
├── source/                       # Core logic
│   ├── io/                       # Excel/compound-data import
│   ├── processing/                # Superseded by CompoundsLibrary methods; kept for obsolete/ scripts
│   ├── calculations/              # Domain models, pure functions
│   │   ├── thermodynamic-models/  # UNIFAC legacy/FV/vdW-FV + shared terms (unifac/)
│   │   ├── diffusion-models/      # Maxwell-Stefan → Fick chain (fick/, maxwellStefan/)
│   │   └── property-models/       # Density/viscosity/molar-volume correlations
│   ├── classes/                   # OOP layer: Libraries, ThermoModel, DiffusionModel
│   └── utils/                     # projectPath, debugMsg, molarfraction
│
├── data/
│   ├── input/                    # Hand-curated/external; code never writes here
│   │   ├── compounds/             # One .xlsx per compound (System A + System B)
│   │   ├── reference/             # Literature correlation tables
│   │   ├── unifac/                # UNIFAC group + interaction-parameter tables
│   │   └── experimental/          # Raw measured data
│   ├── intermediate/             # Generated, consumed downstream; regenerable
│   │   ├── experiments/
│   │   ├── compound-structs/
│   │   └── property-fits/         # Gitignored
│   ├── output/                   # Terminal artifacts (figures, results); gitignored
│   └── test/                      # Literature reference figures (UNIFAC variant papers)
│
├── docs/
│   ├── theory/                    # Math reference (see Theory link above)
│   ├── log/                       # Project log (see Project log link above)
│   ├── next_steps/                # next_steps.md (current ask) + migration-plan.md
│   └── project-state-2026-06-22.md  # Dated snapshot, kept for history
│
├── tests/                        # matlab.unittest suite + fixtures/
├── scripts/                      # System B (PCB77) analysis workflow
├── main/                         # FVP_fit.m -- System B FVP fitting (GA)
├── config/                       # pertractSystem.m (System B constants)
├── obsolete/                     # G1/G2 legacy code, archived not deleted
├── diagrams/                     # Architecture diagrams (drawio)
└── .claude/skills/                # Project-scoped Claude Code skills (e.g. project-log)
```