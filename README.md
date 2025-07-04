# Perstraction model using Maxwell-stefan approuch. 
The perstracion model of my Msc. Thesis.

This code is a revisited version of the code that i made in Matlab for my thesis in 2023. The project structure was improved and new an enehece experimental data was added.

## Initualization procedure
The properties of each component are stored in individual Excel files, with each property organized into separate sheets within the file.

If experimental data for viscosity or density as a function of temperature is available, it is included as a table in the corresponding sheet. When such data is not available, the code estimates viscosity using the Sastri-Rao method and density using the GCVOL60 method.

The function loadCompoundData() imports and organizes the data for each compound into a structured format. Beyond that, the function initCompoundLibrary() builds a master struct that aggregates all the information for the system being modeled.

The function loadUnifacData() imports and organize unifac groups that describe the activity model for the specified system and store tha data as struct for activity coefficients computation. 


# Project structure
- [Source Overview](source/README.md)
- [Data Overview](data/README.md)
- [Scripts Overview](scripts/README.md)











# Recommended project structure folder
/MyProject
├── main.m                        % Entry point script
│
├── /src                          % Core logic
│   ├── /io                       % All import/export functions
│   │   ├── readCSV.m
│   │   └── loadExperimentalData.m
│   │
│   ├── /processing               % Data preprocessing and transformations
│   │   ├── filterOutliers.m
│   │   └── normalizeSignal.m
│   │
│   ├── /calculations             % Scientific or domain-specific models
│   │   ├── computeDiffusion.m
│   │   └── estimateSlope.m
│   │
│   └── /utils                    % General helpers (used everywhere)
│       └── printTable.m
│
├── /data                         % Raw input and reference data
│   └── experiment_01.csv
│
├── /scripts                      % Scripts for running full workflows
│   └── analyzeExperiment.m
│
├── /results                      % Output: figures, tables, etc.
│
└── README.md