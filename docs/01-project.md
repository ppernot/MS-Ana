# Project structure

The project is organized into the following folders structure:

    ├── analysis :   contains the R scripts 
    │               (analysis.R, checkRep.R, quantify.R)
    │               and their auxillary functions
    │
    ├── data :       default folder for the input tables and 
    │               data to be analyzed
    │
    ├── results :    where the outputs of the scripts are stored
    │   │
        ├── figs :    figures
        │
        └── tables :  tables

The scripts should be run from the `analysis` folder, and the paths to
the required folders are defined in the scripts as:

    # Define Data and Results repositories
    dataRepo = '../data/'
    figRepo  = '../results/figs/'
    tabRepo  = '../results/tables/'

The MS and DMS files are expected by default to be in `data`. For
complex projects, the MS files can be placed in sub-folders of `data`,
and their paths are given in the `taskTable` input file.
