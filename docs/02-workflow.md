# Workflow

## 0- Pretreatment

1.  FT-ICR MS files should be compressed using `pretreatFTICR.R` to
    enable analysis

## 1- Analysis

1.  put MS files in `data` folder or a sub-folder

2.  put DMS files in `data`

3.  create `taskTable` and `tgTable` input files in `data` folder

4.  edit and run `analysis.R`

    -   enter the file paths to the `taskTable` and `tgTable` files  
    -   choose the type of peak fit (variable `fit_dim`)

## 2- Quantification

1.  create `quantTable` in `data` folder

2.  edit and run `checkRep.R` and/or `quantify.R`

    -   enter the file paths to the `taskTable` and `quantTable` files  
    -   choose the type of peak fit (variable `fit_dim`) used in the
        analysis step

