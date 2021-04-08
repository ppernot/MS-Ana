# `quantify.R` script

This script is based on the same principle as `checkRep.R` (same input
files), but aims to estimate the quantification parameters, such as the
LOD.

## Control variables

The job is defined by a few parameters.

    taskTable  = 'files_quantification_2019July10.csv'
    quantTable = 'targets_paper_quantification.csv'

    fit_dim = 2
    userTag = paste0('fit_dim_',fit_dim)

-   `taskTable`: (string) file path to the list of experiments to be
    compared

-   `quantTable`: (string) list of the compounds for which the
    comparisons should be done

-   `fit_dim`: (integer) type of peak fit for which the comparisons
    should be done

-   `userTag`: (string) tag to differentiate the outputs. In the present
    case, one wants to compare the repeatability for different peak fit
    approaches.

## Outputs

### Figures

A PDF file is generated, containing the quantification plots for all
compounds.

The name of the file is a concatenation of the date, time, `userTag`,
and ‘\_quantif.pdf’

### Tables

A ‘.csv’ table containing the quantification results.

The name of the file is a concatenation of the date, time, `userTag`,
and ‘\_quantif.csv’

#### Structure

| Name | Int | Slo | Slo0 | LOD |
|------|-----|-----|------|-----|
|      |     |     |      |     |

where

-   `Name` is the name of the compound

-   `Int` is the value of the intercept

-   `Slo` is the value of the slope

-   `Slo0` is the value of the slope with null intercept

-   `LOD` is the estimated limit of detection

