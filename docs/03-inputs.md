# Input files

Three input files (formally named `taskTable`, `tgTable` and
`quantTable`) are used by the set of scripts :

| Script \\ Input | `taskTable` | `tgTable` | `quantTable` |
|-----------------|-------------|-----------|--------------|
| `analysis.R`    | X           | X         |              |
| `checkRep.R`    | X           |           | X            |
| `quantify.R`    | X           |           | X            |

## `taskTable`

This file defines the list of MS-DMS pairs to be analyzed.

It is a “comma” (,) delimited ‘.csv’ file. It can be edited using excel
or Rstudio (safer).

### Structure

| MS\_file                   | DMS\_file                         | t0   | CV0 | dilu | Path                                            |
|----------------------------|-----------------------------------|------|-----|------|-------------------------------------------------|
| C0\_AS\_DV-1800\_1.d.ascii | Fichier\_Dims 20190517-000000.txt | 0.08 | 6   | 0    | Esquire\_MSMS\_Data/2019\_A\_Voir/20190517\_AA/ |

Where:

-   `MS_file` is an ASCII file, extracted using DATAANALYSIS. So far,
    only the ESQUIRE data files extracted using the `profile` option can
    be handled. It is stored in `data` or a sub-folder of `data` defined
    by `Path`.

-   `DMS_file` is the corresponding DMS file. It is expected to be in
    the `data` folder.

-   `t0` and `CV0` are used to convert the ESQUIRE time *t* values into
    DMS *CV* values.

-   `dilu` was initially meant to be the dilution factor of the standard
    metabolites when spiked into a plasma (see checkRep). When you
    perform another type of experiment, you can use  
    `dilu` as an index to specify, for example, the flow-rate of the
    modifier, the day of the experiment, the set of samples…

-   `Path` allows you to organize your data within the `../data/`
    folder. Note that the DMS\_files must be in the `../data` folder. In
    the present example, only the MS file is expected to be found in the
    following folder:
    `../data/Esquire_MSMS_Data/2019_A_Voir/20190517_AA/`.

**Notes**

-   lines starting with “\#” are not processed and treated as comments
    lines

-   the date extracted from `DMS_file` (here ‘20190517’)  
    and the root of the `MS_file` name (here ‘C0\_AS\_DV-1800\_1’) are
    combined to tag the output figures and tables (*e.g.*, ‘20190517\_
    C0\_AS\_DV-1800\_1.results’)

## `tgTable`

This file contains the list of compounds to be analyzed in each MS/DMS
data set.

It is a “commma” (,) delimited ‘.csv’ file. It can be edited using excel
or Rstudio (safer).

### Structure

| Name      | Comments  | m/z\_ref | CV\_ref |
|-----------|-----------|----------|---------|
| \# Gly-AA | C2H5NO2H  | 76       | -10.7   |
| Ala-AA    | 90.054955 | 90.1     | -7.6    |

where:

-   `Name` is the given name of a metabolite,

-   `Comments` is presently not used

-   `m/z_ref` is the expected (approximate) *m/z* value

-   `CV_ref` is the expected *CV* value (it can be omitted)

**Notes**

-   lines starting with “\#” will be considered as comment lines. In the
    present example, Glycine will not be analyzed.

## `quantTable`

This file contains the list of compounds and internal references used
for quantification.

It is a “comma” (,) delimited ‘.csv’ file. It can be edited using excel
or Rstudio (safer).

### Structure

| Name      | IS       | CAA\_Plasma | CAA\_ref | CIS\_ref |
|-----------|----------|-------------|----------|----------|
| \# Gly-AA | Gly-13C2 | 11.75       | 1750     | 26.7     |
| Ala-AA    | Ala-13C2 | 16.15       | 1750     | 20       |

where

-   `Name` is the name of the compound, as used in `tgTable`

-   `IS` is the name of the internal spiking compound. It should also be
    present in `tgTable`

-   `CAA_Plasma` **???**

-   `CAA_ref` **???**

-   `CIS_ref` **???**

