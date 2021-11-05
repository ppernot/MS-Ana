# `pretreatFTICR.R` script

Raw FT-ICR DMS files are very large and are painful to analyse directly.
Two reduction methods can be used and combined:

1.  filtering the data to match as well as possible a predefined regular
    grid

2.  keeping only the data that match a list of targets

Combining both methods might enable to reduce the MS size by nearly 90%,
in the usual applications.

The script is controlled by the following parameters:

    #====================================================
    # User configuration params =========================
    #====================================================

    ## Data and results directories
    origMsDir = 'Test_FTICR_2/FTICR'
    compMsDir = 'Test_FTICR_2/FTICR_compressed'
    ms_type   = 'fticr'

    ## Compress mode
    compMode  = c('grid','targets','grid+targets')[3]

    ## Grid specifications
    mzMin     = 70
    mzMax     = 250
    dmz       = 0.001

    ## Targets specifications
    tgTable   = 'Test_FTICR_2/targets_list.csv'
    dmzTarget = 0.5 # Delta m/z to keep around target

    ## Short run to check ?
    test      = FALSE

where

-   `origMsDir`: (string) is the path to a directory containing the
    files to treat (all files within the directory will be treated).

-   `compMsDir`: (string) is the path to a directory where the
    compressed files will be stored.

-   `ms_type`: (string) is the type of MS to treat. Only `fticr` at the
    moment.

-   `compMode`: (string) compression modality, which can be `grid`,
    `targets` or their combination `grid+targets`.

-   `mzMin`, `mzMax` and `dmz`: (numbers) define the regular grid used
    for grid selection.

-   `tgTable`: (string) is the path to the targets file used for targets
    selection.

-   `dmzTarget`: (number) is the half-width of the m/z interval to keep
    around a target. The interval is centered on `ms_ref` vales in
    `tgTable`. The width should be large enough to not miss the peaks,
    but small enough for efficient compression…

-   `test`: (logical) test run only, to check if data are OK.

**Note**: it is best to store the spectra as compressed .gz files. The
disp space is well reduced and all scripts can handle them
transparently.

