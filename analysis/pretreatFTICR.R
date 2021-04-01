#===============================================
# Pre-treatment of FT-ICR MS to reduce size
#
# 2021_03_29 [PP]
#
#===============================================
#


#=======================================================================
# User configuration params ####
#=======================================================================

## Data and results directories
origMsDir = 'Test_FTICR_2/FTICR'
compMsDir = 'Test_FTICR_2/FTICR_compressed'
ms_type   = 'fticr'

## Compress mode
compMode  = c('grid','targets','grid+targets')[3]

## Grid specifications
mzMin     = 70
mzMax     = 250
dmz       = 0.002

## Targets specifications
tgTable   = 'Test_FTICR_2/targets_list.csv'
dmzTarget = 1 # Delta m/z to keep around target

## Short run to check ?
test      = FALSE

#=======================================================================
# DO NOT MODIFY BELOW...
#=======================================================================

# Load packages and functions ####
source('functions.R')

# Checks and Get ms list ####
assertive::assert_all_are_existing_files(dataRepo)

origMs = file.path(dataRepo, origMsDir)
assertive::assert_all_are_existing_files(origMs)

msFiles = list.files(origMs) # All files in source Dir
if (length(msFiles) == 0) {
  message('>>> No data files to process !!!')
  stop(call. = FALSE)
}

compMs = file.path(dataRepo, compMsDir)
if (!dir.exists(compMs))
  dir.create(compMs)

targets = NA
if(grepl('targets',compMode)) {
  file = file.path(dataRepo, tgTable)
  assertive::assert_all_are_existing_files(file)
  targets = readTargetsFile(file)[,'m/z_ref']
}

# Loop over files ####
for (ms in msFiles)
  compressMS(
    file_in   = file.path(origMs, ms),
    file_out  = file.path(compMs, ms),
    ms_type   = ms_type,
    mzMin     = mzMin,
    mzMax     = mzMax,
    dmz       = dmz,
    test      = test,
    compMode  = compMode,
    targets   = targets,
    dmzTarget = dmzTarget
  )

# END ####
