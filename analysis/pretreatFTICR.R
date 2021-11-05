#=======================================================================
# User configuration params ####
#=======================================================================

## Data and results directories
origMsDir = 'Test_FTICR_2/FTICR'
compMsDir = 'Test_FTICR_2/FTICR_compressed'

## Compress mode
compMode  = 'targets'
tgTable   = 'Test_FTICR_2/targets_list.csv'
dmzTarget = 0.5 # Delta m/z to keep around target

msAnaLib::compressAllMS(
  ms_type   = 'fticr',
  origMsDir = origMsDir,
  compMsDir = compMsDir,
  compMode  = compMode,
  tgTable   = tgTable,
  dmzTarget = dmzTarget,
  test      = FALSE
)
