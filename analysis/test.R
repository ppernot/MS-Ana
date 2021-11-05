#----------------------------------------------------------
# User configuration params -------------------------------
#----------------------------------------------------------

ms_type   = c('esquire','fticr')[1]
taskTable = 'testMeanCV/files_AA_20210614.csv'
tgTable   = 'testMeanCV/targets_10AA.csv'

#----------------------------------------------------------
# Run analysis --------------------------------------------
#----------------------------------------------------------

resu0 = msAnaLib::dmsAnalysis(
  ms_type = ms_type,
  taskTable = taskTable,
  tgTable = tgTable,
  fit_dim = 0
)

