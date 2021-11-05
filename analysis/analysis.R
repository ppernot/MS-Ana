#----------------------------------------------------------
# User configuration params -------------------------------
#----------------------------------------------------------

ms_type   = c('esquire','fticr')[2]
taskTable = 'Test2/files_quantification_2018AA.csv'
tgTable   = 'Test2/targets_paper_renew.csv'

#----------------------------------------------------------
# Run analysis --------------------------------------------
#----------------------------------------------------------

resu0 = msAnaLib::dmsAnalysis(
  ms_type = ms_type,
  taskTable = taskTable,
  tgTable = tgTable,
  fit_dim = 0
)

resu1 = msAnaLib::dmsAnalysis(
  ms_type = ms_type,
  taskTable = taskTable,
  tgTable = tgTable,
  fit_dim = 1
)

resu2 = msAnaLib::dmsAnalysis(
  ms_type = ms_type,
  taskTable = taskTable,
  tgTable = tgTable,
  fit_dim = 2
)

# Results

par(mfrow=c(1,2))

plot(resu1$Area,resu0$Area); abline(a=0,b=1)
plot(resu1$Area,resu2$Area); abline(a=0,b=1)

plot(resu0$`m/z`,resu0$`FWHM_m/z`)
segments(
  resu0$`m/z`-2*resu0$`u_m/z`,resu0$`FWHM_m/z`,
  resu0$`m/z`+2*resu0$`u_m/z`,resu0$`FWHM_m/z`
)
segments(
  resu0$`m/z`,resu0$`FWHM_m/z`-2*resu0$`u_FWHM_m/z`,
  resu0$`m/z`,resu0$`FWHM_m/z`+2*resu0$`u_FWHM_m/z`
)
points(resu2$`m/z`,resu2$`FWHM_m/z`,col=2)
segments(
  resu2$`m/z`-2*resu2$`u_m/z`,resu2$`FWHM_m/z`,
  resu2$`m/z`+2*resu2$`u_m/z`,resu2$`FWHM_m/z`,
  col=2
)
segments(
  resu2$`m/z`,resu2$`FWHM_m/z`-2*resu2$`u_FWHM_m/z`,
  resu2$`m/z`,resu2$`FWHM_m/z`+2*resu2$`u_FWHM_m/z`,
  col=2
)

plot(resu1$CV,resu1$FWHM_CV)
points(resu2$CV,resu2$FWHM_CV,col=2)


