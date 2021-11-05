#
# Quantification
#
#===============================================
# 2020_07_17 [PP]
# - Adapt to new naming conventions in analysis.R
# 2020_07_20 [PP]
# - replaced '=' by '_' in userTag (Windows pb.)
# 2020_09_30 [PP]
# - filtered out results and figs with NA problems
#===============================================
#
# Load packages and functions ####
source('functions.R')

# User configuration params ####
taskTable  = 'files_quantification_2019July10.csv'
quantTable = 'targets_paper_quantification.csv'

fit_dim = 2
userTag = paste0('fit_dim_',fit_dim)

dataTag = format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

# Check sanity of parameters ####
assertive::assert_all_are_existing_files(dataRepo)
assertive::assert_all_are_existing_files(figRepo)
assertive::assert_all_are_existing_files(tabRepo)

file = paste0(dataRepo, taskTable)
assertive::assert_all_are_existing_files(file)

file = paste0(dataRepo, quantTable)
assertive::assert_all_are_existing_files(file)

# Get tasks list & gather results ####
Tasks = readTasksFile(paste0(dataRepo, taskTable))
D     = gatherResults(Tasks, tabRepo, userTag)

# Get targets ####
quant   = readTargetsFile(paste0(dataRepo, quantTable))
targets = quant$Name

# Quantification ####
qres = data.frame(
  Name = NA,
  Int  = NA,
  Slo  = NA,
  Slo0 = NA,
  LOD  = NA)

pdf(
  file=paste0(
    figRepo,
    dataTag,
    ifelse(userTag == '', '', '_'),
    userTag,
    '_quantif.pdf'),
  width = 7,
  height = 10)
par(mfrow=c(3,2),
    mar = mar,
    pty = pty,
    mgp = mgp,
    tcl = tcl)
for(it in 1:length(targets)) {
  AA = targets[it]
  IS = quant[it, 'IS']

  selAA = D$Name == AA
  if(sum(selAA) == 0)
    next

  diluAA  = D[selAA, 'dilu']
  aireAA  = D[selAA, 'Area']
  daireAA = D[selAA, 'u_Area']

  selIS = D$Name == IS

  diluIS  = D[selIS, 'dilu']
  aireIS  = D[selIS, 'Area']
  daireIS = D[selIS, 'u_Area']

  if (any(diluAA != diluIS))
    stop('Pb...')

  cAA = quant[it, 'CAA_Plasma'] +
    quant[it, 'CAA_ref'] * diluAA
  cIS = quant[it, 'CIS_ref']

  ratio = aireAA / aireIS
  dratio = ratio * sqrt(
      daireAA^2/aireAA^2 +
      daireIS^2/aireIS^2
  )
  selr = !is.na(ratio)
  if(sum(selr) == 0)
    next

  cr = cAA[selr]/cIS
  ratio = ratio[selr]
  dratio = dratio[selr]

  if(!is.finite(cr) | !is.finite(ratio))
    next

  plot(
    cr,
    ratio,
    pch = 16,
    col = cols[5],
    log = '',
    xlab = 'cAA/cIS', xlim = c(0,max(cr)),
    ylab = 'AA/IS', ylim = c(0,max(ratio)),
    main = paste(AA,'/',IS)
  )
  segments(cr, ratio - 2 * dratio,
           cr, ratio + 2 * dratio,
           col = cols[5])

  # Linear regression
  io = order(cr)
  xo = cr[io]
  yo = ratio[io]
  reg = lm(yo~xo, weights = 1/(dratio[io]/2)^2)
  x1 = c(0,xo)
  p = predict(reg,
              newdata = list(xo = x1),
              interval = 'conf')
  matlines(x1[!is.na(p[,1])],p,
           col = cols[4],
           lty=c(1,2,2))

  # LInear regression with intercept at 0
  reg0 = lm(yo~0+xo, weights = 1/(dratio[io]/2)^2)
  p0 = predict(reg0,
               newdata = list(xo = x1),
               interval = 'conf')
  matlines(x1[!is.na(p0[,1])],p0,
           col = cols[2],
           lty=c(1,2,2))
  grid()

  # LOD from https://sites.chem.utoronto.ca/chemistry/coursenotes/analsci/stats/LimDetect.html
  Sxy = sqrt(sum(residuals(reg)^2)/length(xo-2))
  Clod = 3*Sxy/coefficients(reg)[2]
  abline(v=Clod , col=cols[3], lty = 2)
  mtext('LOD',side=3,col=cols[3],at=Clod, cex=0.75)

  # Accumulate results
  qres = rbind(
    qres,
    c(AA,
      signif(coefficients(reg),3),
      signif(coefficients(reg0),3),
      signif(Clod,3)
    )
  )
}
dev.off()

write.csv(
  qres[-1,],
  file = paste0(
    tabRepo,
    dataTag,
    ifelse(userTag == '', '', '_'),
    userTag,
    '_quantif.csv'))

# END ####
