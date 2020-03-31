#
# Quantification
#
# User configuration params ####
taskTable  = 'list_of_files_Francis_PP.csv'
quantTable = 'list_of_targets_Francis_quantification.csv'
# filter_results = TRUE
# fwhm_min = 0.5
# fwhm_max = 1.5
# area_min = 10
# const_fwhm = 0.7
# save_figures = TRUE

# Code Setup ####
options(warn=0)

# Install packages if necessary
## CRAN packages
libs <- c('xtable','mixtools','inlmisc',
          'rlist','repmis','assertive')
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(
      lib,
      dependencies = TRUE,
      repos = 'https://cran.univ-paris1.fr'
    )
  }
}
## Load packages
repmis::LoadandCite(libs) # ,file='../article/packages.bib')

# Set graphical params
## For PNG figures
gPars = list(
  cols     = rev(inlmisc::GetColors(8))[1:7],
  cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.2))[1:7],
  cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.5))[1:7],
  pty      = 's',
  mar      = c(4,4,3,1),
  mgp      = c(2,.75,0),
  tcl      = -0.5,
  lwd      = 4.0,
  cex      = 4.0,
  cex.leg  = 0.7,
  reso     = 1200  # (px) base resolution for png figs
)
## For local plots
gParsLoc = gPars
gParsLoc$cex = 1
gParsLoc$lwd = 2

# Expose gPars list
for (n in names(gPars))
  assign(n, rlist::list.extract(gPars, n))

# Define Data and Results repositories
dataRepo = '../data/'
figRepo  = '../results/figs/'
tabRepo  = '../results/tables/'

sink(file ='./sessionInfo.txt')
print(sessionInfo(), locale=FALSE)
sink()

# Functions ####

# Check sanity of parameters ####
assertive::assert_all_are_existing_files(dataRepo)
assertive::assert_all_are_existing_files(figRepo)
assertive::assert_all_are_existing_files(tabRepo)

file = paste0(dataRepo, taskTable)
assertive::assert_all_are_existing_files(file)

file = paste0(dataRepo, quantTable)
assertive::assert_all_are_existing_files(file)

# assertive::assert_is_numeric(fwhm_min)
# if(!assertive::is_positive(fwhm_min))
#   stop(paste0('Erreur: fwhm_min =',fwhm_min,' should be positive'))
#
# assertive::assert_is_numeric(fwhm_max)
# if(!assertive::is_positive(fwhm_max))
#   stop(paste0('Erreur: fwhm_max =',fwhm_max,' should be positive'))
#
# assertive::assert_is_numeric(area_min)
# if(!assertive::is_positive(area_min))
#   stop(paste0('Erreur: area_min =',area_min,' should be positive'))

# Get tasks list & results ####
file = paste0(dataRepo, taskTable)
Tasks = read.table(
  file = file,
  header = TRUE,
  sep = ',',
  check.names = FALSE,
  stringsAsFactors = FALSE
)
D = data.frame()
for(task in 1:nrow(Tasks)) {

  msTable = Tasks[task,'MS_file']
  CVTable = Tasks[task,'DMS_file']
  dilu    = Tasks[task,'dilu']
  # Build tag
  ## Extract date from CVTable
  date =
    strsplit(
      strsplit(
        CVTable,
        split = ' '
      )[[1]][2],
      split = '-'
    )[[1]][1]

  tag = paste0(
    date,'_',
    strsplit(msTable, split='\\.')[[1]][1]
  )

  M = read.csv(
    file = paste0(tabRepo,tag,'_results.csv'),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  M = cbind(M,dilu)
  D = rbind(D,M)
}

# Get targets ####
file = paste0(dataRepo, quantTable)
quant = read.table(
  file = file,
  header = TRUE,
  sep = ';',
  dec = '.',
  check.names = FALSE,
  fill = TRUE,
  stringsAsFactors = FALSE
)
targets = quant$Name

# Quantification ####
qres = data.frame(
  Name = NA,
  Int  = NA,
  Slo  = NA,
  Slo0 = NA,
  LOD  = NA)
pdf(
  file=paste0(figRepo,'quantif.pdf'),
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

  cr = cAA[selr]/cIS
  ratio = ratio[selr]
  dratio = dratio[selr]

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
  p0 = predict(reg0, ,
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

write.csv(qres[-1,],file=paste0(tabRepo,'quantif.csv'))

# END ####
