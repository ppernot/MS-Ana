#
# Check repeatability of analysis
#
# User configuration params ####
taskTable = 'list_of_files_Viet_total.csv'
tgTable   = 'list_of_targets_Plasma_Viet.csv'
quantTable = 'list_of_targets_Quantification.csv'
filter_results = TRUE
fwhm_min = 0.5
fwhm_max = 1.5
area_min = 10
const_fwhm = 0.7
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

file = paste0(dataRepo, tgTable)
assertive::assert_all_are_existing_files(file)

file = paste0(dataRepo, taskTable)
assertive::assert_all_are_existing_files(file)

file = paste0(dataRepo, quantTable)
assertive::assert_all_are_existing_files(file)

assertive::assert_is_numeric(fwhm_min)
if(!assertive::is_positive(fwhm_min))
  stop(paste0('Erreur: fwhm_min =',fwhm_min,' should be positive'))

assertive::assert_is_numeric(fwhm_max)
if(!assertive::is_positive(fwhm_max))
  stop(paste0('Erreur: fwhm_max =',fwhm_max,' should be positive'))

assertive::assert_is_numeric(area_min)
if(!assertive::is_positive(area_min))
  stop(paste0('Erreur: area_min =',area_min,' should be positive'))

# Get targets ####
file = paste0(dataRepo, tgTable)
targets = read.table(
  file = file,
  header = TRUE,
  sep = ';',
  dec = '.',
  check.names = FALSE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

# Get tasks list ####
file = paste0(dataRepo, taskTable)
Tasks = read.table(
  file = file,
  header = TRUE,
  sep = ',',
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check that files exist before proceeding
files = paste0(dataRepo,Tasks[,'MS_file'])
assertive::assert_all_are_existing_files(files)

files = paste0(dataRepo,Tasks[,'DMS_file'])
assertive::assert_all_are_existing_files(files)

# Ckeck results files ####
D = data.frame()
for(task in 1:nrow(Tasks)) {

  msTable = Tasks[task,'MS_file']
  CVTable = Tasks[task,'DMS_file']

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
  conc = strsplit(msTable, split='_')[[1]][2]
  conc = substr(conc,2,nchar(conc))
  conc = 1/as.numeric(conc)
  M = read.csv(
    file = paste0(tabRepo,tag,'_results.csv'),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  M = cbind(M,conc)
  D = rbind(D,M)
}

for(it in 1:nrow(targets)) {
  name = targets[it, 'Name']
  sel = D$Name == name

  par(mfrow = c(1,3), mar = mar)

  # CV vs. conc
  x = D[sel, 'conc']
  y = D[sel, 'CV']
  dy = 2 * D[sel, 'u_CV']
  y0 = D[sel, 'CV_ref']
  plot(
    x,
    y,
    pch = 16,
    col = cols[5],
    log = 'x',
    xlab = 'conc',
    ylab = 'CV',
    main = name
  )
  segments(x, y - dy, x, y + dy, col = cols[5])
  lines(x, y0, lwd = 2, col = cols[3])
  text(min(x),y0[1],'CV_ref',adj=0,cex=0.75)
  cm = sort(unique(x))
  ym = c()
  for (i in 1:length(cm))
    ym[i] = mean(y[x == cm[i]], na.rm = TRUE)
  lines(cm, ym, lty = 2, col = cols[2])
  grid()

  # FWHM vs. conc
  x = D[sel, 'conc']
  y = D[sel, 'FWHM']
  dy = 2 * D[sel, 'u_FWHM']
  plot(
    x,
    y,
    pch = 16,
    col = cols[5],
    log = 'x',
    xlab = 'conc',
    ylab = 'FWHM',
    main = name,
    ylim = c(0.5, 0.9)
  )
  segments(x, y - dy, x, y + dy, col = cols[5])
  rect(0.8*min(x),0.8*const_fwhm,
       1.2*max(x),1.2*const_fwhm,
       col = cols_tr[4], border=NA)
  grid()

  # Area vs. conc
  x = D[sel, 'conc']
  y = D[sel, 'Area']
  dy = 2 * D[sel, 'u_Area']
  plot(
    x,
    y,
    pch = 16,
    col = cols[5],
    log = '',
    xlab = 'conc',
    ylab = 'Area',
    main = name
  )
  segments(x, y - dy, x, y + dy, col = cols[5])
  abline(h=area_min, lwd= 2, col=cols[2])
  io = order(x)
  xo = x[io]
  yo = y[io]
  reg = lm(yo~xo, weights = 1/(dy[io]/2)^2)
  p = predict(reg, interval = 'conf')
  matlines(xo[!is.na(yo)],p,
           col = cols[4],
           lty=c(1,2,2))
  grid()
}

# Quantification ####
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
qres = data.frame(
  Name = NA,
  Int  = NA,
  Slo  = NA,
  LOD  = NA)
par(mfrow=c(1,1),mar = mar)
for(it in 1:length(targets)) {
  AA = targets[it]
  IS = quant[it, 'IS']

  selAA = D$Name == AA
  if(sum(selAA) == 0)
    next

  diluAA  = D[selAA, 'conc']
  aireAA  = D[selAA, 'Area']
  daireAA = D[selAA, 'u_Area']

  selIS = D$Name == IS

  diluIS  = D[selIS, 'conc']
  aireIS  = D[selIS, 'Area']
  daireIS = D[selIS, 'u_Area']

  if (any(diluAA != diluIS))
    stop('Pb...')

  ratio = aireAA / aireIS
  dratio = ratio * sqrt(
      daireAA^2/aireAA^2 +
      daireIS^2/aireIS^2
  )
  selr = !is.na(ratio)
  dilur = diluAA[selr]
  ratio = ratio[selr]
  dratio = dratio[selr]

  plot(
    dilur,
    ratio,
    pch = 16,
    col = cols[5],
    log = '',
    xlab = 'Conc',
    ylab = 'AA/IS',
    main = paste(AA,'/',IS)
  )
  segments(dilur, ratio - 2 * dratio,
           dilur, ratio + 2 * dratio,
           col = cols[5])

  io = order(dilur)
  xo = dilur[io]
  yo = ratio[io]
  reg = lm(yo~xo, weights = 1/(dratio[io]/2)^2)
  p = predict(reg, interval = 'conf')
  matlines(xo[!is.na(yo)],p,
           col = cols[4],
           lty=c(1,2,2))
  grid()

  # LOD

  Sxy = sqrt(sum(residuals(reg)^2)/length(xo-2))
  Clod = 3*Sxy/coefficients(reg)[2]
  abline(v=Clod , col=cols[2], lty = 3)
  mtext('LOD',side=3,col=cols[2],at=Clod)

  qres = rbind(
    qres,
    c(AA,
      signif(coefficients(reg),3),
      signif(Clod,3)
    )
  )

}
print(qres[-1,])
# END ####
