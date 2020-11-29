## ---- eval = FALSE------------------------------------------------------------
#  install.packages('ADMUR')

## ---- eval = FALSE------------------------------------------------------------
#  install.packages('devtools')
#  library(devtools)
#  install_github('UCL/ADMUR')

## ---- message = FALSE---------------------------------------------------------
library(ADMUR)

## ---- eval = FALSE------------------------------------------------------------
#  help(ADMUR)
#  help(SAAD)

## ---- eval = TRUE-------------------------------------------------------------
SAAD[1:5,1:8]

## ---- eval = TRUE-------------------------------------------------------------
citation('ADMUR')

## ---- eval = TRUE, fig.height = 3, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
data <- data.frame( age=c(6562,7144), sd=c(44,51) )
x <- summedCalibratorWrapper(data)

## ---- eval = TRUE, fig.height = 3, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
data <- data.frame( age=c(6562,7144), sd=c(44,51), datingType=c('14C','TL') )
x <- summedCalibratorWrapper(data=data, calcurve=shcal20)

## ---- eval = TRUE, fig.height = 3, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
data <- data.frame( age = c(9144), sd=c(151) )
CalArray <- makeCalArray( calcurve=intcal20, calrange=c(8000,13000) )
cal <- summedCalibrator(data, CalArray)
plotPD(cal)

## ---- eval = TRUE, fig.height = 5, fig.width=7, fig.align = "center", dev='jpeg'----
x <- makeCalArray( calcurve=shcal20, calrange=c(5500,6000), inc=1 )
plotCalArray(x)

## ---- eval = TRUE-------------------------------------------------------------
data <- subset( SAAD, site %in% c('Carrizal','Pacopampa') )
data[,2:7]

## ---- eval = TRUE, fig.height = 3, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
CalArray <- makeCalArray( calcurve=shcal20, calrange=c(2000,6000) )
x <- phaseCalibrator(data=data, CalArray=CalArray)
plotPD(x)

## ---- eval = TRUE-------------------------------------------------------------
SPD <- as.data.frame( rowSums(x) )

# normalise
SPD <- SPD/( sum(SPD) * CalArray$inc )

## ---- eval = TRUE, fig.height = 3, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
SPD <- summedPhaseCalibrator( data=data, calcurve=shcal20, calrange=c(2000,6000) )
plotPD(SPD)

## ---- eval = TRUE-------------------------------------------------------------
set.seed(12345) 
N <- 350

# randomly sample calendar dates from the toy model
cal <- simulateCalendarDates(toy, N)

# Convert to 14C dates. 
age <- uncalibrateCalendarDates(cal, shcal20)
data <- data.frame(age = age, sd = 50, phase = 1:N, datingType = '14C')

# Calibrate each phase, taking care to restrict to the modelled date range with 'remove.external'
CalArray <- makeCalArray(shcal20, calrange = range(toy$year))
PD <- phaseCalibrator(data, CalArray, remove.external = TRUE)

## ---- eval = TRUE-------------------------------------------------------------
print( ncol(PD) )

## ---- eval = TRUE-------------------------------------------------------------
loglik(PD=PD, model=toy)

## ---- eval = TRUE-------------------------------------------------------------
uniform.model <- convertPars(pars=NULL, years=5500:7500, type='uniform')
loglik(PD=PD, model=uniform.model)

## ---- eval = TRUE-------------------------------------------------------------
exp( loglik(PD=PD, model=toy) - loglik(PD=PD, model=uniform.model) )

## ---- eval = TRUE-------------------------------------------------------------
set.seed(12345)
convertPars( pars=runif(11), years=5500:7500, type='CPL' )

## ---- eval = FALSE------------------------------------------------------------
#  library(DEoptimR)
#  best <- JDEoptim(lower = rep(0,5),
#  	upper = rep(1,5),
#  	fn = objectiveFunction,
#  	PDarray = PD,
#  	type = 'CPL',
#  	NP = 100,
#  	Trace = TRUE)

## ---- echo = FALSE------------------------------------------------------------
load('vignette.3CPL.JDEoptim.best.RData')

## ---- eval = TRUE, fig.height = 4, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
CPL <- convertPars(pars=best$par, years=5500:7500, type='CPL')
SPD <- summedPhaseCalibrator( data=data, calcurve=shcal20, calrange=c(5500,7500) )
plotPD(SPD)
lines(CPL$year, CPL$pdf, lwd=2, col='firebrick')
legend(x=6300, y=max(CPL$pdf), cex=0.7, lwd=2, col='firebrick', bty='n', legend='best fitted 3-CPL')
text(x=CPL$year, y=CPL$pdf, pos=3, labels=c('H1','H2','H3','H4'))

## ---- eval = FALSE------------------------------------------------------------
#  chain <- mcmc(PDarray=PD, startPars=best$par, type='CPL', N=100000, burn=2000, thin=5, jumps =0.025)

## ---- eval = FALSE------------------------------------------------------------
#  print(chain$acceptance.ratio)
#  par(mfrow=c(3,2), mar=c(4,3,3,1))
#  col <- 'steelblue'
#  for(n in 1:5){
#  	plot(chain$all.pars[,n], type='l', ylim=c(0,1), col=col, xlab='', ylab='', main=paste('par',n))
#  	}

## ---- eval = FALSE------------------------------------------------------------
#  hinges <- convertPars(pars=chain$res, years=5500:7500, type='CPL')
#  par(mfrow=c(3,2), mar=c(4,3,3,1))
#  c1 <- 'steelblue'
#  c2 <- 'firebrick'
#  lwd <- 3
#  pdf.brk <- seq(0,0.0015, length.out=40)
#  yr.brk <- seq(5500,7500,length.out=40)
#  names <- c('Date of H2','Date of H3','PD of H1','PD of H2','PD of H3','PD of H4')
#  hist(hinges$yr2,border=c1,breaks=yr.brk, main=names[1], xlab='');abline(v=CPL$year[2],col=c2,lwd=lwd)
#  hist(hinges$yr3, border=c1,breaks=yr.brk, main=names[2], xlab='');abline(v=CPL$year[3],col=c2,lwd=lwd)
#  hist(hinges$pdf1, border=c1,breaks=pdf.brk, main=names[3], xlab='');abline(v=CPL$pdf[1],col=c2,lwd=lwd)
#  hist(hinges$pdf2, border=c1,breaks=pdf.brk, main=names[4], xlab='');abline(v=CPL$pdf[2],col=c2,lwd=lwd)
#  hist(hinges$pdf3, border=c1,breaks=pdf.brk, main=names[5], xlab='');abline(v=CPL$pdf[3],col=c2,lwd=lwd)
#  hist(hinges$pdf4, border=c1,breaks=pdf.brk, main=names[6], xlab='');abline(v=CPL$pdf[4],col=c2,lwd=lwd)

## ---- eval = FALSE------------------------------------------------------------
#  require(scales)
#  par( mfrow=c(1,2) , mar=c(4,4,1.5,2), cex=0.7 )
#  plot(hinges$yr2, hinges$pdf2, pch=16, col=alpha(1,0.02), ylim=c(0,0.0005))
#  points(CPL$year[2], CPL$pdf[2], col='red', pch=16, cex=1.2)
#  plot(hinges$yr3, hinges$pdf3, pch=16, col=alpha(1,0.02), ylim=c(0,0.0015))
#  points(CPL$year[3], CPL$pdf[3], col='red', pch=16, cex=1.2)

## ---- eval = FALSE------------------------------------------------------------
#  plot(NULL, xlim=c(7500,5500),ylim=c(0,0.0011), xlab='calBP', ylab='PD', cex=0.7)
#  for(n in 1:nrow(hinges)){
#  	x <- c(hinges$yr1[n], hinges$yr2[n], hinges$yr3[n], hinges$yr4[n])
#  	y <- c(hinges$pdf1[n], hinges$pdf2[n], hinges$pdf3[n], hinges$pdf4[n])
#  	lines( x, y, col=alpha(1,0.005) )
#  	}
#  lines(x=CPL$year, y=CPL$pdf, lwd=2, col=c2)

## ---- eval = TRUE, fig.height = 4, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
N <- 1000
x <- cbind(rep(5100,N),rep(5000,N))
y <- cbind(seq(1,100,length.out=N),seq(100,1,length.out=N))
conventional <- 100 * exp(log(y[,2]/y[,1])/((x[,1]-x[,2])/25))-100
relative <- relativeRate(x,y)
plot(conventional, relative, type='l')

## ---- eval = FALSE------------------------------------------------------------
#  # CPL parameters must be between 0 and 1, and an odd length.
#  CPL.1 <- JDEoptim(lower=0, upper=1, fn=objectiveFunction, PDarray=PD, type='CPL', NP=20)
#  CPL.2 <- JDEoptim(lower=rep(0,3), upper=rep(1,3), fn=objectiveFunction, PDarray=PD, type='CPL', NP=60)
#  CPL.3 <- JDEoptim(lower=rep(0,5), upper=rep(1,5), fn=objectiveFunction, PDarray=PD, type='CPL', NP=100)
#  CPL.4 <- JDEoptim(lower=rep(0,7), upper=rep(1,7), fn=objectiveFunction, PDarray=PD, type='CPL', NP=140)
#  
#  # exponential has a single parameter, which can be negative (decay).
#  exp <- JDEoptim(lower=-0.01, upper=0.01, fn=objectiveFunction, PDarray=PD, type='exp', NP=20)
#  
#  # uniform has no parameters so a search is not required.
#  uniform <- objectiveFunction(NULL, PD, type='uniform')

## ---- echo = FALSE------------------------------------------------------------
load('vignette.model.comparison.RData')

## ---- eval = TRUE-------------------------------------------------------------
# likelihoods
data.frame(L1= -CPL.1$value,
	L2= -CPL.2$value,
	L3= -CPL.3$value,
	L4= -CPL.4$value,
	Lexp= -exp$value,
	Lunif= -uniform)
BIC.1 <- 1*log(303) - 2*(-CPL.1$value)
BIC.2 <- 3*log(303) - 2*(-CPL.2$value)
BIC.3 <- 5*log(303) - 2*(-CPL.3$value)
BIC.4 <- 7*log(303) - 2*(-CPL.4$value)
BIC.exp <- 1*log(303) - 2*(-exp$value)
BIC.uniform <- 0 - 2*(-uniform)
data.frame(BIC.1,BIC.2,BIC.3,BIC.4,BIC.exp,BIC.uniform)

## ---- eval = TRUE, fig.height = 4, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
# convert parameters to model PDs
CPL1 <- convertPars(pars=CPL.1$par, years=5500:7500, type='CPL') 
CPL2 <- convertPars(pars=CPL.2$par, years=5500:7500, type='CPL')  
CPL3 <- convertPars(pars=CPL.3$par, years=5500:7500, type='CPL')  
CPL4 <- convertPars(pars=CPL.4$par, years=5500:7500, type='CPL')  
EXP <- convertPars(pars=exp$par, years=5500:7500, type='exp')  

# Plot SPD and five competing models:
plotPD(SPD)
cols <- c('firebrick','orchid2','coral2','steelblue','goldenrod3')
lines(CPL1$year, CPL1$pdf, col=cols[1], lwd=2)
lines(CPL2$year, CPL2$pdf, col=cols[2], lwd=2)
lines(CPL3$year, CPL3$pdf, col=cols[3], lwd=2)
lines(CPL4$year, CPL4$pdf, col=cols[4], lwd=2)
lines(EXP$year, EXP$pdf, col=cols[5], lwd=2)
legend <- c('1-CPL','2-CPL','3-CPL','4-CPL','exponential')
legend(x=6300, y=max(CPL$pdf), cex=0.7, lwd=2, col=cols, bty='n', legend=legend)

## ---- eval = FALSE------------------------------------------------------------
#  summary <- SPDsimulationTest(data, calcurve=shcal20, calrange=c(5500,7500), pars=CPL.3$par, type='CPL')

## ---- echo = FALSE------------------------------------------------------------
load('vignette.3CPL.SPDsimulationTest.RData')

## ---- eval = TRUE, fig.height = 5, fig.width=7, fig.align = "center", dev='svg', warning=FALSE----
print(summary$pvalue)
hist(summary$simulated.stat, main='Summary statistic', xlab='')
abline(v=summary$observed.stat, col='red')
legend(0.3,6000, bty='n', lwd=c(1,3), col=c('red','grey'), legend=c('observed','simulated'))

## ---- eval = FALSE------------------------------------------------------------
#  summary <- SPDsimulationTest(data, calcurve=shcal20, calrange=c(5500,7500), pars=exp$par, type='exp')

## ---- echo = FALSE------------------------------------------------------------
load('vignette.exp.SPDsimulationTest.RData')

## ---- eval = TRUE, fig.height = 4, fig.width=7, fig.align = "center", dev='svg', warning=FALSE, message=FALSE----
plotSimulationSummary(summary, legend.y=0.0012)

## ---- eval = FALSE------------------------------------------------------------
#  # generate SPDs
#  spd1 <- summedCalibrator(data1, CalArray, normalise='full')
#  spd2 <- summedCalibrator(data2, CalArray, normalise='full')
#  spd3 <- summedCalibrator(data3, CalArray, normalise='full')
#  
#  # calibrate phases
#  CalArray <- makeCalArray(intcal20, calrange = c(5500,7500), inc = 5)
#  PD1 <- phaseCalibrator(data1, CalArray, remove.external = TRUE)
#  PD2 <- phaseCalibrator(data2, CalArray, remove.external = TRUE)
#  PD3 <- phaseCalibrator(data3, CalArray, remove.external = TRUE)
#  
#  # maximum likelihood search, fitting various models to various datasets
#  norm <- JDEoptim(lower=c(5500,1), upper=c(7500,5000), fn=objectiveFunction,
#  		PDarray=PD1, type='norm', NP=40, trace=T)
#  cauchy <- JDEoptim(lower=c(5500,1), upper=c(7500,5000), fn=objectiveFunction,
#  		PDarray=PD1, type='cauchy', NP=40, trace=T)
#  sine <- JDEoptim(lower=c(0,0,0), upper=c(1/1000,2*pi,1), fn=objectiveFunction,
#  		PDarray=PD2, type='sine', NP=60, trace=T)
#  logistic <- JDEoptim(lower=c(0,5000), upper=c(1,8000), fn=objectiveFunction,
#  		PDarray=PD3, type='logistic', NP=40, trace=T)
#  exp <- JDEoptim(lower=c(0), upper=c(1), fn=objectiveFunction,
#  		PDarray=PD3, type='exp', NP=20, trace=T)

## ---- eval = FALSE------------------------------------------------------------
#  # convert parameters to model PDs
#  mod.norm <- convertPars(pars=norm$par, years=5500:7500, type='norm')
#  mod.cauchy <- convertPars(pars=cauchy$par, years=5500:7500, type='cauchy')
#  mod.sine <- convertPars(pars=sine$par, years=5500:7500, type='sine')
#  mod.uniform <- convertPars(pars=NULL, years=5500:7500, type='uniform')
#  mod.logistic <- convertPars(pars=logistic$par, years=5500:7500, type='logistic')
#  mod.exp <- convertPars(pars=exp$par, years=5500:7500, type='exp')
#  
#  # Plot SPDs and various fitted models:
#  par(mfrow=c(3,1), mar=c(4,4,1,1))
#  cols <- c('steelblue','firebrick')
#  
#  plotPD(spd1)
#  lines(mod.norm, col=cols[1], lwd=2)
#  lines(mod.cauchy, col=cols[2], lwd=2)
#  legend(x=7500, y=max(spd1), lwd=2, col=cols, bty='n', legend=c('Gaussian','Cauchy'))
#  
#  plotPD(spd2)
#  lines(mod.sine, col=cols[1], lwd=2)
#  lines(mod.uniform, col=cols[2], lwd=2)
#  legend(x=7500, y=max(spd2),  lwd=2, col=cols, bty='n', legend=c('Sinewave','Uniform'))
#  
#  plotPD(spd3)
#  lines(mod.logistic, col=cols[1], lwd=2)
#  lines(mod.exp, col=cols[2], lwd=2)
#  legend(x=7500, y=max(spd3), lwd=2, col=cols, bty='n', legend=c('Logistic','Exponential'))

