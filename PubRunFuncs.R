



#### RAREFACTION FUNCTIONS (JUST HYPERDOM 50%)
library(sads) # sads used for alpha fitting
source("atdn-estimates-2019-master/functions.R") # source all functions from this file: nDom, ls.rad, etc
# define negative version of LogSeriesLogL log-likelihood function to input into nlm
LogSeriesLogLNeg <- function(gamma, Nbar, S) {
  #   gamma   Parameter of the distribution
  #   Nbar    Observed # of individuals per species
  #   S       Number of species with at least 1 individual
  -S*log(gamma) - (S*Nbar*log(1-exp(-1/gamma)))
}
# define function to subsample n rows (plots) from DF and compute Fisher's Alpha, Richard's Gamma, and hyperdom numbers from SAD (INCLUDING indets)
# input matrix of plot SADs, by default samples all rows, and sample is done without replacement (choose replacement = T to alter)
subsamplestats <- function(df,n=nrow(df),replacement=F) {
  df0 <- df[sample(1:nrow(df),n,replace=replacement),] # subsample n rows from df without replacement
  df1 <- df0[, names(df0)!='Indet'] # remove indets from fitting
  csum <- colSums(df1) # sum SAD matrix to SAD vector across all sampled plots
  csum1 <- csum[csum!=0] # remove 0 entries
  y <- cumsum(sort(csum1,decreasing = T)) # cumulatively sum abundances sorted by decreasing magnitude
  stems <- sum(colSums(df0)) # compute number of stems in subsample, including indets
  h <- min(which (y >= stems*0.5)) # compute number of hyperdominants
  s <- length(csum1) # compute number of species in subsample
  perc <- round(h*100/s,2) # compute percentage hyperdominants
  alphafit <- fitsad(csum1,'ls') # fisher's alpha fit from 'sads' package
  alphaMLE <- as.numeric(alphafit@coef) # Fisher's alpha MLE
  alphaSE <- as.numeric(sqrt(alphafit@vcov)) # sqrt of covariance matrix of logseries fitted alpha for SD of LS alpha fit
  nbar <- mean(csum1) # compute mean number of individuals per species in subsample
  gammafit <- nlm(LogSeriesLogLNeg,0.1,Nbar=nbar,S=s,hessian = T) # Richard's gamma fit using nlm minimisiation of -ve log-likelihood <-> maximisation of +ve
  gammaMLE <- gammafit$estimate # MLE fitted value of gamma
  gammaSE <- 1/sqrt(gammafit$hessian) # SE error associated with gamma MLE
  return(c(n,s,h,perc,alphaMLE,alphaSE,gammaMLE,gammaSE,stems))
}
# define function to calculate SD of vector just like rowMeans in vectorised fashion
rowSDs <- function(x, na.rm=F) {
  sqrt(rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1))
}
# define function to rarefy SAD matrix, calculating TS, H#, H%, Alpha, Gamma values and their SDs for each subsample
# input matrix of plot SADs, by default iterations, and sample is done without replacement (choose replacement = T to alter)
Rarefaction <- function(df,repeats = 100,replacement=F){
  output <- matrix(ncol=12, nrow=nrow(df)) # initiate empty matrix of correct dimensions to be filled by operations below
  for (j in seq(1,nrow(df),by=1)){
    reps <- replicate(repeats,subsamplestats(df,n=j,replacement = replacement)) # replicate function call to masterloop function above with n = j plots 100 times
    SDs <- rowSDs(reps[2:4,]) # take SD across iterations for total species, H#, and H%
    output[j,]<-c(rowMeans(reps),SDs) # take the mean across all iterations of subsample for each measurement aka take the mean of the rows of iterations
    cat(paste0(round(j / length(df) * 100), '%'))
  }
  output1 <-data.frame(output) # convert to DF
  output2 <- output1[complete.cases(output1),] # remove empty rows
  names(output2) <- c('Plots','TS','H#','H%','AlphaMLE','AlphaSE','GammaMLE','GammaSE','Stems','TSSD','H#SD','H%SD') # give columns meaningful names
  return(output2)
}
# define function to replace data in specified set of columns in one DF with data from corresponding columns with same name in other DF
# input original DF to be changed (here WITHOUT replacement) and imputed DF to replace values from original (here WITH replacement)
replace_imputed <- function(original, imputed, toChange = c('TSSD','H#SD','H%SD')){
  toChange <- toChange
  for(i in 1:length(toChange)){
    original[toChange[i]] <- imputed[toChange[i]]
  }
  return(original)
}
# define function to rarefy SAD without and with replacement, calculate TS, H#, H%, Alpha, Gamma values and their SDs for each subsample
# and output 
completeRare <- function(df,repeats = 100){
  without <- Rarefaction(df,repeats=repeats,replacement = F)
  with <- Rarefaction(df,repeats=repeats,replacement = T)
  final <- replace_imputed(without,with)
  return(final)
}

#### RAREFACTION FUNCTIONS (HYPERDOM 10-90%)
# define function to subsample n rows (plots) from DF and compute Fisher's Alpha, Richard's Gamma, and hyperdom numbers from SAD (INCLUDING indets)
# input matrix of plot SADs, by default samples all rows, and sample is done without replacement (choose replacement = T to alter)
subsampleperc <- function(df,n=nrow(df),replacement=F) {
  df0 <- df[sample(1:nrow(df),n,replace=replacement),] # subsample n rows from df without replacement
  df1 <- df0[, names(df0)!='Indet'] # remove indets from fitting
  csum <- colSums(df1) # sum SAD matrix to SAD vector across all sampled plots
  csum1 <- csum[csum!=0] # remove 0 entries
  y <- cumsum(sort(csum1,decreasing = T)) # cumulatively sum abundances sorted by decreasing magnitude
  stems <- sum(colSums(df0)) # compute number of stems in subsample, including indets
  h1 <- min(which (y >= stems*0.1)) # compute number of hyperdominants
  h2 <- min(which (y >= stems*0.2)) # compute number of hyperdominants
  h3 <- min(which (y >= stems*0.3)) # compute number of hyperdominants
  h4 <- min(which (y >= stems*0.4)) # compute number of hyperdominants
  h5 <- min(which (y >= stems*0.5)) # compute number of hyperdominants
  h6 <- min(which (y >= stems*0.6)) # compute number of hyperdominants
  h7 <- min(which (y >= stems*0.7)) # compute number of hyperdominants
  h8 <- min(which (y >= stems*0.8)) # compute number of hyperdominants
  h9 <- min(which (y >= stems*0.9)) # compute number of hyperdominants
  s <- length(csum1) # compute number of species in subsample
  perc1 <- round(h1*100/s,2) # compute percentage hyperdominants
  perc2 <- round(h2*100/s,2) # compute percentage hyperdominants
  perc3 <- round(h3*100/s,2) # compute percentage hyperdominants
  perc4 <- round(h4*100/s,2) # compute percentage hyperdominants
  perc5 <- round(h5*100/s,2) # compute percentage hyperdominants
  perc6 <- round(h6*100/s,2) # compute percentage hyperdominants
  perc7 <- round(h7*100/s,2) # compute percentage hyperdominants
  perc8 <- round(h8*100/s,2) # compute percentage hyperdominants
  perc9 <- round(h9*100/s,2) # compute percentage hyperdominants
  return(c(n,s,h1,h2,h3,h4,h5,h6,h7,h8,h9,perc1,perc2,perc3,perc4,perc5,perc6,perc7,perc8,perc9,stems))
}
# define function to rarefy SAD matrix, calculating TS, H#, H%, Alpha, Gamma values and their SDs for each subsample
# input matrix of plot SADs, by default iterations, and sample is done without replacement (choose replacement = T to alter)
Rarefactionperc <- function(df,repeats = 100,replacement=F){
  output <- matrix(ncol=40, nrow=nrow(df)) # initiate empty matrix of correct dimensions to be filled by operations below
  for (j in seq(1,nrow(df),by=1)){
    reps <- replicate(repeats,subsampleperc(df,n=j,replacement = replacement)) # replicate function call to masterloop function above with n = j plots 100 times
    SDs <- rowSDs(reps[2:20,],na.rm = T) # take SD across iterations for total species, H#, and H% replicates
    output[j,]<-c(rowMeans(reps),SDs) # take the mean across all iterations of subsample for each measurement aka take the mean of the rows of iterations
    cat(paste0(round(j / length(df) * 100), '%'))
  }
  output1 <-data.frame(output) # convert to DF
  output2 <- output1[complete.cases(output1),] # remove empty rows
  names(output2) <- c('Plots','TS','h1','h2','h3','h4','h5','h6','h7','h8','h9','perc1','perc2','perc3','perc4','perc5','perc6','perc7','perc8','perc9','stems','TSSD','h1SD','h2SD','h3SD','h4SD','h5SD','h6SD','h7SD','h8SD','h9SD','perc1SD','perc2SD','perc3SD','perc4SD','perc5SD','perc6SD','perc7SD','perc8SD','perc9SD') # give columns meaningful names
  return(output2)
}
# define function to replace data in specified set of columns in one DF with data from corresponding columns with same name in other DF
# input original DF to be changed (here WITHOUT replacement) and imputed DF to replace values from original (here WITH replacement)
replace_imputedperc <- function(original, imputed, toChange = c('TSSD','h1SD','h2SD','h3SD','h4SD','h5SD','h6SD','h7SD','h8SD','h9SD','perc1SD','perc2SD','perc3SD','perc4SD','perc5SD','perc6SD','perc7SD','perc8SD','perc9SD')){
  toChange <- toChange
  for(i in 1:length(toChange)){
    original[toChange[i]] <- imputed[toChange[i]]
  }
  return(original)
}
# define function to rarefy SAD without and with replacement, calculate TS, H#, H%, Alpha, Gamma values and their SDs for each subsample
# and output 
completeRareperc <- function(df,repeats = 100){
  without <- Rarefactionperc(df,repeats=repeats,replacement = F)
  with <- Rarefactionperc(df,repeats=repeats,replacement = T)
  final <- replace_imputedperc(without,with)
  return(final)
}

# RAREFACTION FOR HYPERDOM IDENTITIES (50%) FUNCTION
# function to calculate the number of times each species classifies as hyperdominant across repetions at different levels of subsampling
hyperSubProps <- function(df,sequence = seq(from = 1, to = 10, by =0.5),repetitions = 100){
  dfsum <- sort(colSums(df),decreasing=T) # sum matrix of SADs to vector SAD across dataset and sort by dereasing abundance
  dfsum1 <- dfsum[names(dfsum) != 'Indet'] # remove Indets from fit
  hypersub <- matrix(nrow = length(dfsum1), ncol = length(sequence)) # initiate matrix in which to store results
  # for each level of subsampling - defined by the divisor applied to the dataset number of plots
  for (i in seq_along(sequence)){ # for i in seq_along(sequence) equivalent of enumerate in Python, i gives integer index, sequence[i] gives entry of interest in sequence
    divisor <- sequence[i] # list[i] gives entry of interest, here entry form seq of divisors
    hypres <- matrix(nrow = length(dfsum1), ncol = repetitions) # initiate matrix with rows = species, columns = iterations, entries = yes / no species is hyperdominant in this iteration
    row.names(hypres) <- names(dfsum1) # set names as species names
    # for each iteration of subsampling within each level of subsampling
    for (j in seq(1,repetitions)){
      sample1 <- df[sample(1:nrow(df),nrow(df)/divisor,replace = F),] # sample number of plots from SAD according to level of subsampling
      sample2 <- sample1[,colSums(sample1!=0)>0] # remove 0 entries (may be unnecessary)
      sample3 <- sample2[, names(sample2)!='Indet'] # remove indets from fitting
      csum <- sort(colSums(sample3),decreasing = TRUE) # sum SAD matrix to SAD vector across all sampled plots
      y <- cumsum(csum) # cumulatively sum abundances sorted by decreasing magnitude
      stems <- sum(colSums(sample2)) # compute number of stems in subsample, including indets
      h <- min(which (y >= stems*0.5)) # compute number of hyperdominants
      csumdf <- data.frame(csum) # convert to data.frame (probably unnecessary)
      rnam <- row.names(csumdf)[1:h] # subset species names to just those of the hyperdominants in this sample
      hypres[,j] <- as.numeric(row.names(hypres) %in% rnam) # use whether species is in rnam as boolean indexer and convert to 0s and 1s with as numeric
    }
    hypersub[,i] <- rowSums(hypres) # sum across iterations to get the number of times a species is considered hyperdominant across all iterations at this sampling level
  }
  hypersub1 <- data.frame(hypersub) # convert matrix to dataframe
  row.names(hypersub1) <- names(dfsum1) # set names of dataframe as those of the species
  colnames(hypersub1) <- sequence # set column names as the sampling level / divisors of the data
  hypersub1$rank <- seq(1,dim(hypersub1)[1]) # add column with observed rank of species by abundance
  # hypersub1$avScore <- rowSums(hypersub1[,2:19])/18 # compute average 'score' number of iterations for which species is hyperdominant at each sampling level
  hypersub1$avScore <- rowSums(hypersub1[,2:length(hypersub1)-1])/(length(hypersub1)-1) # compute average 'score' number of iterations for which species is hyperdominant at each sampling level
  hypersub1$Hrank <- rank(-hypersub1$avScore) # create ranking based on average score, -rank from largest to smallest
  return(hypersub1)
}

#### EMPIRICAL AND UNCORRECTED EXTRAPOLATION FUNCTIONS 
## Compilation of functions to fit and extrapolate Log Series from Ter Steege et al, Scientific reports 10.1 (2020): 1-13.
## Code from https://github.com/piLaboratory/atdn-estimates-2019
#' utility function: incomplete beta function
ibeta <- function(x,a,b, log=FALSE){
  y <- pbeta(x,a,b, log.p=TRUE) + lbeta(a,b)
  if(log)
    return(y)
  else
    exp(y)
}
#' CDF of logseries, using incomplete beta function (https://en.wikipedia.org/wiki/Logarithmic_distribution)
pls2 <-function(x, alpha, N){
  p <- N/(N+alpha)
  1 + ibeta(p, x+1, 1e-12)/log(1-p) 
}
#' Continuous approximation for quantile function for Log-series distribution
qls2 <- function(p, N, alpha, lower=3e-9, upper=3e9){ 
  f2 <- function(target){
    f1 <- function(x) pls2(x,alpha, N) - target
    uniroot(f1, lower=lower, upper=upper, extendInt = "yes")$root
  }
  sapply(p, f2)
}
#' Log-series RAD
#' @description Generates a given number of points of the RAD of a
#'     LS, given the total number of species and the parameters of
#'     the distribution.
rad.ls <- function(S, N, alpha, npoints = round(S), ...){
  S.r <- round(S)
  if(missing(alpha))
    alpha <- fishers.alpha(N = N, S = S)
  pp <- rev(ppoints(S.r))
  x <- seq(1,S.r, length=npoints)
  y <- qls2(p = pp[x], N = N, alpha = alpha, ...)
  data.frame(x, y)
}
## Function to compute dominance results 10% - 90% from SAD vector
#' Finds number of hyperdominants (h), total species (ts), percentage H out of total species (per), total individuals (N.tot)
#' @param x vector of population sizes.
#' @param N.tot total number of individuals. If missing calculations are done on 'sum(x)'.
nDomH <- function(x, N.tot){
  if(missing(N.tot))
    N.tot <- sum(x)
  y <- cumsum(sort(x[x>0], decreasing = TRUE))
  s <- length(unlist(x)) # unlist in case SAD is in df as opposed to vector form
  stems <- N.tot # compute number of stems in subsample, including indets
  h1 <- min(which (y >= stems*0.1)) # compute number of hyperdominants
  h2 <- min(which (y >= stems*0.2)) # compute number of hyperdominants
  h3 <- min(which (y >= stems*0.3)) # compute number of hyperdominants
  h4 <- min(which (y >= stems*0.4)) # compute number of hyperdominants
  h5 <- min(which (y >= stems*0.5)) # compute number of hyperdominants
  h6 <- min(which (y >= stems*0.6)) # compute number of hyperdominants
  h7 <- min(which (y >= stems*0.7)) # compute number of hyperdominants
  h8 <- min(which (y >= stems*0.8)) # compute number of hyperdominants
  h9 <- min(which (y >= stems*0.9)) # compute number of hyperdominants
  # s <- length(csum1) # compute number of species in subsample
  perc1 <- round(h1*100/s,2) # compute percentage hyperdominants
  perc2 <- round(h2*100/s,2) # compute percentage hyperdominants
  perc3 <- round(h3*100/s,2) # compute percentage hyperdominants
  perc4 <- round(h4*100/s,2) # compute percentage hyperdominants
  perc5 <- round(h5*100/s,2) # compute percentage hyperdominants
  perc6 <- round(h6*100/s,2) # compute percentage hyperdominants
  perc7 <- round(h7*100/s,2) # compute percentage hyperdominants
  perc8 <- round(h8*100/s,2) # compute percentage hyperdominants
  perc9 <- round(h9*100/s,2) # compute percentage hyperdominants
  return(c(s,h1,h2,h3,h4,h5,h6,h7,h8,h9,perc1,perc2,perc3,perc4,perc5,perc6,perc7,perc8,perc9,stems))
}
## Empirical and extrapolation results master function 
### input matrix of plot SADs and predicted number of regional stems
extrapH <- function(df,stemsRegion){
  df1 <- df[, names(df)!='Indet'] # excluding indets matrix of SADs
  stemsDet <- sum(colSums(df1)) # number of identified stems (not unidentified)
  obStats <- subsamplestats(df) # output observed stats: n,s,h,perc,alphaMLE,alphaSE,gammaMLE,gammaSE,stems
  obsperc <- subsampleperc(df) # output observed hyperdoms and percs at different H thresholds
  Alpha <- obStats[5] # extract observed fitted alpha from obs stats
  AlphaMin <- Alpha - 1.96*obStats[6] # derive alpha CIs via standard deviations of fitted covariance matrix
  AlphaMax <- Alpha + 1.96*obStats[6] # derive alpha CIs via standard deviations of fitted covariance matrix
  percIdent <-stemsDet/obStats[9] # percent of stems identified to species (identified / unidentified)
  speciesPred <- Alpha*log(1+(stemsRegion/Alpha)) # Compute regional expected number of species given fitted alpha, as alpha*ln(1+(N/alpha)) where N is total number of predicted individuals
  speciesPredMin <- AlphaMin*log(1+(stemsRegion/AlphaMin)) # Propogate Alpha CIs to regional species total prediction
  speciesPredMax <- AlphaMax*log(1+(stemsRegion/AlphaMax)) # Propogate Alpha CIs to regional species total prediction
  obFit <- radpred(sad = "ls", coef =  list(N = stemsDet,alpha = Alpha), S = obStats[2], N = stemsDet)
  empStats <- nDom1(obFit$abund,sum(obFit$abund)/percIdent) # compute h,ts,perc for empirical LS fit NB Multiply total individuals in fitted SAD by identification percent to include indets!
  LSerrorH <- empStats$h/obStats[3] # ratio of observed to fitted number of hyperdominants as elementary proxy for error induced by fit to hyperdominance results
  regDist <- rad.ls(S=speciesPred,N=stemsRegion*percIdent,alpha=Alpha) # generate predicted Log Series RAD with S and alpha from empirical fit, N as predicted total stems multiplied by percent identification (to exclude indets from fit)
  regDistMin <- rad.ls(S=speciesPredMin,N=stemsRegion*percIdent,alpha=AlphaMin) # extrap Log Series RAD CIs with S and alpha minima from empirical fit, N as predicted total stems multiplied by percent identification (to exclude indets from fit)
  regDistMax <- rad.ls(S=speciesPredMax,N=stemsRegion*percIdent,alpha=AlphaMax) # extrap Log Series RAD CIs with S and alpha maxima from empirical fit, N as predicted total stems multiplied by percent identification (to exclude indets from fit)
  regStats <- nDom1(regDist$y,sum(regDist$y)/percIdent) # compute h,ts,perc for regional LS fit NB Multiply total individuals in fitted SAD by identification percent to include indets!
  regStatsMin <- nDom1(regDistMin$y,sum(regDistMin$y)/percIdent) # compute h,ts,perc for regional LS Min NB Multiply total individuals in fitted SAD by identification percent to include indets!
  regStatsMax <- nDom1(regDistMax$y,sum(regDistMax$y)/percIdent) # compute h,ts,perc for regional LS Max NB Multiply total individuals in fitted SAD by identification percent to include indets!
  regStatsH <- nDomH(x = regDist$y,N.tot=sum(regDist$y)/percIdent) # compute h,ts,perc for regional LS fit NB Multiply total individuals in fitted SAD by identification percent to include indets!
  regStatsHMin <- nDomH(x=regDistMin$y,N.tot=sum(regDistMin$y)/percIdent) # compute h,ts,perc for regional LS Min NB Multiply total individuals in fitted SAD by identification percent to include indets!
  regStatsHMax <- nDomH(x=regDistMax$y,N.tot=sum(regDistMax$y)/percIdent) # compute h,ts,perc for regional LS Max NB Multiply total individuals in fitted SAD by identification percent to include indets!
  return(c(obStats[3],obStats[2],obStats[4],Alpha,AlphaMin,AlphaMax,obStats[9],dim(df)[1],percIdent,empStats$h,empStats$ts,empStats$per,empStats$N.tot,regStats$h,regStatsMin$h,regStatsMax$h,regStats$ts,regStatsMin$ts,regStatsMax$ts,regStats$per,regStatsMin$per,regStatsMax$per,regStats$N.tot,stemsRegion,LSerrorH,regStats$h/LSerrorH, regStats$per/LSerrorH,regStatsH,regStatsHMin,regStatsHMax))
}
### Function to bias-correct Log-Series extrapolated values
####NOTE FOR THIS FUNCTION HAVE TO REMOVE INDET FROM DATASET AND PREDATASET, DO NOT REMOVE THEM FROM NON-CORRECTION FUNCTIONS
extrapPradoH <- function(region,predataset,dataset,uncex,minSpecies,maxSpecies,repetitions){
  # tocorrect <- uncorrected[uncorrected$Region == region,]
  tocorrectH <- uncex[uncex$Region == region,]
  Tot.t <- tocorrectH$Tot.t
  Tot.A <- tocorrectH$Tot.a
  Samp.A <- sum(predataset[,4])
  N.plots <- dim(predataset)[1]
  dataset$dens.ha <- dataset$N.ind/Samp.A
  dataset$k <- est.kv(mu=dataset$dens.ha, nzeroes=N.plots-dataset$N.plots, Nplots=N.plots)
  lm.k <-lm(log(k)~log(dens.ha), data=dataset, subset=k<1)
  mc.cores <- detectCores() -1 # set number of cores to use to number of cores computer has (8) -1 to leave one spare
  S1 <- round(runif(repetitions, minSpecies, maxSpecies))
  sim.ls.rad <- mclapply(S1, sim.rad, N = Tot.t, sad = "ls",mc.cores = mc.cores, upper = 1e12)
  sim.ls.samp <- mclapply(sim.ls.rad, sim.radsamp, tot.area  = Tot.A,n.plots = N.plots, lmk.fit = lm.k, mc.cores = mc.cores)
  sim.ls.estS <- data.frame(S = S1,
                            S.est.rnd = unlist( mclapply(sim.ls.samp,
                                                         function(x) ls.estS(x$rnd.samp, N = Tot.t),
                                                         mc.cores = mc.cores)),
                            S.est.clump = unlist( mclapply(sim.ls.samp,
                                                           function(x) ls.estS(x$clump.samp, N = Tot.t),
                                                           mc.cores = mc.cores)))
  lm.S<- lm(S ~ S.est.clump, data = sim.ls.estS)
  # predicted central, upper and lower estimates of species richness for central and alpha-fit defined upper and lower bounds
  S.pred <- predict.lm(lm.S,data.frame(S.est.clump=tocorrectH$S.mean),se.fit = T, interval = 'prediction') # predicted central S
  # S.pred.l <- predict.lm(lm.S,data.frame(S.est.clump=tocorrect$S.min),se.fit = T) # predicted lower S
  # S.pred.u <- predict.lm(lm.S,data.frame(S.est.clump=tocorrect$S.max),se.fit = T) # predicted upper S
  hyper.true <- sapply(sim.ls.rad, nDom, p=0.5)
  sim.ls.rad2 <- mclapply(sim.ls.estS$S.est.clump, sim.rad, N = Tot.t, sad = "ls", mc.cores = mc.cores, upper = 1e12)
  hyper.est.clump <- sapply(sim.ls.rad2, nDom, p=0.5)
  sim.ls.estHyp <-
    cbind(sim.ls.estS, H = hyper.true, H.clump = hyper.est.clump) %>%
    mutate(H.prop = H/S, H.clump.prop = H.clump / S.est.clump)
  lm.H<- lm(H ~ H.clump, data = sim.ls.estHyp)
  H.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$H.mean),se.fit = T, interval = 'prediction')
  # H.pred.l <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$H.min),se.fit = T)
  # H.pred.u <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$H.max),se.fit = T)
  h1.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h1),se.fit = T, interval = 'prediction')
  h2.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h2),se.fit = T, interval = 'prediction')
  h3.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h3),se.fit = T, interval = 'prediction')
  h4.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h4),se.fit = T, interval = 'prediction')
  h5.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h5),se.fit = T, interval = 'prediction')
  h6.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h6),se.fit = T, interval = 'prediction')
  h7.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h7),se.fit = T, interval = 'prediction')
  h8.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h8),se.fit = T, interval = 'prediction')
  h9.pred <- predict.lm(lm.H,data.frame(H.clump=tocorrectH$h9),se.fit = T, interval = 'prediction')
  lm.H.prop <- lm(H.prop ~ H.clump.prop, data = sim.ls.estHyp)
  H.prop.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$Prop.mean/100),se.fit = T, interval = 'prediction') # central predicted Hprop
  perc1.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc1/100),se.fit = T, interval = 'prediction')
  perc2.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc2/100),se.fit = T, interval = 'prediction')
  perc3.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc3/100),se.fit = T, interval = 'prediction')
  perc4.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc4/100),se.fit = T, interval = 'prediction')
  perc5.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc5/100),se.fit = T, interval = 'prediction')
  perc6.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc6/100),se.fit = T, interval = 'prediction')
  perc7.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc7/100),se.fit = T, interval = 'prediction')
  perc8.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc8/100),se.fit = T, interval = 'prediction')
  perc9.pred <- predict.lm(lm.H.prop,data.frame(H.clump.prop=tocorrectH$perc9/100),se.fit = T, interval = 'prediction')
  # All results: H, TS, H%, Centre, min, max
  corrected <- c(H.pred$fit[1],H.pred$fit[2],H.pred$fit[3],
                 h1.pred$fit[1],h1.pred$fit[2],h1.pred$fit[3],
                 h2.pred$fit[1],h2.pred$fit[2],h2.pred$fit[3],
                 h3.pred$fit[1],h3.pred$fit[2],h3.pred$fit[3],
                 h4.pred$fit[1],h4.pred$fit[2],h4.pred$fit[3],
                 h5.pred$fit[1],h5.pred$fit[2],h5.pred$fit[3],
                 h6.pred$fit[1],h6.pred$fit[2],h6.pred$fit[3],
                 h7.pred$fit[1],h7.pred$fit[2],h7.pred$fit[3],
                 h8.pred$fit[1],h8.pred$fit[2],h8.pred$fit[3],
                 h9.pred$fit[1],h9.pred$fit[2],h9.pred$fit[3],
                 S.pred$fit[1],S.pred$fit[2],S.pred$fit[3],
                 H.prop.pred$fit[1]*100,H.prop.pred$fit[2]*100,H.prop.pred$fit[3]*100,
                 perc1.pred$fit[1]*100,perc1.pred$fit[2]*100,perc1.pred$fit[3]*100,
                 perc2.pred$fit[1]*100,perc2.pred$fit[2]*100,perc2.pred$fit[3]*100,
                 perc3.pred$fit[1]*100,perc3.pred$fit[2]*100,perc3.pred$fit[3]*100,
                 perc4.pred$fit[1]*100,perc4.pred$fit[2]*100,perc4.pred$fit[3]*100,
                 perc5.pred$fit[1]*100,perc5.pred$fit[2]*100,perc5.pred$fit[3]*100,
                 perc6.pred$fit[1]*100,perc6.pred$fit[2]*100,perc6.pred$fit[3]*100,
                 perc7.pred$fit[1]*100,perc7.pred$fit[2]*100,perc7.pred$fit[3]*100,
                 perc8.pred$fit[1]*100,perc8.pred$fit[2]*100,perc8.pred$fit[3]*100,
                 perc9.pred$fit[1]*100,perc9.pred$fit[2]*100,perc9.pred$fit[3]*100)
  return(corrected)
}

# function to plot Preston fits for one continent
fitPlot <- function(df,type = 'AIC',color='red'){
  df1 <- df[, names(df)!='Indet'] # remove indets from fitting
  csum <- colSums(df1) # sum SAD matrix to SAD vector across all sampled plots
  og.ls <- fitsad(csum,'ls') # Log Series fitting
  # og.pl <- fitsad(x=csum, sad="poilog") # Poisson Log-Normal fitting
  # mod<-fitdistr(csum,"negative binomial") # initial values for nbinom
  # og.nbin <- fitsad(x=csum, sad="nbinom",start.value = c("size"=mod$estimate[1],"mu"=mod$estimate[2]))
  # predict octaves for SAD plotting based on fitted models
  og.ls.oc <- octavpred(og.ls)
  # og.pl.oc <- octavpred(og.pl)
  # og.nbin.oc <- octavpred(og.nbin)
  plot(octav(csum)) # plot octaves of original data
  lines(og.ls.oc, col=color)
  # lines(og.pl.oc, col="red")
  # lines(og.nbin.oc,col="green")
  # legend("topright",c("Log Series Fit"),lty=1, col=c("red")) # , "Poisson Lognormal","Negative Binomial" ,"red", "green"
}
#function to plot RAD fits for one continent
fitPlot1 <- function(df,type = 'AIC',color='red'){
  df1 <- df[, names(df)!='Indet'] # remove indets from fitting
  csum <- colSums(df1) # sum SAD matrix to SAD vector across all sampled plots
  og.ls <- fitsad(csum,'ls') # Log Series fitting
  # og.pl <- fitsad(x=csum, sad="poilog") # Poisson Log-Normal fitting
  # mod<-fitdistr(csum,"negative binomial") # initial values for nbinom
  # og.nbin <- fitsad(x=csum, sad="nbinom",start.value = c("size"=mod$estimate[1],"mu"=mod$estimate[2]))
  # generate predicted RAD values
  og.ls.rad <- radpred(og.ls)
  # og.pl.rad <- radpred(og.pl)
  # og.nbin.rad <- radpred(og.nbin)
  plot(rad(csum))
  lines(og.ls.rad, col=color)
  # lines(og.pl.rad, col="red")
  # lines(og.nbin.rad, col="green")
  # legend("topright",c("Log Series Fit"),lty=1, col=c("red")) # , "Poisson Lognormal","Negative Binomial" ,"red", "green"
}




