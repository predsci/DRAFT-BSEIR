set.iseed <- function() {

  #### Set the seed for the RNG
  #### 
  #### Create the seed in the R code for the Random Number Generator
  #### @param none required
  #### @return An integer - seed for the RNG
  #### @examples
  #### set.iseed()

  set.seed(seed = NULL)
  iseed <- stats::runif(1) * 10000
  iseed <- round(iseed)

  return(iseed)
}


makeDir <- function(subDir = NULL) {
  #### Create a sub-directory
  ####
  #### Create a sub-directory for a \pkg{DRAFT} run under the current working directory.   If this sub-directory already exists it does NOT overwrite it.
  #### @param subDir String - the sub-directory name for all output files of the run. Default it 'output'
  #### @examples
  #### makeDir{subDir='output'}

  if (is.null(subDir))
    subDir = "output"
  mainDir = getwd()
  if (!file.exists(subDir)) {
    dir.create(file.path(mainDir, subDir))
  }

  return(err = 0)
}

gammaEpi <- function(epi = NULL) {

  #### Calculate the \eqn{Gamma} Function of Disease Number of Cases
  ####
  #### Given a matrix or vector of integers with monthly/weekly/daily number of disease cases calculate for each spatial region, the \eqn{\Gamma} function of this matrix/vector.
  #### This is done to save computational time during the MCMC procedure: the \eqn{\Gamma} function of the number of cases is needed when evaluating the maximum likelihood.
  #### @param epi A nperiods x nregions matrix of integers  with monthly/weekly/daily number of disease cases for each region
  #### @return gamaepi A matrix  of nperiods x nregions with the \eqn{\Gamma} function of monthly/weekly/daily number of disease cases in each region
  #### @examples
  #### gammaEpi{epi = mydata$model$epi$}
  #### gammaEpi{epi = mydata$fit$epi$  }

  copyepi = epi
  gamaepi = epi

  n = dim(epi)[2]
  nperiods = dim(epi)[1]

  # to save computational time we calculate the gamma(epi+1)
  for (i in 1:nperiods) {
    for (j in 1:n) {
      gamaepi[i, j] <- lgamma((epi[i, j] + 1))
    }
  }
  return(gamaepi)
}

trimdata.in <- function(longvec) {
#### Trim an array of raw cases to find the last observed data point
####
#### @param longvec array with raw cases numbers
####
#### @return nperiodsData the number of last observed data point
####
#### @examples
#### nperiodsData <- trimdata.in(longvec = cases)
####

  noobs = length(longvec)
  last  = noobs

  while (is.na(longvec[last]))
    last = last - 1

  return(last)
}


#### Plot The Time Series of the User Provided Incidence
####
#### \code{plot_disease} Uses plotly to plot the user provided disease time series in both png and html formats.
#### The latter is interactive
#### @param mydata A list with all the available data for this \pkg{DRAFT} run
#### @param device - device type for plots 'png' (default) or 'pdf'
#### plot_disease(mydata = mydata,  device = 'pdf')
#### @return err=0 if plot was created
####
#### @export
plot_disease <- function(mydata = NULL,  device = 'png') {
  
  if (is.null(device))
    device = "png"
  if (is.null(mydata))
    return
  
  subDir = mydata$subDir
  
  if (!dir.exists(subDir)) {
    dir.create(subDir)
  }
  
  nperiods = mydata$nperiods
  
  dates = mydata$dates
  ntps  = length(dates)
  x.axis.label = rep(",ntps")
  ind = seq(from = 1, to = ntps, by = 1)
  x.axis.label[ind] = format(dates[ind], "%b-%d")
  x.axis.label[is.na(x.axis.label)] <- ""
  
  myName = mydata$dataName
  myName = gsub(" ","",myName)
  
  if (tolower(device) == "pdf") {
    filename = paste0(subDir,myName,"-incidence.pdf")
  } else  {
    filename = paste0(subDir,myName,"-incidence.png")
  }
  
  cat("\n\n For a Plot of the Results See: ", filename, "\n\n")
  FY = mydata$FY
  
  ylab = " # Cases"
  xlab = " "
  title =  paste0("  User Provided Incidence: ", FY)
  
  pl = ggplot2::ggplot(data = NULL) + ggplot2::scale_x_continuous(name = "Date", limits = c(1, ntps), breaks = 1:ntps, labels = x.axis.label) + 
    ggplot2::scale_y_continuous(name = ylab) + ggplot2::theme(text = ggplot2::element_text(size = 8, color = "gray20", face = "italic"), 
                                                              axis.text.x = ggplot2::element_text(face = "plain", angle = 90, size = 6), axis.text.y = ggplot2::element_text(face = "plain"))
  pl = pl + ggplot2::geom_line(ggplot2::aes(x = 1:ntps, y = mydata$model$raw), col = "darkred", size = 2)	
  
  pl = pl + ggplot2::annotate("text", x = -Inf, y = Inf, label = title, 
                              hjust = 0, vjust = 2.5, col = "black", family = "serif", size = 3.0)
  
  ggplot2::ggsave(filename = filename, plot = pl, device = device, width = 14, height = 9, units = "cm")
  
  err = 0
  return(err)
  
}


#### Plot the results of \pkg{runDRAFT} - Single Region
####
#### Plot the results of \pkg{runDRAFT} for a single region/patch. We show the incidence along with our fits and
#### if appropriate predictions. We show the best result and randomly selected results from the MCMC chain. This is the ggplot2 version.
#### @param rtn A 1D numeric array with the best direct prediction to the region
#### @param profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the data
#### @param tab MCMC history for the parameters
#### @param mydata A dataframe with all the data available for this \code{runDRAFT} run
#### @param device name for plotting: 'png' or 'pdf
#### @return  Returns \eqn{err = 0} if successful
#### @examples
#### \dontrun{
#### plotFitOnePatch(rtn = rtn, profile = profile,
#### mydata = mydata, ireal = ireal, idevice = device)
#### }
#### @export
plot_results <- function(rtn = NULL, profile = NULL, tab = NULL, mydata = NULL, device = 'png') {

  if (is.null(device))
    device = "png"
  if (is.null(profile))
    return
  
  FY = mydata$FY
  
  ndays = mydata$ndays
  dates = mydata$dates
  nperiods = mydata$nperiods
  nperiodsFit = mydata$nperiodsFit
  nperiodsData = mydata$nperiodsData
  reg.model.name = mydata$model$name
  
  obs  = mydata$model$raw
  
  tps = cumsum(ndays)
  
  nRnd = dim(profile)[1]
  
  ## check to see if output directory exists and if not create it
  subDir = mydata$subDir
  
  err = makeDir(subDir = subDir)
  myName = mydata$dataName
  myName = gsub(" ","",myName)
  if (tolower(device) == "pdf") {
    filename = paste0(subDir, "results-", myName, "-", nperiodsFit, ".pdf")
  } else {
    filename = paste0(subDir, "results-", myName, "-", nperiodsFit, ".png")
  } 
  
  ## for plotting - replace zero's with NA
  if (nperiodsFit < nperiods) {
    index <- which(obs[(nperiodsFit+1):nperiods] == 0)
    if (length(index) >= 1)
      obs[index] = NA
  }
  
  
  ## ylab label
  ylab = "# Cases"
  
  
  # Plot the data and fits/predictions
  plotlist = list() # For saving all the plots
  
  model_mean = rep(0, nperiods)
  for (iweek in 1:nperiods) model_mean[iweek] = mean(profile[, iweek])
  ymax = max(rtn[1:nperiodsData], profile[, 1:nperiodsData], obs[1:nperiodsData], na.rm = TRUE)
  
  breaks = seq(from = 1, to = nperiods, by = 4)
  labels = dates[breaks]
  
  plotlist[[1]] = ggplot2::ggplot(data = NULL) + ggplot2::scale_x_continuous(name = "Date", limits = range(tps), breaks = tps[breaks], labels = labels) + ggplot2::theme(text = ggplot2::element_text(size = 10, color = "gray20", face = "italic"), axis.text.x = ggplot2::element_text(face = "plain", angle = 90), axis.text.y = ggplot2::element_text(face = "plain")) #limits = c(1, nperiods)
  
  step = 1 #max(1, nRnd/10)
  irnd.set = seq(from = 1, to = nRnd, by = step)
  dat.rnd = t(profile[irnd.set, 1:nperiodsFit])
  dat.rnd.pred = t(profile[irnd.set, nperiodsFit:nperiods])
  data.rnd = reshape::melt(dat.rnd)
  data.rnd.pred = reshape::melt(dat.rnd.pred)
  
  plotlist[[1]] = plotlist[[1]] + ggplot2::geom_line(ggplot2::aes(x =rep(tps[1:nperiodsFit], nRnd), y = data.rnd[, 3], group = data.rnd[, 2]), col = "#E495A5", size = 2, alpha = 0.4) + 
    ggplot2::geom_line(ggplot2::aes(x = tps[1:nperiodsFit], y = rtn[1:nperiodsFit]), col = "#39BEB1", size = 1) + 
    ggplot2::geom_line(ggplot2::aes(x = tps[1:nperiodsFit], y = model_mean[1:nperiodsFit]), col = "#099DD7", size = 0.8) + 
    ggplot2::geom_line(ggplot2::aes(x = tps[1:nperiodsFit], y = obs[1:nperiodsFit]), col = "black", na.rm = TRUE) + 
    ggplot2::geom_point(ggplot2::aes(x = tps[1:nperiodsFit], y = obs[1:nperiodsFit]), col = "#24576D", size = 1, na.rm = TRUE) #x = data.rnd[, 1]
  
  if (nperiodsFit < nperiods) {
    plotlist[[1]] = plotlist[[1]] + ggplot2::geom_line(ggplot2::aes(x = (rep(tps[nperiodsFit:nperiods], nRnd) + tps[nperiodsFit - 1]), y = data.rnd.pred[, 3], group = data.rnd.pred[, 2]), col = "#E495A5", size = 2, linetype = 2, alpha = 0.4) + 
      ggplot2::geom_line(ggplot2::aes(x = tps[nperiodsFit:nperiods], y = model_mean[nperiodsFit:nperiods]),	col = "#099DD7", size = 0.8, linetype = 2) + 
      ggplot2::geom_line(ggplot2::aes(x = tps[nperiodsFit:nperiods], y = rtn[nperiodsFit:nperiods]), col = "#39BEB1", size = 1, linetype = 2) +  
      ggplot2::geom_line(ggplot2::aes(x = tps[nperiodsFit:nperiods], y = obs[nperiodsFit:nperiods]), col = "black", size = 1, linetype = 2) + 
      ggplot2::geom_rect(ggplot2::aes(xmin = tps[nperiodsFit], xmax = tps[nperiods], ymin = 0, ymax = ymax), fill = "#D497D3", alpha = 0.4)
  }
  
  reg.name = paste0("   ", c("User Data", "Model-Best", "Model-Mean", "Model-Random"))
  plotlist[[1]] = plotlist[[1]] + ggplot2::annotate("text", x = rep(-Inf, 5), y = rep(Inf, 5), label = c(paste0("   ", mydata$FY), reg.name), hjust = rep(0, 5), vjust = seq(from = 2.5, to = 8.5, by = 1.5), col = c("black", "black", "#39BEB1", "#099DD7", "#E495A5"), family = "serif", size = 3.5)
  
  
  
  nlines = dim(tab)[1]
  nlines2 = nlines/2
  
  
  ## Adding the p/q average line in case of bSIR model fitting only 
  
  if (mydata$epi_model == 3) {
    # calculate mean value for dp and dq 
    
    dL = mydata$dL
    
    del_p = del_q = array(data = 0, c(length(nlines2:nlines), nperiods))
    for (k in nlines2:nlines) {
      dp = tab[k, "dp"]
      dq = tab[k, "dq"]
      ## Need also mean value for t0
      t0 = tab[k, "t0"]
      
      del = as.numeric(mydata$ts - mydata$dates[1])
      
      ts_for_run = tps[1] + del
      
      delta = 0.5 * (1 + tanh((tps - ts_for_run)/dL))
      
      del_p[(k-nlines2+1), ] = 1 - (1 - dp) * delta
      del_q[(k-nlines2+1), ] = 1 - (1 - dq) * delta
    }
    p.mean = q.mean = rep(0, nperiods)
    for (i in 1:nperiods) {
      p.mean[i] = mean(del_p[,i])
      q.mean[i] = mean(del_q[,i])
    }
    
    plotlist[[1]] = plotlist[[1]] + ggplot2::geom_line(ggplot2::aes( x = tps, y = p.mean * ymax), col = 'darkred', size = 2, alpha = 0.7)
    plotlist[[1]] = plotlist[[1]] + ggplot2::geom_line(ggplot2::aes( x = tps, y = q.mean * ymax), col = 'orange' , size = 2, alpha = 0.7)
    plotlist[[1]] = plotlist[[1]] + ggplot2::scale_y_continuous(name = ylab, limits = c(0, ymax), sec.axis = ggplot2::sec_axis(~. / ymax, name = "Fractional Contact Rate"))
    
  } else {
    
    plotlist[[1]] = plotlist[[1]] +	ggplot2::scale_y_continuous(name = ylab, limits = c(0, ymax)) 
  }
  
  
  ## Now for a histogram plot of R0 and pC 
  
  if (mydata$epi_model == 1 || mydata$epi_model == 2) {
    hist.params = c("R0", 'pC')
  } else {
    if(mydata$nperiodsData == mydata$nperiods) {
      hist.params = c("R0", "dp", "dq")
    } else {
      hist.params = "R0"
    }
  }
  
  icount = 1
  
  for (par in hist.params) {
    
    icount = icount + 1
    my.col = which(colnames(tab) == par)
    my.mcmc = tab[nlines2:nlines, my.col]
    my.mean = mean(my.mcmc)
    my.sd = stats::sd(my.mcmc)
    digits.mean = 3
    digits.sd = 4
    if (my.mean > 1) {
      digits.mean = 2
      digits.sd = 3
    }
    my.mean = round(my.mean, digits = digits.mean)
    my.sd = round(my.sd, digits = digits.sd)
    
    my.par = par
    if (my.par == 'dp') my.par = 'p'
    if (my.par == 'dq') my.par = 'q'
    data = data.frame(x = my.mcmc)
    my.color = 'cornflowerblue'
    if (my.par == 'p') my.color = 'brown'
    if (my.par == 'q') my.color = 'orange'
    
    # R CMD check --as-cran workaround for ggplot variables '..density..' and 'x'
    x <- ..density.. <- NULL
    # end workaround 
    
    plotlist[[icount]] = ggplot2::ggplot(data = data, ggplot2::aes(x = x)) + 
      ggplot2::geom_histogram(data = data, ggplot2::aes(y = ..density..), fill = my.color, col = "white", alpha = 0.7) + ggplot2::scale_x_continuous(name = my.par) + 
      ggplot2::scale_y_continuous(name = "") + ggplot2::theme(text = ggplot2::element_text(size = 10, color = "gray20", face = "italic"), axis.text.x = ggplot2::element_text(face = "plain"), axis.text.y = ggplot2::element_text(face = "plain"))
    plotlist[[icount]] = plotlist[[icount]] + ggplot2::annotate("text", x = rep(-Inf, 2), y = rep(Inf, 2), label = c(paste0("   ", "Mean ", my.mean), paste0("   ", "SD ", my.sd)), hjust = rep(0, 2), vjust = seq(from = 2.5, to = 4, by = 1.5), col = c("black", "black"), family = "serif", size = 3.5) + 
      ggplot2::annotate("text", x = rep(Inf, 2), y = rep(Inf, 2), label = my.par, hjust = 2, vjust = 2.5, col = "black", family = "serif", size = 4)
    
  }
  
  if (icount == 3) {	
    mp = gridExtra::grid.arrange(grobs = plotlist, widths = c(1, 1), layout_matrix = rbind(c(1, 1), c(2, 3)))
  } else if (icount == 4) {
    mp = gridExtra::grid.arrange(grobs = plotlist, widths = c(1, 1, 1), layout_matrix = rbind(c(1, 1, 1), c(2, 3, 4)))	
  } else if (icount == 2) {
    mp = gridExtra::grid.arrange(grobs = plotlist, widths = c(1, 1), layout_matrix = rbind(c(1, 1), c(2, NA)))	
  } else {
    mp = plotlist[[1]]
  }
  
  cat("\n\n For a Plot of the Results See: ", filename, "\n\n")	
  ggplot2::ggsave(file = filename, plot = mp, device = device, width = 9, height = 7, units = "in")
  
  return(err = 0)
  
}

write.profiles <- function(mydata=NULL, rtn = NULL, profile = NULL) {
  ####
  #### Save the profiles for Model region
  #### Save the profiles for Model region so that incidence plot can be re-created
  #### 
  #### @param mydata A dataframe with all the data available for this \pkg{DRAFT} run	
  #### @param profile A 2D numeric array holding random predictions the incidence
  #### @param rtn A 1D numeric array with the best direct prediction to the model region
  #### @return  Returns \eqn{err = 0} if successful
  #### @export
  #### @examples
  #### write.profiles(mydata = mydata, rtn = rtn, profile = profile)
  
  
  myName = mydata$dataName
  myName = gsub(" ","",myName)
  
  nperiods = mydata$nperiods	
  nperiodsFit = mydata$nperiodsFit
  
  nRnd = dim(profile)[1]
  
  mean_profile = rep(0, nperiods)
  for(i in 1:nperiods) {
    mean_profile[i] = mean(profile[,i])
  }
  
  # Dump all the profiles we have to a file
  dump = list()
  dump = list(mydata = mydata, rtn = rtn, profile = profile, mean_profile = mean_profile)
  
  filename = paste(mydata$subDir, "/profiles-", myName, "-", nperiodsFit,".RData", sep = "")
  save(dump, file = filename) 
  
  cat("\nSaving Profiles to a Binary File : ", filename, "\n")    
  
  return(err = 0)
}

write.mcmc <- function(tab = NULL, opt.list = NULL, run.list = NULL, mydata = NULL, imask = NULL) {
  
  #### Write the MCMC History of a \pkg{DRAFT} run
  ####
  #### Writes an RData file with the MCMC history for an uncoupled \pkg{DRAFT} run.
  ####   The function also calculates and prints to the screen the statistics for all the
  ####   parameters that were optimized and does a Gaussian fit to these parameters.
  ####   The results of these fits are written to separate a csv file.
  #### @param tab The MCMC history of the direct fit of the model data
  #### @param opt.list A logical list of all the parameters \pkg{DRAFT} recognizes and
  ####   can optimize with TRUE/FALSE
  #### @param run.list a list with parameters used for the MCMC procedure
  #### @param mydata A dataframe with all the data available for this run
  #### @param imask An array of integers with +1/-1 values for parameters that are optimized (or not)
  #### @return err   Returns \eqn{err = 0}
  #### @export
  #### @examples
  #### write.mcmc(tab = tab, opt.list = opt.list, run.list = run.list, mydata = mydata, imask = imas)
  
  myName = mydata$dataName
  myName = gsub(" ","",myName)
  nperiods = mydata$nperiods
  nperiodsFit = mydata$nperiodsFit
  
  ithin = run.list$ithin
  nlines = run.list$nlines
  nMCMC = run.list$nMCMC
  
  myColName = names(opt.list)
  nopt = length(which(imask == 1))
  names.opt = myColName[which(imask == 1)]
  
  # convert MCMC output back to matrix and to an MCMC object
  
  if (!is.null(tab)) {
    
    nparam <- dim(tab)[2] - 1
    nparam1 = nparam + 1
    colnames(tab) = c(myColName, "AICc")
    tab[, "AICc"] <- 2 * tab[, "AICc"] + 2 * nopt
    tab[, "AICc"] <- tab[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)
    
  }
  
  
  # how many steps to burn - set to 1/5 here
  iburn <- nlines/5
  
  # This mcmc object is only from 'iburn' and up and used only for the purpose of the statistics
  results.list = list()
  
  if (!is.null(tab))
    results.list = coda::mcmc(data = tab, start = (ithin * iburn), end = nMCMC, thin = ithin)
  
  # check to see if 'mydata' sub-directory exists, if not create it
  subDir = run.list$subDir
  if (is.null(subDir))
    subDir = "output"
  err = makeDir(subDir = subDir)
  
  
  # print the chains statistics to the screen-only for optimized variables and AICc
  
  if (!is.null(tab)) {
    tmp = results.list
    colnames(tmp) = c(myColName, "AICc")
    print(summary(tmp[,c(names.opt,"AICc")], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE), digits = 2))
  }
  
  
  # save the complete chain here
  filename <- paste(subDir, "/mcmc-", myName, "-", nperiodsFit, ".RData", sep = "")
  cat("\n Writing R object Data file for this Chain: ", filename, "\n")
  
  dump = list()
  dump = list(mydata = mydata, tab = tab, opt.list = opt.list, run.list = run.list, imask = imask)
  save(dump, file = filename)
  
  
  success = 0
  return(success)
  
}
