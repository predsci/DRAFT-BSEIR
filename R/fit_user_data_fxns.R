##
## Driver and Functions for fitting and forecasting User provided incidence data
##



#' Main driver for using \code{DRAFT} to fit user provided incidence data.
#'
#' \code{runDRAFT} accepts a user provided incidence dataframe and uses that along with the user provided population and generation time, Tg, (and if relevant the latent period sigma) to fit the incidence data.  Additionally, when specified by the user, \code{runDRAFT} will also generate a forecast for the incidence.
#' Data cadence is arbitrary but at most can be monthly.
#' We support S-I-R and S-E-I-R models for a single population with a fixed (non-time-dependent) force of infection.
#' @param inc_data Dataframe containing incidence data.  Must contain 'date' and 'cases' columns.  The 'date' column must either be Date-class or convert to Date-class using as.Date(inc_data$date, format="\%Y-\%m-\%d").
#' @param Tg Numeric, generation time in days. Default is 3 days
#' @param pop Integer population of the region for which incidence is provided
#' @param sigma inverse of of the latent period in days. Needed only for an SEIR model. Default NULL
#' @param epi_model - integer 1 (SIR), 2 (SEIR), 3 (SIR with behavior terms). Default is 1-SIR.
#' @param dp Proporiton of susceptible contact rate (0-1) following behavior modification. For example: If susceptibles reduce their contacts by 25\%, set \code{dp=0.75}.  Only used if inc_data requires a forecast.
#' @param dq Proportion of infectious contact rate (0-1) following behavior modification.  If infectious cases reduce their contacts by one half, set \code{dq=0.5}.  Only used if inc_data requires a forecast.
#' @param ts Date class.  Start date of behavior modification.
#' @param dL Number of days for behavior modification to completely take effect.
#' @return A list with the input and entire output of the run.  List entries:
#'   \itemize{
#'     \item \emph{$mydata}: A data structure containing the input data.
#'     \item \emph{$rtn}: The best-fit profile.
#'     \item \emph{$profile}: The full MCMC distribution of incidence profiles as a matrix. 
#'     \item \emph{$tab}: The full MCMC distribution of parameters as a matrix.
#'   }
#'   Additional output is written to a subdirectory 'user_data_*' within the current working directory.  File list:
#'   \itemize{
#'     \item \emph{results-user_data-*.png}: This image shows the fit relative to the data as well as critical-parameter distributions. 
#'     \item \emph{user_data-incidence.png}: Incidence data plot.
#'     \item \emph{mcmc-user_data-*.RData}: A list that includes the resulting MCMC parameter distributions.
#'     \item \emph{profiles-user_data-*.RData}: A list that includes resulting distribution of incidence profiles.
#'   }
#' @details Data fitting is done using a Markov Chain Monte Carlo (MCMC) procedure.  While generally discussed in the context of optimization, this procedure results in a mapping of the objective distribution.  Thus the results of runDRAFT() come in the form of a distribution of parameters and the resulting distribution of incidence profiles.
#' @examples
#' # Run an SEIR model using the incidence file and assuming a population of 1 million people.
#' # The generation time and latent period are set to 2.6 days and 3 days respectively
#' head(incidence_data2)
#' \dontrun{output <- runDRAFT(inc_data=incidence_data2, pop = 1e6, 
#' epi_model = 2, Tg = 2.6, sigma = 3.)}
#'
#' # Run an SIR model using the incidence file and assuming a population of 10,000 people.
#' # The generation time is set to 3 days. (No need to define a latent period.)
#' head(incidence_data1)
#' \dontrun{output <- runDRAFT(inc_data=incidence_data2, pop = 1e5, epi_model = 1, Tg = 3.)}
#'
#' @export
runDRAFT <- function(inc_data=NULL,
           pop = 1e4,
           epi_model = 1,
           Tg = 3,
           sigma = NULL,
           dp = NULL,
           dq = NULL,
           ts = NULL,
           dL = NULL) {
    

    ## Set the MCMC parameters

    # nMCMC = 1e6
    nMCMC = 1e4
    # nlines = 1e
    nlines = 1e3
    plot = 1

    device = 'png'

  # check that inc_data is in proper format
    if (is.data.frame(inc_data)) {
      if (all(c("date", "cases") %in% names(inc_data))) {
        if (class(inc_data$date)=="Date") {
          if (all(diff(inc_data$date)>0)) {
            # inc_data is properly formatted, proceed
          } else {
            # re-order data chronologically
            cat("Warning. The dates of 'inc_data' are not in chronological order.  Re-ordering rows now.\n")
            inc_data = inc_data[order(inc_data$date), ]
          }
        } else {
          inc_data$date = as.Date(inc_data$date)
        }
      } else {
        stop("inc_data must have columns 'date' and 'cases'.\n")
      }
    } else {
      stop("inc_data must be of class data.frame.\n")
    }
    
    
    
	## Since the tanh changes on a time span that is ~ x 2 dL we divide by 2 here 
	
	if(!is.null(dL)) dL = dL/2
	
   # Build the mydata list

   mydata <- build.mydata(inc_data=inc_data, pop = pop, epi_model = epi_model, Tg = Tg, sigma = sigma, dp = dp, dq = dq, ts = ts, dL = dL)
   
   #
   # Build a name for the sub-directory that will hold the results
   #
   myName = Sys.Date()
   myName = gsub(" ", "-", myName)
   if (is.null(sigma)) {
   	subDir = paste0(mydata$model$name, "_", myName, "_Tg_", Tg)
	} else {
	subDir = paste0(mydata$model$name, "_", myName, "_Tg_", Tg,"_sigma_",sigma)	
	}
	
	## If directory exists rename it by appending current time to its name
	if(dir.exists(subDir)) {
		oldDir = subDir
		oldDir = paste0(oldDir,"-",Sys.time())
		oldDir = myName = gsub(" ", "-", oldDir)
		err = file.rename(from = subDir, to = oldDir)
		cat("\n Renaming Previous Directory From:", subDir, " to ", oldDir, '\n\n')
		
	}
	subDir = paste0(subDir,'/')
   mydata$subDir = subDir

	## Create a directory for the run
	
   	dir.create(subDir)

   # Plot the incidence - there is no historic data

   err <- plot_disease(mydata = mydata, device = device)

   ##
   ## Mechanistic Modeling
   ## Pack the information for the run

   par_names <- set.user.data.param.list(epi_model = epi_model)

   opt.list <- set.user.data.opt.list(mydata = mydata)

   run.list <- set.run.list(nMCMC = nMCMC, nlines = nlines, device = device, subDir = mydata$subDir)
   ## Fit the data

   output <- fit.user.data(mydata = mydata, par_names = par_names, opt.list = opt.list, run.list = run.list)
   return(output)

}


#### Read incidence data file and Build the DRAFT data list, mydata
####
#### \code{build.mydata} Reads the user provided csv file with two columns: date and cases
#### The date format is: Year-month-day.
#### cases - integer number of cases
#### cadence can be anything including irregular
#### The code builds and populates the mydata DRAFT list
#### @param inc_data Dataframe containing incidence data.  Must contain 'date' and 'cases' columns.  The 'date' column must either be Date-class or convert to Date-class using \code{as.Date(inc_data$date, format="%Y-%m-%d")}.
#### @param pop - Integer, population of the region for which incidence is provided
#### @param epi_model - integer 1 (SIR) or 2 (SEIR), default is 1
#### @param Tg - Numeric, generation time in days. Default is 3 days
#### @param sigma - inverse of of the latent period. Needed only for an SEIR model. Default 5 days
#### @return mydata - a DRAFT list
#### @examples
#### mydata <- build.mydata(filename = "data.csv", pop = 1e5, epi_model = 2, Tg = 3., sigma = 5.)
#### @export
build.mydata <- function(inc_data=NULL,
                         pop = 1e4,
                         epi_model = 1,
                         Tg = 3,
                         sigma = NULL,
                         dp = NULL,
                         dq = NULL,
                         ts = NULL, 
                         dL = NULL) {
  
  
  
  ## If the user chose behavior modification models check that ts and dL are provided
  
  if (epi_model == 3) {
    if (is.null(ts) ||is.null(dL)) {
      cat("\n For behavior Modification Models User must provide date of \n behavior modification start and time it takes to take full effect: \n ts and dL respectively \n\n Code will STOP \n")
      quit()
    }
  }
  
  
  # user.data = utils::read.csv(file = filename, sep = ",")
  user.data = inc_data
  
  raw = user.data$cases
  
  # dates = as.Date(user.data$date, format = '%m/%d/%y')
  dates = user.data$date
  
  cases = raw
  
  cases[is.na(cases)] <- 0
  
  nperiods = length(dates)
  
  ## Find how many data points we have and how many we may need to forecast
  
  nperiodsData = trimdata.in(longvec = raw)
  
  nperiodsDataUse = nperiodsFit = nperiodsData
  
  ## For behavior modification models dp and dq can be fitted ONLY if this is a fit - and not a forecast
  
  
  mydata = list()
  
  model = list()
  
  # mydata$dates = as.Date(dates, format = '%Y-%m-%d')
  mydata$dates = dates
  
  mydata$years = lubridate::year(dates)
  
  mydata$days  = lubridate::yday(dates)  # day of year
  
  mydata$months = lubridate::month(dates)
  
  mydata$weeks = lubridate::epiweek(dates)
  
  mydata$nperiods = nperiods
  
  mydata$nperiodsData = nperiodsData
  
  mydata$nperiodsDataUse = nperiodsDataUse
  
  mydata$nperiodsFit = nperiodsData
  
  ##
  ## Start Day of behavior change - convert from Date to day number relative to first date of incidence
  
  
  if (epi_model == 3) {
    if (as.numeric(ts - mydata$dates[1]) < 0) {
      cat("Start Day of intervention must be after Day 1 of incdience. \n Code will STOP \n\n")
      quit()
    }
  }
  
  if (epi_model == 3 && nperiodsData < nperiods) {
    if(is.null(dp) || is.null(dq)) {
      cat("\n For a Forecast with Behavior Modification Model user Must provide dp and dq:\n Fractional reduction in mixing rate  of susceptible and infectious populations. \n Code will STOP\n\n")
      quit()
    }
  }
  # build the cadence between the dates
  # the first one is the median of all the other ones
  
  
  ndays = as.numeric(diff(dates, lag = 1))
  
  ndays0 = stats::median(ndays)
  
  mydata$ndays = c(ndays0, ndays)
  
  mydata$season = lubridate::year(dates)[1]
  
  year.end =  lubridate::year(dates)[nperiods]
  
  year.start = lubridate::year(dates)[1]
  
  if (year.end > year.start) {
    mydata$FY = paste0(year.start, '-', year.end)
  } else {
    mydata$FY = year.start
  }
  
  # create an array of zeros - will be used later
  
  zero.vec = rep(0, nperiods)
  
  # determine cadence
  if (all(mydata$ndays == 1)) {
    cadence = "Daily"
  } else if (all(mydata$ndays == 7)) {
    cadence = 'Weekly'
  } else if (all(mydata$ndays >= 28 && mydata$ndays <= 31)) {
    cadence = 'Monthly'
  } else {
    cadence = 'nonuniform'
  }
  
  mydata$cadence = cadence
  
  mydata$disease = 'unknown'
  
  mydata$dataName = 'user_data'
  
  mydata$data_source = 'user'
  
  mydata$data_desc = NULL
  
  mydata$method = 'mech'
  
  if ((nperiods * max(cases)) > 1e6) {
    Temp = 100
  } else if ((nperiodsData * max(cases)) >= 1e5 &&
             (nperiodsData * max(cases)) < 1e6) {
    Temp = 10
  } else {
    Temp = 1
  }
  mydata$Temp = Temp
  
  mydata$epi_model = epi_model
  
  mydata$single = 1
  
  mydata$imodel = 4
  
  mydata$prior = 0
  
  mydata$da = 0
  
  mydata$fit_level = mydata$mod_level = 'unknown'
  
  mydata$Tg = Tg
  
  mydata$sigma = sigma
  
  
  ## Behavior Changing parameters - dp and dq will be used only if this is a forecasting calculation
  ## otherwise they are fitted
  
  mydata$ts = ts
  
  mydata$dL = dL
  
  mydata$dp = dp
  
  mydata$dq = dq
  
  ## End of behavior changine models
  
  ## Build the model list
  
  model$level = 'unknown'
  
  model$name = 'user_data'
  
  model$attr = NULL
  
  model$pop = pop
  
  model$raw_units = '# of cases'
  
  model$factor = 1
  
  model$wght = zero.vec
  
  model$wght[1:nperiodsFit] = 1.0
  
  model$raw   = raw
  
  model$cases = cases
  
  model$epi   = cases
  
  model$gamaepi = lgamma((model$epi + 1))
  
  model$sh = zero.vec
  
  model$temp = zero.vec
  
  model$precip = zero.vec
  
  model$school = zero.vec
  
  model$attr = NULL
  
  mydata$model = model
  
  fit = list()
  
  fit$names = 'user_data'
  
  fit$nregions = 1
  
  mydata$fit = fit
  
  
  return(mydata)
}


##
## General setup of parameter lists routines
##

set.user.data.param.list <- function(epi_model = 1) {
  #### Short Parameter List - Fixed Force of Infection Case
  ####
  ####\code{set.user.data.param.list} creates a list with parameters that the \pkg{DRAFT} code recognizes.
  #### @param epi_model Integer mechanistic model type: SIR (default), SEIR behavior modification SIR (1, 2, and 3)
  #### @return An array with parameter names - the order of parameters will set the order for the min/max arrays also
  #### and the 'mask' for which parameters are optimized (or not)
  #### examples
  #### set.user.data.param.list()

  par_names <-
    c("NH", "Tg", "R0", "sigma",  "pC", "t0", "seed", "e_bckgrnd")

 if (epi_model == 3) {
 	par_names <- c("NH", "Tg", "R0", "sigma",  "pC", "t0", "seed", "e_bckgrnd","dp", "dq", "ts", "dL")
 }
  return(par_names)

}

set.user.data.opt.list <- function(mydata = NULL) {
  #### Create a Logical List  with TRUE/FALSE values for parameter optimization
  ####
  #### \pkg{DRAFT} has a list of model parameters it recognizes, for both uncoupled and coupled runs.
  #### This function sets the values of these parameters to either TRUE or FALSE for the simple
  #### case of user provided data
  #### @param mydata - The \code{DRAFT}  data list
  #### @examples
  #### set.opt.list{mydata = mydata}
  ####

 opt.list = list(NH = FALSE, Tg = FALSE, R0 = TRUE, sigma = FALSE, pC = TRUE, t0 = TRUE, seed = TRUE,
 	e_bckgrnd = TRUE)

 # Behavior modification model
 #
 if (mydata$epi_model == 3) {

 	## None of the behavior modification parameters are fitted in the case of forecast
 	
 	opt.list = list(NH = FALSE, Tg = FALSE, R0 = TRUE, sigma = FALSE, pC = FALSE, t0 = TRUE, seed = TRUE, e_bckgrnd = TRUE, dp = FALSE, 
 		dq = FALSE, ts = FALSE, dL = FALSE)

 	## dp and dq are fitted if NO forecast
 	if (mydata$nperiodsData == mydata$nperiods) {
 		opt.list = list(NH = FALSE, Tg = FALSE, R0 = TRUE, sigma = FALSE, pC = FALSE, t0 = TRUE, seed = TRUE, e_bckgrnd = TRUE, 
 			dp = TRUE, dq = TRUE, ts = FALSE, dL = FALSE)
 	}
 }
 return(opt.list)

}


set.user.data.model <- function(mydata = NULL,
           par_names = NULL,
           opt.list = NULL) {
    #### Setup of Parameters for an MCMC procedure
    ####
    #### \code{set.user.data.model} Creates the arrays for modeling of user
    #### provided incidence data with min/max values, step size,
    #### and initial guess/default values for all the parameters.
    #### @param par_names A character array with parameter names
    #### @param opt.list A Logical list with TRUE or FALSE values for all the parameters
    ####   supported by \pkg{DICE}.
    #### @return A list with min, max, step size and initial guess/default values for the sir parameters
    #### @examples
    #### set.user.data.model(mydata = mydata, par_names = par_names,
    #### run.list = run.list, opt.list = opt.list)


   nparam = length(par_names)

   parmin = parmax = pardx = par = logvec = rep(0, length = nparam)
   names(parmin) = names(parmax) = names(pardx) = names(par) = names(logvec) = par_names
   par_opt = opt.list

   nopt = length(par_opt[par_opt == TRUE])

   NH = mydata$model$pop
   R0 = 1.4
   t0 = 7
   pC = 0.005
   seed = 1
   sigma = mydata$sigma
   Tg = mydata$Tg
	

   if (is.null(Tg))
   	Tg = 3
   if (is.null(sigma))
   	sigma = 5
   ## Population, Tg and sigma

   par["NH"] = parmin["NH"] = parmax["NH"] = NH
   par["Tg"] = parmin["Tg"] = parmax["Tg"] = Tg
   par["sigma"] = parmin["sigma"] = parmax["sigma"] = sigma
   ## R0

   parmin["R0"] = 1.1
   parmax["R0"] = 4.
   par["R0"] = stats::runif(1, 1.2, 2)

   ## pC
   parmin["pC"] = 1e-06
   parmax["pC"] = 1
   par["pC"] = pC * stats::runif(1, 0.8, 1.2)

   nperiodsData = mydata$nperiodsData
   nperiodsData2 = round(nperiodsData/2)
   sum.days = sum(mydata$ndays[1:(nperiodsData)])
   ## Time of first infection
   parmin["t0"] = 1
   parmax["t0"] = sum.days
   
   # If this is a bSIR model need to restrict 't0' differently
   if (!is.null(mydata$ts)) {
   	parmax["t0"] =  mydata$ts - 1
   }
   par["t0"] =round(t0 * stats::runif(1,min(7,parmax['t0']), min(28,parmax['t0'])))

   ## Initial Number of cases

   parmin["seed"] = 1
   parmax["seed"] = 0.001 * mydata$model$pop
   par["seed"] = round(10 * stats::runif(1, 0.8, 1.2))

   est_bckgrnd = mean(mydata$model$cases[1:3], na.rm = TRUE)
   parmin["e_bckgrnd"] = round(est_bckgrnd * 1)
   parmin["e_bckgrnd"] = max(parmin["e_bckgrnd"], 1)
   parmax["e_bckgrnd"] = round(est_bckgrnd * 10)
   parmax["e_bckgrnd"] = max(parmax["e_bckgrnd"], 10)
   par["e_bckgrnd"] = est_bckgrnd

   ##
   dx = 0.01
   pardx[par_names] = dx

   if (mydata$epi_model == 3) {
   	
   	parmin["pC"] = parmax["pC"] = par['pC'] = 1.0
   	
   	parmin["dp"] = parmin["dq"] = 0.001
   	parmax["dp"] = parmax["dq"] = 1.0

   	# ts and dL are never optimized so we just put place holders here 
   	# parmin['ts'] = parmax['ts'] = mydata$ts
   	parmin['dL'] = parmax['dL'] = mydata$dL
   	
   	## Now need to convert ts in absolute date to a ts that is the same number of days away from start of data, tps
   	tps = cumsum(mydata$ndays)
   	del = as.numeric(mydata$ts - mydata$dates[1])
   	ts_for_run = tps[1] + del
   	
   	par['ts'] = ts_for_run
   	parmin['ts'] = parmax['ts'] = ts_for_run
   	
   	
   	par['dL'] = mydata$dL
   	if (mydata$nperiodsData == mydata$nperiods) {
   		par['dp'] = stats::runif(1, 0.5,1.)
   		par['dq'] = stats::runif(1, 0.5,1.)
   	} else {
   		par['dp'] = parmin['dp'] = parmax['dp'] = mydata$dp
   		par['dq'] = parmin['dq'] = parmax['dq'] = mydata$dq
   	}
   }
   
   logbase = 10 #use log base 10 when needed
   logvec[1:nparam] = 1
   logvec["t0"] = 0 #0 #linear for t0



   setup = list(parmin = parmin, parmax = parmax, pardx = pardx, par = par, logbase = logbase, logvec = logvec, nparam = nparam, nopt = nopt, par_names = par_names)

   return(setup)

  }

set.run.list <- function(nMCMC = 1e+05, nlines = NULL, device = "png", subDir = "output") {

  #### Create a List of Parameters For MCMC Procedure
  ####
  #### Pack various parameters related to the MCMC procedure and plotting into a single list.
  ####   These include: The number of MCMC chains, the number of MCMC steps in each chain, the number
  ####   of instances saved for each chain and the device name for plotting the results.
  #### @param nMCMC Integer - the number of steps in each MCMC chain (default is 1e5)
  #### @param nlines Integer - the number of instances saved for the history of the chain (default is to save every 100 steps and not to exceed 1e4)
  #### @param device String  - the device name for plotting.  Default is 'png' but we also support 'pdf'
  #### @param subDir String  - the sub-directory name for all output files of the run. Default it 'output'
  #### @return A list packed with these parameters
  #### @examples
  #### set.run.list{nreal = 1, nMCMC = 1e+05, nlines = 1e3, device = 'png'}
  #### set.run.list{nreal = 3, nMCMC = 5e+06, nlines = 1e4, device = 'pdf'}

  run.list = list()
  run.list$nMCMC = nMCMC

  if (is.null(nlines)) {
    nlines = round(nMCMC/100)
    nlines = max(nlines, 100)
    nlines = min(nlines, 10000)
    if (nlines < nMCMC)
      nlines = nMCMC
  }

  run.list$nlines = nlines
  run.list$device = device
  run.list$subDir = subDir
  ithin <- round(nMCMC/nlines)
  ithin = max(1, ithin)

  run.list$ithin = ithin

  return(run.list)
}
set.imask <- function(par_names = NULL, opt.list = NULL) {

  #### Set a Mask for Parameters
  ####
  #### Given a named logical list of all the parameters that \pkg{DRAFT} supports
  ####   prepare an integer list of the same length with +1/-1 for parameters
  ####   that are/are not optimized
  #### @param par_names - A list with all the \pkg{DRAFT} parameter names
  #### @param opt.list A named logical list with TRUE/FALSE values for each parameter
  #### @return imask An integer list of length nparam with +1 or -1 values
  #### @examples
  #### set.imask{par_names = par_names, opt.list = opt.list}

  nparam = length(par_names)

  imask = rep(-1, nparam)
  names(imask) = par_names

  for (i in 1:nparam) {
    if(opt.list[i] == TRUE) imask[i] = 1
  }

  return(imask)
}


#### Fit User Provided Incidence
####
#### \code{fit.user.data} fits and forecasts user provided incidence using an MCMC procedure
#### and a compartmental S-I-R orr S-E-I-R model
####
#### @param mydata The \pkg{DRAFT} data list
#### @param par_names Array with the \pkg{DRAFT} parameter names
#### @param opt.list Array with TRUE/FALSE for each f the parameters
#### @param run.list A list with MCMC parameters
####
#### @return results A list with all the MCMC fitting/forecasting results and the input mydata list
#### @export
####
#### @examples
#### \dontrun{
#### results <- fit.user.data(mydata = mydata, par_names = par_names,
#### opt.list = opt.list, run.list = run.list)
#### }
#### @export
fit.user.data <- function(mydata = NULL,
           par_names = NULL,
           opt.list = NULL,
           run.list = NULL) {


	tps = cumsum(mydata$ndays)

	device = run.list$device
	subDir = run.list$subDir

	nperiods = mydata$nperiods
	nperiodsFit = mydata$nperiodsFit

	Temp = mydata$Temp

	# seed for RNG

	iseed = set.iseed() %% .Machine$integer.max

	# Here we set the random number generation seed for R
	set.seed(iseed)

	setup <- set.user.data.model(par_names = par_names, mydata = mydata, opt.list = opt.list)
  
	par_names = setup$par_names
	pmin = setup$parmin
	pmax = setup$parmax

	dx = setup$pardx
	par = setup$par
  
	nparam = setup$nparam
	nopt = setup$nopt
	logbase = setup$logbase
	logvec = setup$logvec

	tab = setup$tab

	nreal = run.list$nreal
	ithin = run.list$ithin
	nlines = run.list$nlines
	nMCMC = run.list$nMCMC

	imask <- set.imask(par_names = par_names, opt.list = opt.list)

	pois = 0

	cases = mydata$model$epi
	gamaepi = mydata$model$gamaepi
	wght = mydata$model$wght

	nRnd = 1000
	profile = array(0, c(nRnd, nperiods))

   nlines = run.list$nlines

   tab <- matrix(data = 0, nrow = nlines, ncol = (nparam + 1))

	## Need nperiods + 2 because we patch the data at the beginning and end

	ndays = tps[nperiods] + mydata$ndays[1] + mydata$ndays[nperiods]

	## MCMC Fitting procedure
	##
	if (mydata$epi_model == 1) {
		cat("\n\n Fitting S-I-R model to User Data \n\n")
		out <- .Fortran("fitsir", epi = as.double(cases), gamaepi = as.double(gamaepi), wght = as.double(wght), nparam = as.integer(nparam), par = as.double(par), parmin = as.double(pmin),
			parmax = as.double(pmax), step = as.double(dx), ilog = as.integer(logvec),
			Temp = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(tps), rtn = as.double(rep(0, nperiods)), pois = as.double(pois), ndays = as.integer(ndays), nRrnd = as.integer(nRnd),
			tab = as.single(tab), profile = as.single(profile), PACKAGE="DRAFT")

	} else if (mydata$epi_model == 2) {
		cat("\n\n Fitting S-E-I-R model to User Data \n\n")
		out <- .Fortran("fitseir", epi = as.double(cases), gamaepi = as.double(gamaepi), wght = as.double(wght), nparam = as.integer(nparam), par = as.double(par), parmin = as.double(pmin),
			parmax = as.double(pmax), step = as.double(dx), ilog = as.integer(logvec),
			Temp = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(tps), rtn = as.double(rep(0, nperiods)), pois = as.double(pois), ndays = as.integer(ndays), nRrnd = as.integer(nRnd),
			tab = as.single(tab), profile = as.single(profile), PACKAGE="DRAFT")

	}  else if (mydata$epi_model == 3) {
		cat("\n\n Fitting Behavior Modification S-I-R model to User Data \n\n")
		out <- .Fortran("fitbsir", epi = as.double(cases), gamaepi = as.double(gamaepi), wght = as.double(wght), nparam = as.integer(nparam), par = as.double(par), parmin = as.double(pmin),
			parmax = as.double(pmax), step = as.double(dx), ilog = as.integer(logvec),
			Temp = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(tps), rtn = as.double(rep(0, nperiods)), pois = as.double(pois), ndays = as.integer(ndays), nRrnd = as.integer(nRnd),
			tab = as.single(tab), profile = as.single(profile), PACKAGE="DRAFT")

	} else {
		cat("\n\n Fitting S-I-R model to User Data \n\n")
		out <- .Fortran("fitsir", epi = as.double(cases), gamaepi = as.double(gamaepi), wght = as.double(wght), nparam = as.integer(nparam), par = as.double(par), parmin = as.double(pmin),
			parmax = as.double(pmax), step = as.double(dx), ilog = as.integer(logvec),
			Temp = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(tps), rtn = as.double(rep(0, nperiods)), pois = as.double(pois), ndays = as.integer(ndays), nRrnd = as.integer(nRnd),
			tab = as.single(tab), profile = as.single(profile), PACKAGE="DRAFT")
	}


	# Random fits
	profile = array(out$profile, c(nRnd, nperiods))
	# Best Fit
	rtn = out$rtn

	# History of MCMC 
	tab = matrix(out$tab, ncol = (nparam + 1))

	colnames(tab) = c(par_names, "AICc")
	
	# plot the results
	
	plotlist <- plot_results(rtn = rtn, profile = profile, tab = tab, mydata = mydata, device = device)

	## Dump all the profiles we have to a file
	##
	
	err <- write.profiles(mydata = mydata, rtn = rtn, profile = profile)

	## Dump the MCMC parameters
	
	err <- write.mcmc(mydata = mydata, tab = tab, opt.list = opt.list, run.list = run.list, imask = imask)
	
	## Build and return the results list
	##
	results = list(mydata = mydata, rtn = rtn, profile = profile, tab = tab)

	return(results)

  }
