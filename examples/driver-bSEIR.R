rm(list=ls())

library(bSEIR)

##


## Fit example - SD
## -------------------
epi_model = 1 # or 2, 3

pop = 3e6

## Generation time in days
##
Tg = 2.6

## latent period - needed only for SEIR model, epi_model = 2
##
sigma = NULL

## Parameters for behavior modification model only

dp = dq = ts = dL = NULL

filename = 'data.csv'


results <- runbSEIR(filename = filename, epi_model = epi_model, pop = pop, Tg = Tg, dp = dp, dq = dq, ts = ts, dL = dL)


## Forecast example SD
## ---------------------

filename = "sandiego13forecast.csv"

results <- runbSEIR(filename = filename, epi_model = epi_model, pop = pop, Tg = Tg, dp = dp, dq = dq, ts = ts, dL = dL)




## Fit example bSIR model - dp and dq are fitted, ts and dL are provided
## -------------------------------------------------------------------------
ts = "2017-09-08"

filename = 'lr_ebola.csv'

epi_model = 3

Tg = 12

ts = "2014-07-30"
dL = 10  

dp = dq = NULL

results <- runbSEIR(filename = filename, epi_model = epi_model, pop = pop, Tg = Tg, dp = dp, dq = dq, ts = ts, dL = dL)

