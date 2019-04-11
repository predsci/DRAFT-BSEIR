rm(list=ls())

library(DRAFT)

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

