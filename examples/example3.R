rm(list=ls())

library(DRAFT)

##

## Fit example bSIR model - dp and dq are fitted, ts and dL are provided
## -------------------------------------------------------------------------

pop = 4035777 
 

filename = 'lr_ebola.csv'

epi_model = 3

Tg = 12

sigma = NULL

ts = "2014-07-30"
dL = 20  

dp = dq = NULL

results <- runbSEIR(filename = filename, epi_model = epi_model, pop = pop, Tg = Tg, dp = dp, dq = dq, ts = ts, dL = dL)

