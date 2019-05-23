rm(list=ls())

library(DRAFT)

##



## Forecast example bSIR model -dp, dq ts and dL are all  provided
## ----------------------------------------------------------------

pop = 4035777
 

filename = 'lr_ebola.csv'

epi_model = 1

Tg = 12

sigma = NULL

ts = NULL
dL = NULL 

dp = NULL

dq =  NULL

results <- runbSEIR(filename = filename, epi_model = epi_model, pop = pop, Tg = Tg, dp = dp, dq = dq, ts = ts, dL = dL)

