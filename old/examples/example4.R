rm(list=ls())

library(DRAFT)

##



## Forecast example bSIR model -dp, dq ts and dL are all  provided
## ----------------------------------------------------------------

pop = 4035777 
 
filename = 'lr_ebola_forecast.csv'

epi_model = 3

Tg = 12

sigma = NULL

ts = "2014-07-30"
dL = 20  

dp = 0.002

dq =  0.45

results <- runbSEIR(filename = filename, epi_model = epi_model, pop = pop, Tg = Tg, dp = dp, dq = dq, ts = ts, dL = dL)

