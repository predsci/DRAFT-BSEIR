rm(list = ls())

library("DICE")

mod_name = c(NAME_2="LR")

data = get.DICE.data(data_source=NULL, mod_level=2, mod_name=mod_name, fit_level=3, year=2014, db_opts=list(DICE_db="predsci", CDC_server=FALSE), disease="ebola")

mydata = data$mydata

dates = mydata$dates

cases = mydata$model$raw

data.to.save = data.frame(dates = dates, cases = cases)

write.csv(data.to.save, 'lr_ebola_forecast.csv')

write.csv(data.to.save, 'lr_ebola.csv')
