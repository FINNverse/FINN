


cohorts_{t+1} = Growth_process(Light, Temp, Prec, DBH) + Mort_process(Light, Temp, Prec, DBH) + Reg_process(Light, Temp, Prec, DBH)


cohorts_[t+1] = Growth_NN(Light, Temp, Prec, DBH) + Mort_process(Light, Temp, Prec, DBH) + Reg_process(Light, Temp, Prec, DBH)



finn(
  data = observed_data,
  env = climate_data,
  disturbance = disturbance_regime,
  init = initial_cohorts,
  regenerationProcess = createProcess(~ Temp^2 + Prec + Temp:Prec, func = regeneration),
  growthProcess = createProcess(~ Temp^2 + Prec + Temp:Prec, func = growth),
  mortalityProcess = createProcess(~ Temp^2 + Prec + Temp:Prec, func = mortality)
  )


env_dt = fread("data/calibration-data/BCI-1h-patch-V2/env_dt.csv")
stand_dt = fread("data/calibration-data/BCI-1h-patch-V2/stand_dt.csv")
obs_df = fread("data/calibration-data/BCI-1h-patch-V2/obs_df.csv")

env_dt[,.(year, siteID = 1,Temp = T_mean, Prec)]


stand_dt[,.(year = census, siteID, species, ba = round(ba, 4), trees, dbh = round(dbh_cm, 4), g = round(g,4), r, m = round(m,3))][order(year, siteID, species)]
