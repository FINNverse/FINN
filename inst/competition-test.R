library(FINN)
library(data.table)
library(ggplot2)

# check sensitivity of competition function

# sensitivity to dbh vs parGlobal

dbhTest = seq(0,200, length.out = 200)
parGlobalTest = seq(0.0001, 1, length.out = 9)

out_dt <- data.table()
for(dbh_i in dbhTest){
  for(parGlobal_i in parGlobalTest){
    i_cohort_df =
      data.frame(
        Species = rep(c(0),10),
        dbh = c(rep(dbh_i,5),rep(dbh_i+1,5)),
        nTree = 1
      )
    i_cohort_list = cohort_df2arrays(i_cohort_df)
    i_dbh2height = parGlobal_i
    compOut <- FINN::competition(
      dbh = i_cohort_list$dbh, Species = i_cohort_list$Species, parGlobal = i_dbh2height
    )
    out_dt <- rbind(
      out_dt,
      data.table(dbh = dbh_i,
                 parGlobal = parGlobal_i,
                 AL = as.vector(compOut),
                 heights = c(rep("small",5), rep("big", 5)))
    )

  }
}

ggplot(out_dt, aes(y = AL, x = dbh, color = heights))+
  geom_point()+
  facet_wrap(~factor(parGlobal))
