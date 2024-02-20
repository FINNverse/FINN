################################################################################
#' This scripts tests the competition function and contains functions that
#' help to identify plausible parameter ranges
################################################################################

library(FINN)
library(data.table)
library(ggplot2)

# check sensitivity of competition function

# sensitivity to dbh vs parGlobal
dbhTest = seq(0,200, length.out = 200)
parGlobal_i = c(0.5,0.6,0.7)
Ncohorts_per_Sp = 4
Nspecies = 3
out_dt <- data.table()
for(dbh_i in dbhTest){
  i_cohort_df =
    data.frame(
      Species = rep(0:(Nspecies-1),each = Ncohorts_per_Sp),
      parGlobal = rep(i_dbh2height, each = Ncohorts_per_Sp),
      heights = rep(c(rep("small",Ncohorts_per_Sp/2), rep("big (dbh + 20)", Ncohorts_per_Sp/2)), Nspecies),
      dbh = rep(c(rep(dbh_i,Ncohorts_per_Sp/2),rep(dbh_i+20,Ncohorts_per_Sp/2)),Nspecies),
      nTree = 1
    )
  i_cohort_list = cohort_df2arrays(i_cohort_df)
  compOut <- FINN::competition(
    dbh = i_cohort_list$dbh, Species = i_cohort_list$Species, parGlobal = i_dbh2height
  )
  out_dt <- rbind(
    out_dt,
    data.table(
      constelation = constelation,
      dbh = dbh_i,
      parGlobal = rep(i_dbh2height, each = Ncohorts_per_Sp),
      AL = as.vector(compOut),
      heights = i_cohort_df$heights
      )
  )
}

ggplot(out_dt, aes(y = AL, x = dbh, color = heights))+
  geom_point()+
  facet_wrap(~factor(parGlobal), ncol = 1)

# potential test cases can be

# do species with lower parGlobal values have less light available than species with higher?
sum(out_dt[parGlobal == parGlobal_i[1]]$AL) < sum(out_dt[parGlobal == parGlobal_i[2]]$AL) &
sum(out_dt[parGlobal == parGlobal_i[2]]$AL) < sum(out_dt[parGlobal == parGlobal_i[3]]$AL)

# Does a tree with a dbh of 500 has a lower height than 150 meter?
# The tallest tree on earth was 116 m with a dbh of 450 cm.
# FIXME this test ist not working
# max_parGlobal = 0.8
# cohort_list = FINN::cohort_df2arrays(cohort_df = data.frame(dbh = 500, nTree = 1, Species = 0))
# FINN:::pkg.env$FINN$height_P(cohort_list$dbh, array(max_parGlobal, dim = 1))
