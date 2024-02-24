################################################################################
#' This scripts tests the competition function and contains code that
#' helps to identify plausible parameter ranges
################################################################################

library(FINN)
library(data.table)
library(ggplot2)

# check sensitivity of competition function

# sensitivity to dbh vs parGlobal
dbhTest = seq(0,200, length.out = 200)
i_dbh2height = c(0.5,0.6,0.7)
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
      nTree = 3
    )
  i_cohort_list = cohort_df2arrays(i_cohort_df)
  compOut <- FINN::competition(
    dbh = i_cohort_list$dbh,
    nTree = i_cohort_list$nTree,
    Species = i_cohort_list$Species,
    parGlobal = i_dbh2height
    )
  out_dt <- rbind(
    out_dt,
    data.table(
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
sum(out_dt[parGlobal == i_dbh2height[1]]$AL) < sum(out_dt[parGlobal == i_dbh2height[2]]$AL) &
  sum(out_dt[parGlobal == i_dbh2height[2]]$AL) < sum(out_dt[parGlobal == i_dbh2height[3]]$AL)

# Does a tree with a dbh of 500 has a lower height than 150 meter?
# The tallest tree on earth was 116 m with a dbh of 450 cm.
# FIXME this test ist not working
# max_parGlobal = 0.8
# cohort_list = FINN::cohort_df2arrays(cohort_df = data.frame(dbh = 500, nTree = 1, Species = 0))
# FINN:::pkg.env$FINN$height_P(cohort_list$dbh, array(max_parGlobal, dim = 1))


################################################################################
## test regeneration
################################################################################


# sensitivity to dbh vs parGlobal
dbhTest = seq(0,200, length.out = 500)
envTest = seq(0,10, length.out = 10)
Ncohorts_per_Sp = 4
Nspecies = 6
i_parReg = seq(0.0001,1, length.out = Nspecies)
i_dbh2height = rep(0.7,Nspecies)
out_dt <- data.table()
for(env_i in envTest){

for(dbh_i in dbhTest){
  i_cohort_df =
    data.frame(
      Species = rep(0:(Nspecies-1),each = Ncohorts_per_Sp),
      parGlobal = rep(i_dbh2height, each = Ncohorts_per_Sp),
      dbh = rep(rep(dbh_i,Ncohorts_per_Sp),Nspecies),
      nTree = 1
    )
  # stop()
  i_cohort_list = cohort_df2arrays(i_cohort_df)
  compOut <- FINN::competition(
    dbh = i_cohort_list$dbh,
    nTree = i_cohort_list$nTree,
    Species = i_cohort_list$Species,
    parGlobal = i_dbh2height,
    h = 0
  )
  regOut <- FINN::reg(
    dbh = i_cohort_list$dbh,
    nTree = i_cohort_list$nTree,
    Species = i_cohort_list$Species,
    parGlobal = i_dbh2height,
    parReg = i_parReg,
    pred = c(env_i)
  )
  out_dt <- rbind(
    out_dt,
    data.table(
      dbh = dbh_i,
      parGlobal = rep(i_dbh2height, each = Ncohorts_per_Sp),
      parReg = i_parReg,
      regN = as.vector(regOut),
      AL = as.vector(compOut),
      heights = i_cohort_df$heights,
      env = env_i
    )
  )
}
}


regfun = function(parReg, AL = seq(0,1,0.01)) (plogis((AL + (1-parReg) - 1)/1e-1))
plot(AL, regfun(1), type = "l", col = "red", ylim = c(0,1))
lines(AL, regfun(0.5), type = "l", col = "green")
lines(AL, regfun(0.01), type = "l", col = "blue")

ggplot(out_dt, aes(y = regN, x = cut(AL, breaks = 10)))+
  geom_boxplot()+
  facet_grid(factor(parReg)~factor(env))

################################################################################
## test growth
################################################################################

# sensitivity to dbh vs parGlobal
dbhTest = seq(0,200, length.out = 200)
envTest = seq(0,1, length.out = 10)
Ncohorts_per_Sp = 4
Nspecies = 3
i_parGrowth = seq(0,5, length.out = Nspecies)
# i_dbh2height = rep(0.7,Nspecies)
i_dbh2height = seq(0.5,0.7, length.out = Nspecies)
out_dt <- data.table()
for(env_i in envTest){

for(dbh_i in dbhTest){
  i_cohort_df =
    data.frame(
      Species = rep(0:(Nspecies-1),each = Ncohorts_per_Sp),
      parGlobal = rep(i_dbh2height, each = Ncohorts_per_Sp),
      dbh = rep(rep(dbh_i,Ncohorts_per_Sp),Nspecies),
      nTree = 1
    )
  # stop()
  i_cohort_list = cohort_df2arrays(i_cohort_df)
  compOut <- FINN::competition(
    dbh = i_cohort_list$dbh,
    nTree = i_cohort_list$nTree,
    Species = i_cohort_list$Species,
    parGlobal = i_dbh2height
  )
  growthOut <- FINN::growth(
    dbh = i_cohort_list$dbh,
    nTree = i_cohort_list$nTree,
    Species = i_cohort_list$Species,
    parGlobal = i_dbh2height,
    parGrowth = matrix(c(i_parGrowth,rep(5,length(i_parGrowth))), ncol = 2),
    pred = rep(env_i, length(i_parGrowth))
  )
  out_dt <- rbind(
    out_dt,
    data.table(
      dbh = dbh_i,
      parGlobal = rep(i_dbh2height, each = Ncohorts_per_Sp),
      parGrowth = i_parGrowth,
      growth = as.vector(growthOut),
      AL = as.vector(compOut),
      heights = i_cohort_df$heights,
      env = env_i
    )
  )
}
}

AL = seq(0,1,0.01)
shadeF = function(AL, parGrowth1) 1-(plogis((AL**2)*parGrowth1*10)-0.5)*2
plot(AL, shadeF(AL, 2), type = "l", col = "red", ylim = c(0,2))
lines(AL, shadeF(AL, 0.1), type = "l", col = "green")
lines(AL, shadeF(AL, 10), type = "l", col = "blue")
lines(AL, shadeF(AL, 100), type = "l", col = "blue")

AL = 0
env = 1
AL*(AL+env)

pred = seq(0,1,length.out = 100)
# regfun = function(parReg, AL = seq(0,1,0.01)) (plogis((AL + (1-parReg) - 1)/1e-1))
growfun = function(parGrowth, pred = seq(0,1,length.out = 100)) (1.- (1.- pred)^4.0) * parGrowth
plot(pred, growfun(4), type = "l", col = "red", ylim = c(-0.1,4.1))
lines(pred, growfun(3), type = "l", col = "green")
lines(pred, growfun(0), type = "l", col = "blue")

parMort1 = 30
mortF = function(dbh, parMort) 0.1*(dbh/(parMort*10))^(2.3)
plot(1:300, mortF(1:300,parMort1), type = "l", col = "red")

ggplot(out_dt, aes(y = growth, x = env, color = factor(parGlobal)))+
  geom_jitter(alpha = 0.3)+
  # geom_boxplot()+
  # facet_grid(factor(parGlobal)~factor(env))
  facet_grid(factor(parGrowth)~cut(AL, breaks = 4))


