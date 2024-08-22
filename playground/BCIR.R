library(data.table)
df = fread("~/Downloads/stand_dt.csv")
head(df)

cohort = fread("~/Downloads/obs_df.csv")
obs_array <- FINN::obsDF2arrays(cohort)
cohort1 <- FINN::CohortMat$new(dbh = obs_array$dbh, trees = obs_array$trees, species = obs_array$species)
cohort1$trees$max()

env = fread("/Users/maximilianpichler/Downloads/env_dt.csv")

head(df)
df
df$AL = NA
colnames(df)[c(1,7, 9, 11)] = c("year", "growth", "reg", "mort")
head(df)
data = df

# fill up missing sites
# fill sites

empty_res =
  lapply(unique(data$species), function(sp) {
    empty_sites =
      lapply(unique(data$year), function(yy) {
          expand.grid(siteID = setdiff(data$siteID |> unique(), data[species==sp & year == yy]$siteID),
                      year = yy,
                      species = sp
          )
      })
    rbindlist(empty_sites, fill = TRUE)
  })
data2 = rbindlist(list(data, rbindlist(empty_res, fill = TRUE)), fill = TRUE)
data = data2

missing_years =
lapply(1982:1984, function(year) {
  tmp = env[2,]
  tmp$year = year
  return(tmp)
  }) |> rbindlist()

env = rbindlist(list(missing_years, env[-1,]))
env = env[!year %in% 2016:2019]

env =
  lapply(unique(data$siteID), function(site) {
    tmp = env
    tmp$siteID = site
    return(tmp)
  }) |> rbindlist()

head(env)
head(data)
env$Prec = scale(env$Prec)
env$T_mean = scale(env$T_mean)
env$T_max = scale(env$T_max)
env$T_min = scale(env$T_min)
m = finn(env = env, data = data)
