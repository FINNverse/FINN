library(ProfoundData)
library(data.table)
library(sf)
library(mapview)

myDB <- path.expand("data/calibration-data/PROFOUND/ProfoundData.sqlite")
setDB(myDB)

overview <- data.table(browseData())

summarizeData(dataset = "CLIMATE_LOCAL", site = "bily_kriz", mode = "overview")

sites_dt <- data.table(getData(dataset = "SITES"))
sitedescription_dt <- data.table(getData(dataset = "SITEDESCRIPTION"))

stand_dt = data.table()
for(i_site in overview[STAND == 1]$site){
  stand_dt <- rbind(
    stand_dt,
    data.table(getData(dataset = "STAND", site = i_site)),
    fill = T
  )
}

tree_dt = data.table()
for(i_site in overview[TREE == 1]$site){
  tree_dt <- rbind(
    tree_dt,
    data.table(getData(dataset = "TREE", site = i_site)),
    fill = T
  )
}

unique(tree_dt[,.(site_id, size_m2, year)]) %>% View

tree_dt[,periodID := as.integer(factor(year, ordered = T)), by = .(site_id)]
stand_dt[,periodID := as.integer(factor(year, ordered = T)), by = .(site_id)]

obs_dt = tree_dt[periodID == 1, .(siteID = site_id, dbh = dbh1_cm, species)]

library(ggplot2)

plot_dt = obs_dt[,.(count = .N),
                 by = .(
                   siteID,
                   species,
                   dbh = as.numeric(as.character(cut(dbh,breaks = seq(0,100,2), labels = seq(2,100,2)))))]

plot_dt <- merge(plot_dt, unique(stand_dt[periodID == 1,.(
  siteID = site_id,
  dbh_mean = dbhArith_cm,
  dbh_weightedmean = dbhBA_cm,
  species)]), by = c("siteID","species"))
ggplot(plot_dt, aes(x = dbh, y = count))+
  # add stacked bars for count of each species in dbh class
  geom_bar(aes(fill = species), stat = "identity", position = "stack")+
  # add horizontal error bars and mean estimate for each species
  geom_point(aes(x = dbh_mean, y = 200, color = "mean"))+
  geom_point(aes(x = dbh_weightedmean, y = 200, color = "weighted mean"))+
  facet_grid(species~siteID)

FINN::CohortMat$new(tree_dt[periodID == 1, .(siteID = site_id, dbh = cut(), species)])
CohortMat$new(tree_dt, stand_dt, sites_dt, sitedescription_dt)

summary(table(tree_dt$record_id))


tree_dt = getData(dataset = "TREE")


names(sites_dt)
shp_sites <- st_as_sf(sites_dt, coords = c("lon", "lat"), crs = 4326)
#mapview(shp_sites, zcol = "site")


