# Load necessary library
library(ggplot2)
library(data.table)
library(torch)
library(FINN)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## create pdf ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# Open the existing PDF file in append mode
# pdf("function_vis.pdf")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## height ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
test_dt =
  data.table(
    expand.grid(
      list(
        dbh = seq(0, 500,1),
        parHeight = seq(0,1,0.1)
      )
    )
  )
test_dt$rowID = 1:nrow(test_dt)

test_dt[, height := height(dbh, parHeight), by = rowID]

p_height <-
  ggplot(test_dt, aes(x = dbh, y = height, group = parHeight,color = parHeight))+
    ggtitle("height",paste0(deparse(height),collapse = "\n"))+
    ylab("height(dbh,parHeight)")+
    xlab("dbh [cm]")+
    geom_line()+
    geom_label(aes(label = parHeight, color = parHeight), test_dt[,.(height = max(height),dbh = max(dbh)),by=parHeight])
print(p_height)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## BA_stem ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# Generate test data
test_dt <- data.table(
  dbh = seq(0, 500, 1)
)

# Calculate basal area
test_dt[, basal_area := BA_stem(dbh)]

# Plot
p_BA_stem <- ggplot(test_dt, aes(x = dbh, y = basal_area)) +
  ggtitle("BA_stem", paste0(deparse(BA_stem), collapse = "\n")) +
  ylab("BA_stem(dbh)") +
  xlab("dbh [cm]") +
  geom_line()
print(p_BA_stem)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## BA_stand ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
trees_vec = c(0:500,10^(seq(2,4, length.out = 20)))
dbh_vec = dbh = seq(1,300, 1)
patch_size_vec = seq(0.04, 0.2, 0.02)
patch_size_vec = c(0.01,0.05,0.08,seq(0.1,1,0.2),1)

test_dt_out <- data.table()
i=patch_size_vec[1]
for(i in patch_size_vec){
  # Generate test data
  cohort_df1 =
    data.frame(
      patchID = 1,
      cohortID = 1,
      species = 1,
      expand.grid(
        trees_ha = trees_vec,
        patch_size_ha = i,
        dbh = dbh_vec
      )
    )
  cohort_df1$siteID = 1:nrow(cohort_df1)
  cohort_df1$trees = round(cohort_df1$trees_ha*i)
  cohort = CohortMat$new(obs_df = cohort_df1)

  basal_area = BA_stand(cohort$dbh, cohort$trees, patch_size_ha = i)
  cohort_df1$basal_area = torch::as_array(basal_area)[,1,1]
  cohort_df1$patch_size_ha = i
  test_dt_out <- rbind(test_dt_out,cohort_df1)
}

max_ba = 60
test_dt_out[basal_area > max_ba, basal_area := max_ba,]
p_BA_stand <- ggplot(test_dt_out[trees_ha <=500],aes(y = factor(trees_ha), x = dbh, fill = basal_area))+
  geom_tile()+
  ylab("trees [N/ha]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="basal area\n[m^2/ha]"))+
  ggtitle("BA_stand",paste0(deparse(BA_stand),collapse = "\n"))+
  scale_fill_viridis_c(breaks = seq(0,max_ba,10), labels = c(seq(0,max_ba-10,10),paste0(">",max_ba)))+
  scale_y_discrete(breaks = seq(0,500,100))+
  facet_wrap(~patch_size_ha, labeller = labeller(patch_size_ha = function(x) paste0("patch size = ",x, " ha")))
print(p_BA_stand)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## competition ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
patch_size = 0.1

trees_vec = c(1:99,10^(seq(2,4, length.out = 16)))
dbh_vec = dbh = seq(1,300, 1)

cohort_df1 =
  data.frame(
    patchID = 1,
    cohortID = 1,
    species = 1,
    expand.grid(
      trees = round(trees_vec*patch_size),
      dbh = dbh_vec
    )
  )
cohort_df1$siteID = 1:nrow(cohort_df1)

cohort = CohortMat$new(obs_df = cohort_df1)

dbh = cohort$dbh
species = cohort$species
trees = cohort$trees

comp = competition(cohort$dbh, cohort$species, cohort$trees,
                   parHeight = torch::torch_tensor(0.5), h=0)

cohort_df1$light = torch::as_array(comp)[,1,1]

ggplot(cohort_df1,aes(y = factor(trees/0.1), x = dbh, fill = light*100))+
  geom_tile()+
  ylab("trees/ha [N]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
  ggtitle("competition",paste0(deparse(competition),collapse = "\n"))

#
# Close the PDF device
dev.off()
