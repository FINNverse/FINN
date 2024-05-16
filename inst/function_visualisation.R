# Load necessary library
library(ggplot2)
library(data.table)
################################################################################
## create pdf file
################################################################################
# Open the existing PDF file in append mode
pdf("function_vis.pdf")

################################################################################
## sensitivity of height_P function
################################################################################
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

test_dt[, height := height_P(dbh, parHeight), by = rowID]

p_height_P <-
  ggplot(test_dt, aes(x = dbh, y = height, group = parHeight,color = parHeight))+
    ggtitle("height_P",paste0(deparse(height_P),collapse = "\n"))+
    ylab("height_P(dbh,parHeight)")+
    xlab("dbh [cm]")+
    geom_line()+
    geom_label(aes(label = parHeight, color = parHeight), test_dt[,.(height = max(height),dbh = max(dbh)),by=parHeight])
print(p_height_P)

################################################################################
## sensitivity of BA_P function
################################################################################

# Generate test data
test_dt <- data.table(
  dbh = seq(0, 500, 1)
)

# Calculate basal area
test_dt[, basal_area := BA_P(dbh)]

# Plot
p_BA_P <- ggplot(test_dt, aes(x = dbh, y = basal_area)) +
  ggtitle("BA_P", paste0(deparse(BA_P), collapse = "\n")) +
  ylab("BA_P(dbh)") +
  xlab("dbh [cm]") +
  geom_line()
print(p_BA_P)

#
# ################################################################################
# ## sensitivity of compF_P function
# ################################################################################
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

ggplot(cohort_df1,aes(y = factor(trees/0.1), x = dbh, fill = light))+
  geom_tile()+
  ylab("trees/ha [N]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
  ggtitle("competition",paste0(deparse(competition),collapse = "\n"))

#
# Close the PDF device
dev.off()
