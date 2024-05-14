# Load necessary library
library(ggplot2)

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
#
#   cohort = CohortMat$new(obs_df = data.frame(
#     siteID = 1,
#     patchID = 1,
#     cohortID = 1,
#     species = 1,
#     nTree = 1,
#     dbh = 100
#   ))
#
# cohort$asDF()
#
#   data.frame(
#     dbh = runif(Ncohorts,1,200),
#     nTree = sample.int(10,Ncohorts,replace = T)
#   )
#   array(1:12,dim = c(2,3,2))
#
#
#
#   test_cohorts = CohortMat$new(dbh = 1:3, nTree = 4:7,Species = rep(1,3) )
#   # dims:  sites, patches, cohorten
#   test_cohorts$dbh
#   test_cohorts$nTree
#   test_cohorts$Species
#   compF_P(torch_tensor(dbh), torch_tensor(parHeight),
#           torch_tensor(nTree), parHeight = torch_tensor(0.5))
# }
# CohortMat(1)
# compF_P()
# test_dt[, nTree := 100]
# test_dt[, parHeight := 0.5]
# test_dt$rowID = 1:nrow(test_dt)
# test_dt[, AL := as.numeric(compF_P(torch_tensor(dbh), torch_tensor(Species), torch_tensor(nTree), torch_tensor(parHeight))), by = rowID]
#
#
#
# # Plot
# p_compF_P <- ggplot(test_dt, aes(x = dbh, y = AL, group = Species, color = as.factor(Species))) +
#   ggtitle("compF_P", paste0(deparse(compF_P), collapse = "\n")) +
#   ylab("compF_P(dbh,Species,nTree,parHeight)") +
#   xlab("dbh [cm]") +
#   geom_line() +
#   geom_label(aes(label = Species, color = as.factor(Species)),
#              test_dt[, .(AL = max(AL), dbh = max(dbh)), by = Species])
# print(p_compF_P)
#
# Close the PDF device
dev.off()
