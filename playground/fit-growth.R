# Create the data frame for the table
data <- data.frame(
  Alter = c(30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150),
  Stammzahl = c(5094, 2348, 1384, 922, 665, 488, 368, 291, 236, 193, 162, 139, 120),
  Mittelhoehe = c(10.4, 15.5, 20.3, 24.2, 27.2, 30.1, 32.6, 34.8, 36.7, 38.5, 39.9, 41.2, 42.2),
  Oberhoehe = c(12.2, 17.3, 21.9, 25.6, 28.4, 31.1, 33.4, 35.3, 37.0, 38.6, 39.9, 41.2, 42.2),
  Grundflaeche = c(15.2, 19.2, 22.4, 25.2, 27.7, 29.6, 30.8, 31.7, 32.4, 32.9, 33.4, 33.9, 34.5),
  Mitteldurchmesser = c(6.2, 10.2, 14.4, 18.6, 23.1, 27.8, 32.6, 37.3, 41.8, 46.6, 51.2, 55.8, 60.5),
  Formzahl = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  Derbholz_Ausersch = c(279, 389, 441, 461, 470, 473, 477, 480, 483, 486, 490, 493, 496),
  Derbholz_Herz = c(44, 116, 200, 281, 355, 421, 479, 530, 576, 615, 653, 689, 722),
  Derbholz_Summe = c(7, 30, 48, 61, 70, 76, 78, 79, 80, 77, 75, 73, 73),
  Derbholz_GWL = c(7.9, 11.4, 12.9, 13.5, 13.6, 12.9, 12.5, 11.9, 11.5, 11.1, 10.6, 10.6, 10.8),
  Derbholz_dGz = c(7, 37, 85, 146, 216, 292, 370, 449, 529, 606, 681, 754, 754),
  Alter_2 = c(30, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, NA),
  Derbholz_dGz2 = c(123, 237, 366, 501, 637, 771, 900, 1025, 1144, 1259, 1370, 1476, NA),
  Alter_3 = c(3.1, 4.7, 6.1, 7.2, 8.0, 8.6, 9.0, 9.3, 9.5, 9.7, 9.8, 9.8, NA)
)

cohort = CohortMat$new(
  obs_df = data.frame(
    siteID = c(1,2),
    patchID = c(1,1),
    cohortID = c(1,1),
    species = c(1,1),
    nTree = c(5094,120),
    dbh = c(6.2,60.5)
  ))

years = 0:150
plot(Mitteldurchmesser~Alter,data, xlim = c(0,150), ylim = c(0,65))
fm_dbh<-lm(log(Mitteldurchmesser) ~ Alter, data = data)
lines(x = years, y = exp(predict.lm(fm_dbh, newdata = data.frame(Alter = years))))
dbh = predict.lm(fm, newdata = data.frame(Alter = years))

plot(Stammzahl~Alter,data, xlim = c(0,150), ylim = c(0,10000))
fm_nTree<-lm(log(Stammzahl) ~ Alter, data = data)
lines(x = years, y = exp(predict.lm(fm_nTree, newdata = data.frame(Alter = years))))
nTree = exp(predict.lm(fm_nTree, newdata = data.frame(Alter = years)))

data.frame(
  year = years[-1],
  growthRate = diff(dbh),
  mortRate = -diff(nTree),
  regRate = NA,
  dbh = dbh[-1],
  nTree = nTree[-1]
)
