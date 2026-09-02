# env_l2 sweep fit (fresh process) -> overshoot + census fit. Patched FINN (~/Rlib-finnfix).
local({ l <- file.path(path.expand("~"), "Rlib-finnfix"); if (dir.exists(l)) .libPaths(c(l, .libPaths())) })
suppressWarnings(suppressMessages({ library(data.table); library(FINN) }))
source("dev/pft-bci/finn_membership.R"); set.seed(1); FINN.seed(1)
D <- "dev/pft-bci"; RES <- file.path(D, "results"); dir.create(RES, showWarnings = FALSE)
COND   <- Sys.getenv("COND", "pft5")
ENV_L2 <- as.numeric(Sys.getenv("ENV_L2", "0"))
EPOCHS <- as.integer(Sys.getenv("EPOCHS", "2000"))
PATCH  <- as.integer(Sys.getenv("PATCHES", "25"))
is5 <- COND == "pft5"
obs <- fread(file.path(D, if (is5) "data/pft5/obs_dt.csv" else "data/obs_species.csv"))
env <- fread(file.path(D, if (is5) "data/pft5/env_dt.csv" else "data/env.csv"))
coh <- fread(file.path(D, if (is5) "data/pft5/initial_cohorts1985.csv" else "data/initial_cohorts1985.csv"))
Nsp <- max(obs$species); cens <- c(5,10,15,20,25,30)
m <- finn_membership(mm_K=5L, N_species=Nsp, recruits_dbh=1.0,
  regeneration_saturation=list(K_init=200,shared=TRUE,bounds=c(1,3000)),
  competition_process=createProcess(~0,FINN::competition,optimizeSpecies=TRUE),
  growth_process=createProcess(~.,FINN::growth,optimizeSpecies=TRUE,optimizeEnv=TRUE),
  regeneration_process=createProcess(~.,FINN::regeneration,optimizeSpecies=TRUE,optimizeEnv=TRUE),
  mortality_process=createProcess(~.,FINN::mortality,optimizeSpecies=TRUE,optimizeEnv=TRUE))
if (ENV_L2 > 0) m$env_l2 <- ENV_L2
ch <- FINN::CohortMat(obs_df=coh, sp=Nsp)
fit(m, data=obs, env=env, init_cohort=ch, device="cpu", epochs=EPOCHS,
    batchsize=min(15L, uniqueN(obs$siteID)), patches=PATCH, lr=0.01,
    optimizer=torch::optim_adam, env_autoscale=TRUE, plot_progress=FALSE,
    weights=c(0.1,10,1,1,1,1), loss=c(dbh="mse",ba="mse",trees="nbinom",growth="mse",mortality="mse",regeneration="nbinom"))
m$eval()
w <- as.data.table(predict(m, env=env, init_cohort=ch, patches=PATCH, device="cpu")$long$site)
ba <- w[variable=="ba", .(ba=sum(value,na.rm=T)), by=.(siteID,year)][,.(sim=mean(ba)),by=year]
ov <- max(ba[!year%in%cens]$sim)/max(ba[year%in%cens]$sim)
ob <- obs[,.(obs=sum(ba,na.rm=T)),by=.(siteID,year)][,.(obs=mean(obs)),by=year]
cmp <- merge(ba[year%in%cens], ob, by="year")
line <- sprintf("COND=%s ENV_L2=%s overshoot=%.2f peakBA=%.2f censusBA_cor=%.2f\n",
                COND, ENV_L2, ov, max(ba$sim), suppressWarnings(cor(cmp$sim, cmp$obs)))
cat("RESULT", line)
cat(line, file = file.path(RES, "envl2_sweep.txt"), append = TRUE)
