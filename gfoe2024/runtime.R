library(FINN)
library(data.table)

sites = seq(200, 2000, by = 200)
years = 200
patches = 100L
sp = 10L
NGPUs = c(1, 4, 0)
test_frame = expand.grid(sites = sites, years = years, patches = patches, sp = sp, NGPUs = NGPUs, mean_time = NA, sd_time = NA)

for(K in 1:nrow(test_frame)) {
  tmp = test_frame[K,]
  env = data.table(siteID = rep(1:tmp$sites, each = tmp$years),
                   year = rep(1:tmp$years, tmp$sites),
                   env1 = runif(tmp$sites*tmp$years, -1, 1))

  batchsize = 15L
  parallel = 20L
  device = "cpu"
  if(tmp$NGPUs > 0.5) {
    device = "gpu"
  } else {
    batchsize = round(tmp$sites/2)
    parallel = 2L
  }
  if(tmp$NGPUs == 1) parallel = 5L

  results_time =
    sapply(1:3, function(i) {
      time =
        tryCatch({

          system.time({

              simulateForest(env,
                             sp = sp,
                             patches=tmp$patches,
                             device = device,
                             parallel = parallel,
                             NGPU = tmp$NGPUs,
                             batchsize = batchsize)
          })
        }, error = function(e) e)
      torch::cuda_empty_cache()
      gc()
      torch::cuda_empty_cache()
      if(inherits(time, "error")) return(NA)
      else time[3]
    })
    torch::cuda_empty_cache()
    gc()
    torch::cuda_empty_cache()
    test_frame$mean_time[K] = mean(results_time, na.rm=TRUE)
    test_frame$sd_time[K] = sd(results_time, na.rm=TRUE)
    print(test_frame[K,])
}
saveRDS(test_frame, file = "gfoe2024/results_runtime_21-8-24.RDS" )
