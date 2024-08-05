library(torch)
library(FINN)

test_gradients_regeneration = function(
      n_species = 10,
      n_sites = 10,
      n_patches = 5,
      n_cohort = 3,
      multiplier_parReg = 5,
      multiplier_pred = 5,
      multiplier_AL = 1) {

  torch::with_detect_anomaly({

    species = torch::torch_randint(1, n_species, size = c(n_sites, n_patches, n_cohort), dtype = torch::torch_int64())
    parReg = torch::torch_randn(size = c(n_species), requires_grad = TRUE)
    parReg_multiplied = parReg*multiplier_parReg
    pred = torch::torch_rand(size = c(n_sites, n_species), requires_grad = TRUE)
    pred_multiplied = pred*multiplier_pred
    light = torch::torch_randn(size = c(n_sites, n_patches, 1))*multiplier_AL

    regeneration(species, parReg_multiplied, pred_multiplied, light, patch_size_ha= 0.1)$sum()$backward()
    grad_par = as.matrix(parReg$grad)
    grad_pred = as.matrix(pred$grad)

  })
  return(list(grad_par, grad_pred))
}

grads = test_gradients_regeneration(n_species = 50L, n_cohort = 20L, multiplier_parReg = 100.001, multiplier_pred = 0.1)
hist(as.vector(grads[[1]]))



test_gradients_mortality = function(
    n_species = 10,
    n_sites = 10,
    n_patches = 5,
    n_cohort = 3,
    lambda_trees = 1,
    multiplier_dbh = 100,
    multiplier_parMort = 5,
    multiplier_AL = 1) {

  torch::with_detect_anomaly({

    species = torch::torch_randint(1, n_species,
                                   size = c(n_sites, n_patches, n_cohort),
                                   dtype = torch::torch_int64())
    trees = torch::torch_tensor(array(rpois(prod(c(n_sites, n_patches, n_cohort)),
                                            lambda = lambda_trees),
                                      dim = c(n_sites, n_patches, n_cohort)) )
    dbh = torch::torch_rand(size = c(n_sites, n_patches, n_cohort))*multiplier_dbh
    parMort = torch::torch_randn(size = c(n_species, 2), requires_grad = TRUE)
    parMort_multiplied = parMort*multiplier_parMort
    pred = torch::torch_rand(size = c(n_sites, n_species), requires_grad = TRUE)
    light = torch::torch_rand(size = c(n_sites, n_patches, n_cohort))*multiplier_AL

    mortality(dbh, species,trees, parMort_multiplied, pred, light)$sum()$backward()

    grad_par = as.matrix(parMort$grad)
    grad_pred = as.matrix(pred$grad)

  })
  return(list(grad_par, grad_pred))
}

grads = test_gradients_mortality(lambda_trees = 1, multiplier_parMort = 50, multiplier_AL = 0.0)
hist(as.vector(grads[[1]]))



test_gradients_growth = function(
    n_species = 10,
    n_sites = 10,
    n_patches = 5,
    n_cohort = 3,
    multiplier_dbh = 100,
    multiplier_parMort = 5,
    multiplier_parGrowth = 5,
    multiplier_AL = 1) {

  torch::with_detect_anomaly({

    species = torch::torch_randint(1, n_species,
                                   size = c(n_sites, n_patches, n_cohort),
                                   dtype = torch::torch_int64())
    dbh = torch::torch_rand(size = c(n_sites, n_patches, n_cohort))*multiplier_dbh
    parMort = torch::torch_randn(size = c(n_species, 2), requires_grad = TRUE)
    parMort_multiplied = parMort*multiplier_parMort
    parGrowth = torch::torch_randn(size = c(n_species, 2), requires_grad = TRUE)
    parGrowth_multiplied = parGrowth*multiplier_parGrowth
    pred = torch::torch_rand(size = c(n_sites, n_species), requires_grad = TRUE)
    light = torch::torch_rand(size = c(n_sites, n_patches, n_cohort))*multiplier_AL

    growth(dbh, species,  parGrowth_multiplied, parMort_multiplied, pred, light)$sum()$backward()

    grad_par = cbind(as.matrix(parMort$grad), as.matrix(parGrowth$grad))
    grad_pred = as.matrix(pred$grad)

  })
  return(list(grad_par, grad_pred))
}

grads = test_gradients_growth( multiplier_parMort = 5,multiplier_parGrowth = 1, multiplier_AL = 1.0)
hist(as.vector(grads[[1]][,4]))



test_gradients_competition = function(
    n_species = 10,
    n_sites = 10,
    n_patches = 5,
    n_cohort = 3,
    lambda_trees = 5,
    multiplier_dbh = 100,
    multiplier_parHeight = 5,
    h = 1) {

  torch::with_detect_anomaly({

    species = torch::torch_randint(1, n_species,
                                   size = c(n_sites, n_patches, n_cohort),
                                   dtype = torch::torch_int64())
    trees = torch::torch_tensor(array(rpois(prod(c(n_sites, n_patches, n_cohort)),
                                            lambda = lambda_trees),
                                      dim = c(n_sites, n_patches, n_cohort)) )
    dbh = torch::torch_rand(size = c(n_sites, n_patches, n_cohort))*multiplier_dbh
    parHeight = torch::torch_randn(size = c(n_species), requires_grad = TRUE)
    parHeight_multiplied = parHeight*multiplier_parHeight

    competition(dbh, species, trees, parHeight, h, patch_size_ha = 0.1)$sum()$backward()

    grad_par = cbind(as.matrix(parHeight$grad))

  })
  return(list(grad_par))
}

grads = test_gradients_competition(multiplier_dbh = 100,lambda_trees = 100, h = NULL)
hist(as.vector(grads[[1]]))

