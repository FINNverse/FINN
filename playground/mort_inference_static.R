library(torch)
library(coro)

load(file = "mort_data.RData")

net <- nn_module(
  "MonthlyRNN",
  initialize = function(input_size, hidden_size, num_layers, output_size) {
    self$fc1 = nn_linear(8, 50)
    self$fc2 = nn_linear(50, 50)
    self$fc3 = nn_linear(50, 169)
  },
  forward = function(x) {
    out <- nnf_selu(self$fc1(x))
    out <- nnf_selu(self$fc2(out))
    out <- nnf_sigmoid(self$fc3(out))
    return(out)
  }
)

time_frame = 685:744


tmp = cohort_array$dbh
tmp[is.na(tmp)] = 1

fill_na = function(x, fill = 1) {
  x[is.na(x)] = fill
  return(x)
}

dbh = torch::torch_tensor(fill_na(cohort_array$dbh))
trees = torch::torch_tensor(fill_na(cohort_array$trees, 0))
species = torch::torch_tensor(fill_na(cohort_array$species, 1), dtype = torch::torch_int64())

climate_array[,,1] = (climate_array[,,1] - mean(climate_array[,,1], na.rm = TRUE) ) / sd(climate_array[,,1], na.rm = TRUE)
climate_array[,,2] = (climate_array[,,2] - mean(climate_array[,,2], na.rm = TRUE) ) / sd(climate_array[,,2], na.rm = TRUE)
climate_array[is.na(climate_array)] = 0

dim(climate_array)
climate_array2 = cbind(apply(climate_array[,time_frame,], c(1, 3), mean), apply(climate_array[,time_frame,], c(1, 3), max), apply(climate_array[,time_frame,], c(1, 3), min), apply(climate_array[,time_frame,], c(1, 3), sd))

env = torch_tensor(climate_array2)
Y = torch::torch_tensor(fill_na(cohort_array$y))

model = net()
model$cuda(device = "cuda:0")

parHeightT = torch::torch_tensor(fill_na(parHeight, 1), device = "cuda:0")
parMort = torch::torch_tensor(matrix(c(0.25, 100), as.integer(max(species)), 2, byrow = TRUE), requires_grad = TRUE, device = "cuda:0")


dataset = torch::tensor_dataset(dbh, trees, species, env, Y)
dataloader = torch::dataloader(dataset, batch_size = 2000L, shuffle = TRUE, pin_memory = TRUE)

opt = optim_adagrad(params = c(model$parameters), lr = .01)

batch = coro::as_iterator(dataloader)()

mortality_pred = function(dbh, species, trees, parMort, pred, light, debug = F) {
  shade = 1-torch_sigmoid((light + (1-parMort[,1][species]) - 1)/(1/10^(1.5 + torch_abs(light-0.5))))
  environment = (index_species(pred, species))
  gPSize = 0.1*(dbh/torch_clamp((parMort[,2][species]*100), min = 0.00001))$pow(2.3)
  gPSize = 2*(torch_sigmoid(gPSize)-0.5)
  predM = torch_sigmoid(environment*(shade+gPSize) + shade*gPSize + shade + gPSize)
  return(predM)
}

library(FINN)

for(epoch in 1:100) {
  counter = 1
  b_loss = rep(0, 1e6)
  coro::loop(for (batch in dataloader) {
    dbh_b = batch[[1]]$to(device = "cuda:0")
    trees_b = batch[[2]]$to(device = "cuda:0")
    species_b = batch[[3]]$to(device = "cuda:0")
    env_b = batch[[4]]$to(device = "cuda:0")
    Y_b = batch[[5]]$to(device = "cuda:0")

    opt$zero_grad()
    pred_env = model(env_b)$squeeze()
    AL = FINN::competition(dbh = dbh_b, trees = trees_b, species = species_b, parHeight = parHeightT, h = NULL, patch_size_ha = 0.4)
    pred_prob = mortality_pred(dbh_b, species_b, trees_b, parMort, pred_env, AL)
    loss = (torch::distr_bernoulli(probs = pred_prob)$log_prob(Y_b)*(trees_b$ge_(0.5)))$sum(dim=3)$mean(dim=1)$negative()
    loss$backward()
    opt$step()
    b_loss[counter] = loss$item()
    counter = counter + 1
  })
  print(sum(b_loss))
}
as.matrix(pred_env[,7]$cpu())


df = data.frame()
for(j in c(1, 23, 30, 34)) {
  for(i in 1:20419){
    tmp_true = fill_na(cohort_array$species[i,,], 0) == j
    if(sum(tmp_true)> 0 ) {
      cohort_array$trees[i,,][tmp_true]
      died = cohort_array$y[i,,][tmp_true]
      dbh = cohort_array$dbh[i,,][tmp_true]
      df = rbind(df,  cbind(data.frame(died = died, site = i, species = j, dbh = dbh), matrix(climate_array2[i,], nrow = sum(tmp_true), ncol = 8L, byrow = TRUE)))
    }
  }
}

colnames(df)[5:12] = LETTERS[1:8]
colnames(df)[3] = "species"
colnames(df)[4] = "dbh"
library(dplyr)
tmp_df =
df %>% group_by(site, species) %>% summarise(k = length(died), n = sum(died), DBH = mean(dbh), A=mean(A), B = mean(B), C = mean(C), D = mean(D), E = mean(E), F = mean(F), G = mean(G), H = mean(H))



tmp_df
library(lme4)
mm3 = glm(cbind(n, k-n)~A+B+C+D+E+F+G+H+DBH, data = tmp_df %>% filter(species == 25), family = binomial())
summary(mm3)
car::Anova(mm3)


mm2= glm(died~A+B+C+D+E+F+G+H+dbh, data = df %>% filter(species == 25), family = binomial())
summary(mm2)





sp_43 = df %>% filter(species == 30)
table(sp_43$died)
dbh = scale(sp_43$dbh)

sub_climate = climate_array[sp_43$site,time_frame,]

net <- nn_module(
  "MonthlyRNN",
  initialize = function(input_size, hidden_size, num_layers, output_size) {
    self$rnn <- nn_lstm(input_size, hidden_size, num_layers, batch_first = TRUE)
    self$fc <- nn_linear(hidden_size, output_size)
  },
  forward = function(x) {
    out <- self$rnn(x)
    out <- self$fc(out[[1]])[,37,]
    out$sigmoid()
  }
)


env = torch_tensor(sub_climate)
dbh = torch_tensor(matrix(dbh, ncol = 1L))
Y = torch::torch_tensor(sp_43$died, dtype = torch::torch_float32())



library(FINN)

net <- nn_module(
  "MonthlyRNN",
  initialize = function(input_size, hidden_size, num_layers, output_size) {
    self$rnn <- nn_lstm(input_size, hidden_size, num_layers, batch_first = TRUE, bidirectional = FALSE)
    self$fc1 <- nn_linear(hidden_size+1, hidden_size)
    self$fc2 <- nn_linear(hidden_size, output_size)

  },
  forward = function(x, dbh) {
    out <- self$rnn(x)[[1]][,x$shape[2],]
    out2 = torch_cat(list(out, dbh), dim = 2L)
    out = nnf_selu(self$fc1(out2))
    out = self$fc2(out)
    return(torch_sigmoid(out))
  }
)




model = net(3L, 100L, 2L, 1L)
model$cuda(device = "cuda:0")

ww = (1/(table(as.numeric(Y$cpu())) / sum(table(as.numeric(Y$cpu())) ))[as.character(as.numeric(Y$cpu()))])


dim(env)
env_de = as.array(env)
for(i in 1:dim(env_de)[1]) {
  for(j in 1:2) {
    env_de[i,,j] = detrend(env_de[i,,j])
  }
}

env_de[,,1] = env_de[,,1] - mean(env_de[,,1])
env_de[,,2] = env_de[,,2] - mean(env_de[,,2])
env = torch_tensor(env_de)

dbh_time = dbh$expand(list(dbh$shape[1], 60))$unsqueeze(3)
env = torch_cat(list(env, dbh_time), dim = 3L)

dataset = torch::tensor_dataset(torch_tensor(env), Y, dbh, wwT)
dataloader = torch::dataloader(dataset, batch_size = 500L, shuffle = TRUE, pin_memory = TRUE)

opt = optim_rmsprop(params = c(model$parameters), lr = .001)

batch = coro::as_iterator(dataloader)()

for(epoch in 1:1000) {
  counter = 1
  b_loss = rep(0, 1e6)
  coro::loop(for (batch in dataloader) {
    env_b = batch[[1]]$to(device = "cuda:0")
    Y_b = batch[[2]]$to(device = "cuda:0")
    dbh_b = batch[[3]]$to(device = "cuda:0")
    weights = batch[[4]]$to(device = "cuda:0")

    opt$zero_grad()
    pred_prob = model(env_b, dbh_b)$squeeze()
    loss = (torch::distr_bernoulli(probs = pred_prob)$log_prob(Y_b))*weights
    loss = loss$negative()$mean()
    loss$backward()
    opt$step()
    b_loss[counter] = loss$item()
    counter = counter + 1
  })
  cat("Epoch: ", epoch, "Loss: ", sum(b_loss), "AUC: ", Metrics::auc( as.matrix(Y_b$cpu()),as.matrix(pred_prob$cpu())), "\n")
}
Metrics::auc( as.matrix(Y$cpu()),as.matrix(model(env$to(device = "cuda:0"), dbh$to(device = "cuda:0"))$cpu()))


predictions = as.matrix(model(env$to(device = "cuda:0"), dbh$to(device = "cuda:0"))$cpu())
predictions
Y_hat = as.matrix(Y$cpu())



pp = as.matrix(model(torch::torch_zeros_like( env$to(device = "cuda:0") ), dbh$to(device = "cuda:0"))$cpu())
plot(as.matrix(dbh)* attr(scale(sp_43$dbh), "scaled:scale") +  attr(scale(sp_43$dbh), "scaled:center"), pp, xlab = "dbh")


E = as.array(env)
d = as.matrix(dbh)

dim(E)
E2 = cbind(matrix(t(E[,,1]), ncol = 1L), matrix(t(E[,,2]), ncol = 1L), matrix(t(E[,,3]), ncol = 1L))
dbh2 = rep(d[,1], each = length(time_frame))

new_df = cbind(E2, dbh2)
head(new_df)
predict_f = function(model, newdata) {
  e_tmp = newdata[,c(1:3)]
  dbh_tmp = newdata[,4]

  e_tmp = (abind::abind(list(matrix(e_tmp[,1], ncol = length(time_frame), byrow = TRUE), matrix(e_tmp[,2], ncol = length(time_frame), byrow = TRUE), matrix(e_tmp[,3], ncol = length(time_frame), byrow = TRUE)), along = 3L))
  dbh_tmp = matrix(dbh_tmp, ncol = length(time_frame), byrow = TRUE)[,1,drop=FALSE]
  pp = as.matrix(model( torch::torch_tensor(e_tmp, device = "cuda:0") , torch::torch_tensor(dbh_tmp, device = "cuda:0"))$cpu())
  return(pp[,1])
}

res = marginalEffectsGeneric(model, new_df, predict_func = predict_f)
mes = res[[1]]
dim(mes)
fields::image.plot( t(matrix(mes[,2,2], ncol = length(time_frame), byrow = TRUE)[4:7,]) )
Y_hat[392:400]
predictions[392:400]






indices = order(as.matrix(dbh$cpu()), decreasing = TRUE)

head(cbind(Y_hat, predictions, as.matrix(dbh$cpu()[,1]))[indices,], n = 250)
ind = indices[c(37:40)]
Y_hat[ind, ]
predictions[ind,]



par(mfrow = c(1,4))

for(K in 1:4) {
  site = ind[K]
  E_tmp = E
  #E_tmp[,,1] = 0
  preds = sapply(1:length(time_frame), function(k) as.numeric(model( torch::torch_tensor(E_tmp, device = "cuda:0")[site,1:k,,drop=FALSE] , dbh$to( device = "cuda:0")[site,,drop=FALSE])$cpu()))
  plot(1:length(time_frame), (E[site,,1]), col = "red", type = "l", main = Y_hat[ind, ][K], xlab = "Month", ylim = c(-2, 3))
  abline(v = seq(0, 60, by = 12), col = "grey")
  abline(h = 0)
  #points(1:length(time_frame), matrix(mes[,1,1], ncol = length(time_frame), byrow = TRUE)[site,]* E[site,,1], type = "l", col = "blue")
  points(1:length(time_frame), (E[site,,2]), col = "blue", type = "l", main = "Died", xlab = "Month", ylim = c(-6, 6))
  points(1:length(time_frame), 1*( preds ), type = "l", col = "green")


}

site = 2019
preds = sapply(1:length(time_frame), function(k) as.numeric(model( torch::torch_tensor(E_tmp, device = "cuda:0")[site,1:k,,drop=FALSE] , dbh$to( device = "cuda:0")[site,,drop=FALSE])$cpu()))
plot(1:length(time_frame), E[site,,2], col = "red", type = "l", main = "Died", xlab = "Month", ylim = c(-6, 6))
points(1:length(time_frame), matrix(mes[,2,2], ncol = length(time_frame), byrow = TRUE)[site,]* E[site,,2], type = "l", col = "blue")
points(1:length(time_frame), binomial()$linkfun( preds ), type = "l", col = "green")
legend("topright", col = c("blue", "red", "green"), legend = c("Mort increase", "Temp", "Mort Prob"), pch = 15, bty = "n")



par(mfrow = c(1,2))

E_tmp = E
#E_tmp[,,1] = 0
site = 892
preds = sapply(1:length(time_frame), function(k) as.numeric(model( torch::torch_tensor(E_tmp, device = "cuda:0")[site,1:k,,drop=FALSE] , dbh$to( device = "cuda:0")[site,,drop=FALSE])$cpu()))
plot(1:length(time_frame), detrend(E[site,,1]), col = "red", type = "l", main = "Survived", xlab = "Month", ylim = c(-1, 3))
abline(v = seq(0, 60, by = 12), col = "grey")
abline(h = 0)
#points(1:length(time_frame), matrix(mes[,1,1], ncol = length(time_frame), byrow = TRUE)[site,]* E[site,,1], type = "l", col = "blue")
points(1:length(time_frame), detrend(E[site,,2]), col = "blue", type = "l", main = "Died", xlab = "Month", ylim = c(-6, 6))
points(1:length(time_frame), ( preds ), type = "l", col = "green")

site = 2019

preds = sapply(1:length(time_frame), function(k) as.numeric(model( torch::torch_tensor(E_tmp, device = "cuda:0")[site,1:k,,drop=FALSE] , dbh$to( device = "cuda:0")[site,,drop=FALSE])$cpu()))
plot(1:length(time_frame), detrend(E[site,,1]), col = "red", type = "l", main = "Died", xlab = "Month", ylim = c(-1, 3))
abline(v = seq(0, 60, by = 12), col = "grey")
abline(h = 0)

points(1:length(time_frame), detrend(E[site,,2]), col = "blue", type = "l", main = "Died", xlab = "Month", ylim = c(-6, 6))
#points(1:length(time_frame), matrix(mes[,1,1], ncol = length(time_frame), byrow = TRUE)[site,]* E[site,,2], type = "l", col = "blue")
points(1:length(time_frame), binomial()$linkfun( preds ), type = "l", col = "green")
legend("topright", col = c("blue", "red", "green"), legend = c("Precip", "Temp", "Mort Prob"), pch = 15, bty = "n")

detrend = function(env) {
 sin_model <- lm(precipitation ~ sin(2 * pi * month_num / 12) + cos(2 * pi * month_num / 12), data = data.frame(month_num = c(rep(1:12, length(env)/12)), precipitation = env))
 fitted <- predict(sin_model)
 return( env - fitted + coef(sin_model)[1] )
}






breaks = seq(0, 1, length.out = 6)
breaks = cut(preds, breaks = breaks)
cols = colorRampPalette(c("blue", "red"))(6)
plot(E[site,,1], E[site,,2], col = cols[as.integer(breaks)], pch = 15)

