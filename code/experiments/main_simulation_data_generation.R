###############################################################################

## The generation of the regression and classification simultations in
## CDNFS: A New Filter Method for Nonlinear Feature Selection in High Dimensional Data

###############################################################################


for (iteration in 1:10){
  setwd(paste0("../../simulation_data/"))
  
  ###############################################################################
  ## Regression Simulation 1
  ###############################################################################
  simulation = "regression_simulation"
  sigmoid <- function(x){1 / (1+exp(-x))}
  
  n = 5000
  p = 100
  g = 10
  z = matrix(rnorm(n*(p+g+1), mean = 0, sd = 1), nrow = n)
  eps = z[, (p+g+1)]
  e = z[, (p+g)]
  z0 = z[, 1:p]
  z1 = (e + z0) / 2
  x <- z1
  lambda = 1e-2
  for (idex in 1:(g-2)){
    x = cbind(x, z1+z[, (p+idex)]*lambda)
  }
  for (idex in 1:p){
    mask <- rep(0, times = n)
    mask[sample(n, n*0.001, replace = FALSE)] <- 1
    noise <- mask*rnorm(n, mean = 0, sd = 1e-3)
    x = cbind(x, noise)
  }
  colnames(x) <- paste0("X", 1:ncol(x))
  y = x[,1]*x[,2] + sigmoid(4*x[,3]) + cos(pi*x[,4]*x[,5]) + 0.1*eps
  data <- data.frame(y, x)
  colnames(data) <- c("y", paste0("X", 1:ncol(x)))
  write.csv(data, paste0(simulation,"_data_", iteration, ".csv"), row.names = FALSE)
  
  
  ###############################################################################
  ## Classification Simulation 1
  ###############################################################################
  simulation = "classification_simulation"
  n = 5000
  p = 100
  g = 10
  z = matrix(rnorm(n*(p+g+1), mean = 0, sd = 1), nrow = n)
  eps = z[, (p+g+1)]
  e = z[, (p+g)]
  z0 = z[, 1:p]
  z1 = (e + z0) / 2
  x <- z1
  lambda = 1e-2
  for (idex in 1:(g-2)){
    x = cbind(x, z1+z[, (p+idex)]*lambda)
  }
  
  for (idex in 1:p){
    mask <- rep(0, times = n)
    mask[sample(n, n*0.001, replace = FALSE)] <- 1
    noise <- mask*rnorm(n, mean = 0, sd = 1e-3)
    x = cbind(x, noise)
  }
  colnames(x) <- paste0("X", 1:ncol(x))

  y <- tanh(2*tanh(2*x[,1]-x[,2]))+2*tanh(tanh(x[,3]-2*x[,4])-tanh(2*x[,5])) + 0.1*eps 
  # plot(density(y))
  class <- ifelse(y>-2, 
                  ifelse(y>-1, 
                         ifelse(y>0, 
                                ifelse(y>1, 
                                       ifelse(y>2, 6, 5), 4), 3), 2), 1)
  
  data <- data.frame(y = class, x)
  colnames(data) <- c("y", paste0("X", 1:ncol(x)))
  data$y <- as.factor(data$y)
  write.csv(data, paste0(simulation,"_data_", iteration, ".csv"), row.names = FALSE)
  
}