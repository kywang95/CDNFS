###############################################################################

## Evaluate the feature selection performance on simulations in
## CDNFS: A New Filter Method for Nonlinear Feature Selection in High Dimensional Data

###############################################################################


cases <- "classification_simulation" # "regression_simulation"

x.real = c(1:5)
p = 1000
g = 10
s = 50

repeat_times <- 10
start_times <- 1
iteration <- c(start_times:repeat_times)
method <- c("CDNFS", "FOCI", "DC_SIS", "DCOL", "mRMR", "NMIFS")


prob_evaluation <- function(x.model, x.real = c(1:5), p = 1000, g = 10, s = 50){
  n = length(x.real)
  p.g <- p/g
  g.non0 <- g-1
  p.model <- rep(0, times = n+1)
  x.selected <- x.model[1:min(length(x.model), s)]
  if (any(x.selected>=p.g*g.non0))
    x.selected <- x.selected[-which(x.selected >= p.g*g.non0)]
  x.selected <- x.selected%%p.g
  
  k = 0
  m = rep(FALSE, times = n)
  for (i in 1:n){
    judge.p <- (x.real[i] %in% x.selected)
    k = k + judge.p
    p.model[i] <- p.model[i] + judge.p
  }
  if (k == n){
    p.model[n+1] <- 1
  }
  
  return(p.model)
}

min.model_evaluation <- function(x.model, x.real = c(1:5), 
                                 p = 1000, g = 10){
  p.g <- p/g
  g.non0 <- g-1
  x.selected <- x.model
  if (any(x.selected>=p.g*g.non0))
    x.selected[which(x.selected >= p.g*g.non0)] = 0
  x.selected <- x.selected%%p.g
  x.subset <- x.real
  if (!is.null(length(x.selected))){
    for (i in 1:length(x.selected)){
      if (x.selected[i] %in% x.subset){
        min.x.model <- i
        x.subset <- x.subset[-which(x.subset == x.selected[i])]
      }
    }
  }
  if (length(x.subset)) min.x.model <- NA
  
  return(min.x.model)
}

result <- list()
for (i_model in 1:length(method)){
  performance <- data.frame(t(c(x.real, 1, 1, 0)))
  for (iter in iteration){
    data <- read.csv(paste0("../FS_performance/", cases, "/",
                            cases, "_variables_selected_", method[i_model], "_", iter, ".csv"))
    feature <- data[,1]
    time <- read.csv(paste0("../FS_performance/", cases, "/",
                            cases, "_variable_selection_time_", method[i_model], "_", iter, ".csv"))
    if (method[i_model]=="DCOL"){
      feature <- feature[feature!=0]
    }
    # Evaluation
    probX <- prob_evaluation(feature, x.real, p, g, s=20)
    minX <- min.model_evaluation(feature, x.real, p, g)
    performance <- rbind(performance, as.numeric(t(c(probX, minX, time))))
  }
  performance <- performance[-1,]
  colnames(performance) <- c(x.real, "All", "min", "time")
  result[[i_model]] <- performance
}

result

p.model <- matrix(0, nrow = length(method), ncol = length(x.real)+1)
min.model <- matrix(0, nrow = length(method), ncol = 6)

for (i_med in 1:length(method)){
  i_result <- result[[i_med]]
  p.model[i_med, ] <- colMeans(i_result[, 1:(length(x.real)+1)])
  min.model[i_med, 1:5] <- quantile(i_result[, ncol(i_result)-1], 
                                  probs = seq(0, 1, 0.25), na.rm = TRUE)
  min.model[i_med, 6] <- mean(i_result[, ncol(i_result)])
}

rownames(p.model) <- method
colnames(p.model) <- c(paste0("G", x.real), "All")
rownames(min.model) <- method
colnames(min.model) <- c(seq(0,1,0.25), "time")
write.csv(as.data.frame(p.model), paste0("../evaluation_p_", cases, ".csv"))
write.csv(as.data.frame(min.model), paste0("../evaluation_min_", cases, ".csv"))

