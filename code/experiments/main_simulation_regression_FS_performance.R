###############################################################################

## The feature selection performance on regression simulations in
## CDNFS: A New Filter Method for Nonlinear Feature Selection in High Dimensional Data

###############################################################################

is.regre = TRUE
simulation = "regression_simulation"

num.hidden.node = 100  
tanh <- function(x) {(exp(x)-exp(-x))/(exp(x)+exp(-x))}
softmax = custom <- function(x) {log(1+exp(x))}
sigmoid = function(x){1 / 1+exp(-x)}
active_f = "logistic" # tanh", "logistic", sigmoid, softmax

source("function_CDNFS.R")

i_repeat = 1
for (iteration in 1:10){
  path = "../../data/simulation_data/"
  data <- read.csv(paste0(path, simulation, "_data_", iteration, ".csv"))
  
  
  outputDir <- paste0("../FS_performance")
  repeatDir <- simulation
  if (!file.exists(outputDir)){
    dir.create(file.path(outputDir))
  }
  if (!file.exists(file.path(outputDir, repeatDir))){
    dir.create(file.path(outputDir, repeatDir))
  }
  setwd(file.path(outputDir, repeatDir))
  
  
  data_x <- data.frame(data[,-1])
  data_y <- data.frame(y = data[, 1])
  train <- cbind(data_y, data_x)
  
  ##########################################################################
  ## Feature selection
  d.DC <- 1000
  d.MI <- 50
  
  for (index_method in c(1:6)){
    train_y <- data.frame(y = train[,1])
    train_x <- train[,-1]
    if (index_method == 1){
      # my method
      # Reduce the dimension with Gram-Schmidt transformation. Train the parameters using only the training set;
      result <- CDNFS_main(train, class_col = 1, d = NULL)
      col_remain <- result$feature$col_remain
      model = "CDNFS"
      print(col_remain)
    }else
      if (index_method == 2){
        # FOCI - A simple measure of conditional dependence
        if(all(apply(train_x,2,var)!=0)){
          result <- foci(as.numeric(train_y$y), train_x, stop = TRUE, numCores = 1)
          col_remain <- result$feature$col_remain
          model = "FOCI"
          print(col_remain)
        } else {
          next
        }
      }else
        if (index_method == 3){
          # DC-SIS - Feature Screening via Distance Correlation Learning
          # d.DC <- ceiling(ncol(train_x)/log(ncol(train_x)))
          result <- DC_SIS(train_y$y, train_x, d = d.DC)
          col_remain <- result$feature$col_remain
          model = "DC_SIS"
          print(col_remain)
        }else
          if (index_method == 4){
            # DCOL - nonlinear variable selection with continuous outcome: A fully nonparametric incremental forward stagewise approach
            if(all(apply(train_x,2,var)!=0)){
              core_time <- Sys.time()
              library(nlnet)
              result <- nlnet::stage.forward(X = t(as.matrix(train_x)), 
                                             y = as.numeric(train_y$y),
                                             stop.alpha = 0.01,
                                             stop.var.count = ncol(train_x))
              col_remain <- result$found.pred
              model = "DCOL"
              if (length(result$found.pred)<length(result$ssx.rec)){
                result$found.pred <- c(result$found.pred, 
                                       rep(0, times = length(result$ssx.rec)-length(result$found.pred)))
              }
              feature <- as.data.frame(result)
              core_time <- difftime(Sys.time(), core_time, units = 'sec')
              result <- list(feature = feature, time = core_time)
              print(col_remain)
            } else {
              next
            }
          }else
            if (index_method == 5){
              # Mutual information - Using Mutual Information for Selecting Features in Supervised Neural Net Learning
              result <- NMIFS(train_x, train_y, num_features = d.MI, method = "emp")
              col_remain <- result$feature$feature
              model = "NMIFS"
              print(col_remain)
            }else
              if (index_method == 6){
                # Mutual information - Using Mutual Information for Selecting Features in Supervised Neural Net Learning
                result <- mRMR(train_x, train_y, num_features = d.MI, method = "emp")
                col_remain <- result$feature$feature
                model = "mRMR"
                print(col_remain)
              }
    
    # save the variable selection result of each method 
    filename = paste0(simulation,"_variables_selected_", model, "_", iteration, ".csv")
    write.csv(result$feature, filename, row.names = FALSE)
    filename = paste0(simulation,"_variable_selection_time_", model, "_", iteration, ".csv")
    col.time <- data.frame(time.sec = result$time)
    write.table(col.time, filename, row.names = FALSE)
    
  }
  setwd("../../experiments/")
}


