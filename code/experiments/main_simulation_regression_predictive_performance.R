###############################################################################

## The predictive result on regression simulations in
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
  
  
  outputDir <- paste0("../", simulation)
  repeatDir <- paste0("repeat_", iteration)
  if (!file.exists(outputDir)){
    dir.create(file.path(outputDir))
  }
  if (!file.exists(file.path(outputDir, repeatDir))){
    dir.create(file.path(outputDir, repeatDir))
  }
  setwd(file.path(outputDir, repeatDir))
  
  
  # 5-fold cross-validation
  library(plyr)
  CVgroup <- function(k, datasize){
    cvlist <- list()
    # set.seed(123)
    n <- rep(1:k, ceiling(datasize/k))[1:datasize]    
    temp <- sample(n,datasize)   
    x <- 1:k
    dataseq <- 1:datasize
    cvlist <- lapply(x,function(x) dataseq[temp==x])  
    return(cvlist)
  }
  
  k <- 5
  datasize <- nrow(data)
  cvlist <- CVgroup(k = k, datasize = datasize)
  
  for (i_fold in c(1:k)){
    # i_fold = 1
    train <- data[-cvlist[[i_fold]],]
    test <- data[cvlist[[i_fold]],]
    write.csv(test, paste0(simulation,"_data_test_", iteration, "_", i_fold, ".csv"))
    
    train_x <- data.frame(train[,-1])
    test_x <- data.frame(test[,-1])
    train_y <- data.frame(y = train[, 1])
    test_y <- data.frame(y = test[, 1])
    train <- cbind(train_y, train_x)
    test <- cbind(test_y, test_x)
    
    
    x_index <- colnames(train[,-1])
    f <- as.formula(paste("y ~", paste(x_index, collapse = " + ")))
    model = "full"
    # SVM
    output_SVM_full <- SVM_regression(train, test, f, scaled=NULL, model, simulation, iteration, i_fold)
    # CART
    output_CART_full <- CART_regression(train, test, f, scaled=NULL, model, simulation, iteration, i_fold)
    # lightGBM
    output_lightGBM_full <- lightGBM_regression(train, test, scaled=NULL, model, simulation, iteration, i_fold)
    # 2-layer MLP
    output_NN2_full <- NN_regression(train, test, f, active_f, num.hidden.node, scaled=NULL, model, simulation, iteration, i_fold)
    
    
    ##########################################################################
    ## Feature selection
    d.DC <- 50
    
    for (index_method in c(1:6)){
      train_y <- data.frame(y = train[,1])
      train_x <- train[,-1]
      test_y <- data.frame(y = test[,1])
      test_x <- test[,-1]
      if (index_method == 1){
        # my method
        # Reduce the dimension with Gram-Schmidt transformation. Train the parameters using only the training set;
        result <- CDNFS_main(train, class_col = 1, d = NULL)
        col_remain <- result$feature$col_remain
        model = "CDNFS"
        print(col_remain)
        d.DC <- length(col_remain)
      }else
        if (index_method == 2){
          # FOCI - A simple measure of conditional dependence
          if(all(apply(train_x,2,var)!=0)){
            result <- foci(as.numeric(train_y$y), train_x, stop = TRUE, numCores = 1, 
                           num_features = d.DC)
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
                                               stop.var.count = d.DC)
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
                result <- NMIFS(train_x, train_y, num_features = d.DC, method = "emp")
                col_remain <- result$feature$feature
                model = "NMIFS"
                print(col_remain)
              }else
                if (index_method == 6){
                  # Mutual information - Using Mutual Information for Selecting Features in Supervised Neural Net Learning
                  result <- mRMR(train_x, train_y, num_features = d.DC, method = "emp")
                  col_remain <- result$feature$feature
                  model = "mRMR"
                  print(col_remain)
                }
      
      # save the variable selection result of each method 
      filename = paste0(simulation,"_variables_selected_", model, "_", iteration, "_", i_fold, ".csv")
      write.csv(result$feature, filename, row.names = FALSE)
      filename = paste0(simulation,"_variable_selection_time_", model, "_", iteration, "_", i_fold, ".csv")
      col.time <- data.frame(time.sec = result$time)
      write.table(col.time, filename, row.names = FALSE)
      
      
      # Obtain the features extracted from the test set using the variable selection methods.
      train_x <- data.frame(train_x[, col_remain])
      train_VS <- cbind(train_y, train_x)
      colnames(train_VS) <- c("y", paste0("X", col_remain))
      test_x <- data.frame(test_x[, col_remain])
      test_VS <- cbind(test_y, test_x)
      colnames(test_VS) <- c("y", paste0("X", col_remain))
      
      # model construction
      x_index <- paste0("X", col_remain)
      f <- as.formula(paste("y ~", paste(x_index, collapse = " + ")))
      # SVM
      output_SVM_FS <- SVM_regression(train, test, f, scaled=NULL, model, simulation, iteration, i_fold)
      # CART
      output_CART_FS <- CART_regression(train, test, f, scaled=NULL, model, simulation, iteration, i_fold)
      # lightGBM
      output_lightGBM_FS <- lightGBM_regression(train_VS, test_VS, scaled=NULL, model, simulation, iteration, i_fold)
      # 2-layer MLP
      output_NN2_FS <- NN_regression(train, test, f, active_f, num.hidden.node, scaled=NULL, model, simulation, iteration, i_fold)


      output <- cbind(output_SVM_FS, output_CART_FS, output_lightGBM_FS, output_NN2_FS)
      rownames(output) <- c("insample_MSE", "outsample_MSE")
      colnames(output) <- paste0(c("SVM", "CART", "lightGBM", "NN2"), "_", rep(model, each = 4))
      write.csv(output, paste0("result_", simulation,"_", model,"_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
    }
  }
  setwd("../../experiments")
}


