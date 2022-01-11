###############################################################################

## Evaluate the predictive performance on the real datasets in
## CDNFS: A New Filter Method for Nonlinear Feature Selection in High Dimensional Data

###############################################################################

evaluation <- c("MSE", "ACC")
cases <- c("regression_simulation", "classification_simulation")

method <- c("CDNFS", "FOCI", "DC_SIS", "DCOL", "mRMR", "NMIFS")
start_times <- 1
repeat_times <- 10
k <- 5 # k-fold analysis
iteration <- c(1:k)

outsample_result_svm <- matrix(0, nrow = k*repeat_times+2, ncol = length(method))
outsample_result_CART <- matrix(0, nrow = k*repeat_times+2, ncol = length(method))
outsample_result_lightGBM <- matrix(0, nrow = k*repeat_times+2, ncol = length(method))
outsample_result_NN2 <- matrix(0, nrow = k*repeat_times+2, ncol = length(method))
for (iter in 1:length(cases)){
  for (i_repeat in start_times:repeat_times){
    pathname <- paste0("../", cases[iter], "/repeat_", i_repeat, "/")
    for (i_fold in 1:k){
      for (imethod in 1:length(method)){
        filename <- paste0(pathname, "result_", cases[iter], "_", method[imethod], 
                           "_", i_repeat, "_", i_fold, ".csv")
        fullname <- paste(pathname, "result_", cases[iter], "_full", 
                          "_", i_repeat, "_", i_fold, ".csv")
        if (file.exists(filename)){
          data <- read.csv(filename)
        } else{
          data <- read.csv(fullname)
        }
        outsample_result_svm[i_fold+(i_repeat-1)*k, imethod] <- data[2,2]
        outsample_result_CART[i_fold+(i_repeat-1)*k, imethod] <- data[2,3]
        outsample_result_lightGBM[i_fold+(i_repeat-1)*k, imethod] <- data[2,4]
        outsample_result_NN2[i_fold+(i_repeat-1)*k, imethod] <- data[2,5]
      }
    }
  }
  outsample_result_svm[(k*repeat_times)+1,] <- colMeans(outsample_result_svm[(k*(start_times-1)+1):(k*repeat_times),])
  outsample_result_CART[(k*repeat_times)+1,] <- colMeans(outsample_result_CART[(k*(start_times-1)+1):(k*repeat_times),])
  outsample_result_lightGBM[(k*repeat_times)+1,] <- colMeans(outsample_result_lightGBM[(k*(start_times-1)+1):(k*repeat_times),])
  outsample_result_NN2[(k*repeat_times)+1,] <- colMeans(outsample_result_NN2[(k*(start_times-1)+1):(k*repeat_times),])
  
  outsample_result_svm[(k*repeat_times)+2,] <- apply(outsample_result_svm[(k*(start_times-1)+1):(k*repeat_times),], 2, sd)
  outsample_result_CART[(k*repeat_times)+2,] <- apply(outsample_result_CART[(k*(start_times-1)+1):(k*repeat_times),], 2, sd)
  outsample_result_lightGBM[(k*repeat_times)+2,] <- apply(outsample_result_lightGBM[(k*(start_times-1)+1):(k*repeat_times),], 2, sd)
  outsample_result_NN2[(k*repeat_times)+2,] <- apply(outsample_result_NN2[(k*(start_times-1)+1):(k*repeat_times),], 2, sd)
  
  
  outsample_result_CART <- data.frame(outsample_result_CART)
  colnames(outsample_result_CART) <- method
  rownames(outsample_result_CART) <- c(paste0(rep("outsample_", k*repeat_times),
                                              c(paste0(rep(1:repeat_times, each = k),"_", rep(1:k,times = repeat_times)))),
                                       "mean", "sd")
  
  outsample_result_svm <- data.frame(outsample_result_svm)
  colnames(outsample_result_svm) <- method
  rownames(outsample_result_svm) <- c(paste0(rep("outsample_", k*repeat_times),
                                             c(paste0(rep(1:repeat_times, each = k),"_", rep(1:k,times = repeat_times)))),
                                      "mean", "sd")
  
  outsample_result_lightGBM <- data.frame(outsample_result_lightGBM)
  colnames(outsample_result_lightGBM) <- method
  rownames(outsample_result_lightGBM) <- c(paste0(rep("outsample_", k*repeat_times),
                                                  c(paste0(rep(1:repeat_times, each = k),"_", rep(1:k,times = repeat_times)))),
                                           "mean", "sd")
  
  outsample_result_NN2 <- data.frame(outsample_result_NN2)
  colnames(outsample_result_NN2) <- method
  rownames(outsample_result_NN2) <- c(paste0(rep("outsample_", k*repeat_times),
                                             c(paste0(rep(1:repeat_times, each = k),"_", rep(1:k,times = repeat_times)))),
                                      "mean", "sd")
  
  write.csv(outsample_result_svm, paste0("../", evaluation[iter], "_outsample_SVM_", cases[iter], ".csv"))
  write.csv(outsample_result_CART, paste0("../", evaluation[iter], "_outsample_CART_", cases[iter], ".csv"))
  write.csv(outsample_result_lightGBM, paste0("../", evaluation[iter], "_outsample_lightGBM_", cases[iter], ".csv"))
  write.csv(outsample_result_NN2, paste0("../", evaluation[iter], "_outsample_MLP_", cases[iter], ".csv"))
  
}

