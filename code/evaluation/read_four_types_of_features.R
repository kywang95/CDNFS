###############################################################################

## Evaluate the feature selection performance on the real datasets in
## CDNFS: A New Filter Method for Nonlinear Feature Selection in High Dimensional Data

###############################################################################

library(stringr)
repetition <- 10
simulation <- c("warpAR10P", "Isolet1", "ORL", "slice_localization_data")


for (i_simu in 1:length(simulation)){
  path = "../../data/real_data/"
  data <- read.csv(paste0(path, simulation[i_simu], ".csv"), header = TRUE)
  if (!i_simu==length(simulation)){
    data_x <- data.frame(data[, -ncol(data)])
  }else{
    data_x <- data.frame(data[, -c(1, ncol(data))])
  }
  colnames(data_x) <- c(paste0("X", 1:ncol(data_x)))
  data_y <- data.frame(y = data[, ncol(data)])
  data <- data.frame(data_x, y = data_y)
  
  
  ##############################################################################
  obs1 <- data[1, ]
  obs1.x <- as.numeric(obs1[1:(length(obs1)-1)])
  obs1.y <- obs1[length(obs1)]
  
  result <- matrix(0, nrow = repetition, ncol = 4)
  for (iteration in 1:repetition){
    nonfeature <- read.csv(paste0("../visualization/", simulation[i_simu], "/", 
                                  simulation[i_simu], "_delete_lack_info_",
                                  iteration, ".csv"), header = TRUE)
    nonsignal <- as.numeric(str_remove(string = nonfeature$zero, pattern = "V"))
    
    CODECfeature <- read.csv(paste0("../visualization/", simulation[i_simu], "/", 
                                    simulation[i_simu], "_variables_selected_codec_GS_",
                                    iteration, ".csv"), header = TRUE)
    nfeature <- nrow(CODECfeature)
    relevant <- as.numeric(str_remove(string = CODECfeature[, 1], pattern = "V"))
    
    redundant <- NULL
    for (ifeature in 1:nfeature){
      GSfeature <- read.csv(paste0("../visualization/", simulation[i_simu], "/", 
                                   simulation[i_simu], "_variables_deleted_GS_", iteration, 
                                   "_", ifeature, ".csv"), header = TRUE)
      if (nrow(GSfeature)){
        redundant <- c(redundant, 
                       as.numeric(str_remove(string = GSfeature[, 1], 
                                             pattern = "V")))
      }
    }
    
    result[iteration, ] <- c(length(nonsignal), 
                             length(relevant), 
                             length(redundant),
                             length(obs1.x)-length(nonsignal)-length(relevant)-length(redundant)
    )
  }
  
  sd.result <- apply(result, 2, sd)
  result <- rbind(result, colMeans(result))
  result <- rbind(result, sd.result)
  colnames(result) <- c("nonsignal", "relevant", "redundant", "remanent")
  rownames(result) <- c(paste0("rep_", c(1:repetition)), "mean", "sd")
  print(result)
  write.csv(result, paste0("../", simulation[i_simu], "_CDNFS_selected_variables.csv"))
}

