###############################################################################

## The interpretation of feature selection results on real datasets in
## CDNFS: A New Filter Method for Nonlinear Feature Selection in High Dimensional Data

###############################################################################

source("function_CDNFS.R")

## real datasets for regression
for (iteration in 1:10){
  class_type <- "regression"
  simulation <- "slice_localization_data"
  

  path = "../../data/real_data/" 
  data <- read.csv(paste0(path, simulation, ".csv"), header = TRUE)
  
  path = paste0("../", simulation, "/repeat_", i_repeat, "/")
  data_index <- read.csv(paste0(path, "realdata_regression_data_sampled_2675_slice_localization_data.csv"),
                         header = TRUE)
  
  
  outputDir <- paste0("../visualization/")
  repeatDir <- simulation
  if (!file.exists(outputDir)){
    dir.create(file.path(outputDir))
  }
  if (!file.exists(paste0(outputDir, repeatDir))){
    dir.create(file.path(outputDir, repeatDir))
  }
  setwd(paste0(outputDir, repeatDir))
  
  
  index <- data_index[,1]
  data_x <- data.frame(data[index, -c(1, ncol(data))])
  colnames(data_x) <- c(paste0("X", 1:ncol(data_x)))
  data_y <- data.frame(y = data[index, ncol(data)])
  
  if (any(is.na(as.numeric(data_y$y)))){
    unique_y = unique(data_y$y)
    df <- data.frame(y = unique_y, result = 1:length(unique_y))
    Leaders <- data.frame(y = data_y$y, ID = seq(1:nrow(data_y)))
    new_y <- merge(df, Leaders, by = "y", all = TRUE)
    new_y <- new_y[order(match(new_y$ID, Leaders$ID)),]
    data_y <- data.frame(y = new_y$result)
  }
  data <- cbind(data_y, data_x)
  
  result <- CDNFS_visual(data, class_col = 1, d = NULL)
  col_remain <- result$feature$col_remain
  model = "CDNFS"
  print(col_remain)
  filename = paste0(simulation,"_variables_selected_", model, "_", iteration, ".csv")
  write.csv(result$feature, filename, row.names = FALSE)
  
  setwd("../../experiments")
}

## real datasets for classification
for (iteration in 1:10){
  class_type <- "classification"
  simulation <- c("warpAR10P", "Isolet1", "ORL")
  
  for (i_simu in 1:length(simulation)){
    
    path = "../../data/real_data/" 
    data <- read.csv(paste0(path, simulation[i_simu], ".csv"), header = TRUE)
    
    
    outputDir <- paste0("../visualization/")
    repeatDir <- simulation
    if (!file.exists(outputDir)){
      dir.create(file.path(outputDir))
    }
    if (!file.exists(paste0(outputDir, repeatDir))){
      dir.create(file.path(outputDir, repeatDir))
    }
    setwd(paste0(outputDir, repeatDir))
    
    
    data_x <- data.frame(data[, -ncol(data)])
    colnames(data_x) <- c(paste0("X", 1:ncol(data_x)))
    data_y <- data.frame(y = data[, ncol(data)])
    
    if (any(is.na(as.numeric(data_y$y)))){
      unique_y = unique(data_y$y)
      df <- data.frame(y = unique_y, result = 1:length(unique_y))
      Leaders <- data.frame(y = data_y$y, ID = seq(1:nrow(data_y)))
      new_y <- merge(df, Leaders, by = "y", all = TRUE)
      new_y <- new_y[order(match(new_y$ID, Leaders$ID)),]
      data_y <- data.frame(y = new_y$result)
    }
    data <- cbind(data_y, data_x)
    
    result <- CDNFS_visual(data, class_col = 1, d = NULL)
    col_remain <- result$feature$col_remain
    model = "CDNFS"
    print(col_remain)
    filename = paste0(simulation[i_simu],"_variables_selected_", model, "_", iteration, ".csv")
    write.csv(result$feature, filename, row.names = FALSE)
    
    setwd("../../experiments")
  }
}
