###############################################################################

## Interpret the feature selection performance on ORL dataset in
## CDNFS: A New Filter Method for Nonlinear Feature Selection in High Dimensional Data

###############################################################################

simulation = "ORL"

path = "../../data/real_data/"
data <- read.csv(paste0(path, simulation, ".csv"), header = TRUE)

class <- c(1:40)
n = 10
img_row <- 32
obs <- data[rep((class-1)*10, each = n)+rep(1:n, times = length(class)), ]

for (i_class in class){
  for (i_obs in 1:n){
    
    ######################################################################################
    # images
    
    ######################################################################################
    
    obs1 <- obs[i_obs + n*(i_class-1), ]
    obs1.x <- as.numeric(obs1[1:(length(obs1)-1)])
    obs1.y <- obs1[length(obs1)]
    mat1.x <- matrix(obs1.x, nrow = img_row, byrow = FALSE)

    
    library (ggplot2)
    library (reshape2)
    data1 <- as.data.frame(mat1.x)
    colnames(data1) <- paste0("C", c(1:ncol(data1)))
    data1 <- cbind(Row = paste0("R", c(1:nrow(data1))), data1)
    data2 <- melt (data1, id="Row") # transform into long vector
    data2$variable= factor(data2$variable, levels=paste0("C", c(1:ncol(data1))))   # set the sequence  
    data2$Row= factor(data2$Row, levels=paste0("R", c(nrow(data1):1)))  # set the sequence 
    img <- ggplot(data2,aes(x=variable, y=Row, fill=value)) + # draw the heatmap
      geom_raster()+ scale_fill_gradient(low="black", high="white") + # fill in different colors
      theme(rect = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.title=element_blank(),
            legend.position="none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA), # remove the backgroup of the plot  
            panel.border = element_blank(), # remove the frame
      )
    
    # ggsave(filename = paste0(simulation, "_origin_class_", i_class, "_obs_", i_obs, ".png"),
    #        plot = img,
    #        units="in", dpi=300, width=4, height=4, device="png", bg = "transparent")
    
    ######################################################################################
    # features
    
    ######################################################################################
    library(stringr)
    path  = paste0("../visualization/", simulation, "/")
    iteration = 10 # 1,...,10
    nonfeature <- read.csv(paste0(path, simulation, "_delete_lack_info_",
                                  iteration, ".csv"), header = TRUE)
    nonsignal <- as.numeric(str_remove(string = nonfeature$zero, pattern = "V"))
    
    CODECfeature <- read.csv(paste0(path, simulation, "_variables_selected_codec_GS_",
                                    iteration, ".csv"), header = TRUE)
    nfeature <- nrow(CODECfeature)
    relevant <- as.numeric(str_remove(string = CODECfeature[, 1], pattern = "V"))
    
    redundant <- NULL
    for (ifeature in 1:nfeature){
      GSfeature <- read.csv(paste0(path, simulation, "_variables_deleted_GS_", iteration, 
                                   "_", ifeature, ".csv"), header = TRUE)
      if (nrow(GSfeature)){
        redundant <- c(redundant, 
                       as.numeric(str_remove(string = GSfeature[, 1], 
                                             pattern = "V")))
      }
    }
    
    mat.zero <- rep(0, times = length(obs1.x))
    mat.feature <- mat.zero
    if (length(nonsignal)){
      mat.feature[nonsignal] = -1
    }
    if (length(relevant)){
      mat.feature[relevant] = 1
    }
    if (length(redundant)){
      mat.feature[redundant] = 2
    }
    mat.feature <- matrix(mat.feature, nrow = img_row, byrow = FALSE)
    
    library (ggplot2)
    library (reshape2)
    data1 <- as.data.frame(mat.feature)
    colnames(data1) <- paste0("C", c(1:ncol(data1)))
    data1 <- cbind(Row = paste0("R", c(1:nrow(data1))), data1)
    data2 <- melt(data1, id="Row") # transform into long vector
    
    color.map <- function(cl) {
      if(cl=="1") "darkred" else if (cl=="-1") "darkorange" else if (cl=="2") "darkblue" else NA
    }
    feature.color <- unlist(lapply(data2$value, color.map))
    data2 <- cbind(data2, feature.color)
    colnames(data2) <- c(colnames(data2)[1:3],"color")
    
    data2$variable= factor(data2$variable, levels=paste0("C", c(1:ncol(data1))))   # set the sequence 
    data2$Row= factor(data2$Row, levels=paste0("R", c(nrow(data1):1)))   # set the sequence
    ht.feature <- ggplot(data2, aes(x=variable, y=Row)) + # draw the heatmap
      geom_tile(aes(fill=color), alpha = 1) +
      scale_fill_identity() + 
      theme(rect = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.title=element_blank(),
            legend.position="none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(), # delete the backgroup of the plot 
            panel.border = element_blank(), # remove the frame
      )
    ht.feature
    
    
    library(cowplot)
    img.combine <- ggdraw() +
      draw_plot(img, scale = 1) + 
      draw_plot(ht.feature)
    
    
    ggsave(filename = paste0("../", simulation, "_OCDFS_class_", i_class, "_obs_", i_obs,
                             "_iter_", iteration, ".png"), 
           plot = img.combine, 
           units="in", dpi=300, width=4, height=4, device="png", bg = "transparent")
  }
}

