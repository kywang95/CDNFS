###############################################################################

## All of the procedure code necessary for the paper
## CDNFS: A New Filter Method for Nonlinear Feature Selection in High Dimensional Data

###############################################################################


###############################################################################

## element: CODEC
# Estimate the conditional dependence coefficient (CODEC)
# Azadkia, M. and Chatterjee, S. 
# A simple measure of conditional dependence
# Annuals of Statistics

# MAIN FUNCTION:
# codec:  Performs the conditional dependence coefficient (CODEC)
# HELPER FUNCTIONS for codec
# .estimateConditionalQ
# .estimateConditionalS
# .estimateConditionalT
# .estimateQ
# .estimateS
# .estimateT
# randomNNs

###############################################################################

codec <- function(Y, Z, X = NULL, na.rm = TRUE){
  if(is.null(X)) {
    # if inputs are not in proper matrix format change if possible
    # otherwise send error
    if(!is.vector(Y)) {
      Y = as.vector(Y)
    }
    if(!is.matrix(Z)) {
      Z = as.matrix(Z)
    }
    if((length(Y) != nrow(Z))) stop("Number of rows of Y and X should be equal.")
    if (na.rm == TRUE) {
      # NAs are removed here:
      ok = complete.cases(Y,Z)
      Z = as.matrix(Z[ok,])
      Y = Y[ok]
    }
    
    n = length(Y)
    if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
    
    q = ncol(Z)
    return(.estimateT(Y, Z))
  }
  # if inputs are not in proper matrix format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Z)) {
    Z = as.matrix(Z)
  }
  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  if((length(Y) != nrow(Z))) stop("Number of rows of Y and Z should be equal.")
  if((nrow(Z) != nrow(X))) stop("Number of rows of Z and X should be equal.")
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok = complete.cases(Y,Z,X)
    Z = as.matrix(Z[ok,])
    Y = Y[ok]
    X = as.matrix(X[ok,])
  }
  
  n = length(Y)
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
  
  q = ncol(Z)
  p = ncol(X)
  
  return(.estimateConditionalT(Y, Z, X))
}

# codec_S -------------------------------------------------------------------------
#' Estimate the conditional dependence coefficient (CODEC)
codec_S <- function(Y, X = NULL){
  if(is.null(X)) {
    # if inputs are not in proper matrix format change if possible
    # otherwise send error
    if(!is.vector(Y)) {
      Y = as.vector(Y)
    }
    
    n = length(Y)
    if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
    
    return(.estimateS(Y))
  }
  # if inputs are not in proper matrix format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  
  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  
  n = length(Y)
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
  
  return(.estimateConditionalS(Y, X))
}

codec_Q <- function(Y, Z, X = NULL, na.rm = TRUE){
  if(is.null(X)) {
    # if inputs are not in proper matrix format change if possible
    # otherwise send error
    if(!is.vector(Y)) {
      Y = as.vector(Y)
    }
    if(!is.matrix(Z)) {
      Z = as.matrix(Z)
    }
    if((length(Y) != nrow(Z))) stop("Number of rows of Y and X should be equal.")
    if (na.rm == TRUE) {
      # NAs are removed here:
      ok = complete.cases(Y,Z)
      Z = as.matrix(Z[ok,])
      Y = Y[ok]
    }
    
    n = length(Y)
    if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
    
    q = ncol(Z)
    return(.estimateQ(Y, Z))
  }
  # if inputs are not in proper matrix format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Z)) {
    Z = as.matrix(Z)
  }
  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  if((length(Y) != nrow(Z))) stop("Number of rows of Y and Z should be equal.")
  if((nrow(Z) != nrow(X))) stop("Number of rows of Z and X should be equal.")
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok = complete.cases(Y,Z,X)
    Z = as.matrix(Z[ok,])
    Y = Y[ok]
    X = as.matrix(X[ok,])
  }
  
  n = length(Y)
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
  
  q = ncol(Z)
  p = ncol(X)
  
  return(.estimateConditionalQ(Y, X, Z))
}


# .estimateConditionalQ -------------------------------------------------------------------------
# Estimate Q(Y, Z | X)
#
# Estimate Q(Y, Z | X), the numinator of the measure of conditional dependence of Y on Z given X
#
# @param X: Matrix of predictors (n by p)
# @param Z: Matrix of predictors (n by q)
# @param Y: Vector (length n)
#
# @return estimation \eqn{Q(Y, Z|X)}
.estimateConditionalQ <- function (Y, X, Z) {
  
  id <- group <- rnn <- NULL
  
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Z)) {
    Z = as.matrix(Z)
  }
  
  n = length(Y)
  
  W = cbind(X, Z)
  
  # compute the nearest neighbor of X
  nn_X = RANN::nn2(X, query = X, k = 3)
  nn_index_X = nn_X$nn.idx[, 2]
  # handling repeated data
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  
  nn_index_X[repeat_data] = df_X$rnn
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  # compute the nearest neighbor of W
  nn_W = RANN::nn2(W, query = W, k = 3)
  nn_index_W = nn_W$nn.idx[, 2]
  repeat_data = which(nn_W$nn.dists[, 2] == 0)
  
  df_W = data.table::data.table(id = repeat_data, group = nn_W$nn.idx[repeat_data])
  df_W[, rnn := .randomNN(id), by = "group"]
  
  nn_index_W[repeat_data] = df_W$rnn
  # nearest neighbors with ties
  ties = which(nn_W$nn.dists[, 2] == nn_W$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_W[ties] = sapply(ties, helper_ties)
  }
  
  # estimate Q
  R_Y = rank(Y, ties.method = "max")
  Q_n = sum(apply(rbind(R_Y, R_Y[nn_index_W]), 2, min),
            -apply(rbind(R_Y, R_Y[nn_index_X]), 2, min)) / (n^2)
  return(Q_n)
}




# .estimateConditionalS -------------------------------------------------------------------------
# Estimate S(Y, X)
#
# Estimate S(Y, X), the denuminator of the measure of dependence of Y on Z given X
#
# @param X: Matrix of predictors (n by p)
# @param Y: Vector (length n)
#
# @return estimation \eqn{S(Y, X)}
.estimateConditionalS <- function (Y, X){
  
  id <- group <- rnn <- NULL
  
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  n = length(Y)
  
  # compute the nearest neighbor of X
  nn_X = RANN::nn2(X, query = X, k = 3)
  nn_index_X = nn_X$nn.idx[, 2]
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  
  nn_index_X[repeat_data] = df_X$rnn
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  # estimate S
  R_Y = rank(Y, ties.method = "max")
  S_n = sum(R_Y - apply(rbind(R_Y, R_Y[nn_index_X]), 2, min)) / (n^2)
  
  return(S_n)
}


# estimateConditionalT -------------------------------------------------------------------------
# Estimate T(Y, Z | X)
#
# Estimate T(Y, Z | X), the measure of dependence of Y on Z given X
#
# @param Y: Vector (length n)
# @param Z: Matrix of predictors (n by q)
# @param X: Matrix of predictors (n by p)
#
# @return estimation of \eqn{T(Y, Z|X)}.
.estimateConditionalT <- function(Y, Z, X){
  S = .estimateConditionalS(Y, X)
  
  # happens only if Y is constant
  if (S == 0) {
    return(1)
  } else {
    return(.estimateConditionalQ(Y, X, Z) / S)
  }
}




# .estimateQ -------------------------------------------------------------------------
# Estimate Q(Y, X)
#
# Estimate Q(Y, X), the numinator of the measure of dependence of Y on X
#
# @param X: Matrix of predictors (n by p).
# @param Y: Vector (length n).
#
# @return estimation of \eqn{Q(Y, X)}.
.estimateQ <- function(Y, X) {
  
  id <- group <- rnn <- NULL
  
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  
  n = length(Y)
  nn_X = RANN::nn2(X, query = X, k = 3)
  # remove the first nearest neighbor for each x which is x itself in case of no repeat data
  # when there is repeated data this is wrong but for repeated data we find the nearest
  # neighbors separately.
  nn_index_X = nn_X$nn.idx[, 2]
  
  # find all data points that are not unique
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  
  # for the repeated data points, choose one of their identicals at random and set its index
  # as the index of the nearest neighbor
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  nn_index_X[repeat_data] = df_X$rnn
  
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }
    
    nn_index_X[ties] = sapply(ties, helper_ties)
  }
  
  
  # estimate Q
  R_Y = rank(Y, ties.method = "max")
  L_Y = rank(-Y, ties.method = "max")
  Q_n = sum(apply(rbind(R_Y, R_Y[nn_index_X]), 2, min) - (L_Y^2)/n) / (n^2)
  
  return(Q_n)
}



# .estimateS -------------------------------------------------------------------------
# Estimate S(Y)
#
# Estimate S(Y) , the denuminator of the measure of dependence of Y on X
#
# @param Y: Vector (length n).
# @return estimation of \eqn{S(Y)}.
.estimateS <- function (Y) {
  n = length(Y)
  L_Y = rank(-Y, ties.method = "max")
  S_n = gmp::asNumeric(sum(gmp::as.bigz(L_Y) * gmp::as.bigz(n - L_Y))) / (n^3)
  return(S_n)
}



# .estimateT -------------------------------------------------------------------------
# Estimate T(Y, X)
#
# Estimate T(Y, X), the measure of dependence of Y on X
#
# @param X: Matrix of predictors (n by p)
# @param Y: Vector (length n)
# @return estimation of \eqn{T(Y, X) = Q(Y, X) / S(Y)}.
.estimateT <- function(Y, X){
  
  S = .estimateS(Y)
  # happens only if Y is a constant vector.
  if (S == 0) {
    return(1)
  } else {
    return(.estimateQ(Y, X) / S)
  }
}


# randomNN -------------------------------------------------------------------------
# Find the random nearest neighbors.
#
# Gives the indices of nearest neighbors of points that are not unique in the data set.
# For each point we know that to which cluster of points it belongs (repeated points),
# we chose one of those indices that is not equal to the same index of our given point
# at random as its nearest neighbor.
#
# @param ids: a vector of ids of groups that each index is a member of.
#
# @return a vector of indices.
.randomNN <- function(ids) {
  m <- length(ids)
  
  x <- sample(x = (m - 1), m, replace = TRUE)
  x <- x + (x >= (1:m))
  
  return(ids[x])
}

###############################################################################

## Baseline: foci
# FOCI is a variable selection algorithm based on the measure of conditional dependence
# Azadkia, M. and Chatterjee, S. 
# A simple measure of conditional dependence
# Annuals of Statistics

# MAIN FUNCTIONS:
# foci_main: the core function for implementing the foci algorithm
# foci:  Performs feature ordering by conditional independence

###############################################################################

foci_main <- function(Y, X, num_features = NULL, stop = TRUE, numCores = parallel::detectCores()){
  
  namesX <- colnames(X)
  if (is.null(num_features)) num_features = dim(X)[2]
  n = length(Y)
  p = ncol(X)
  Q = rep(0, num_features)
  index_select = rep(0, num_features)
  # select the first variable
  if (is.null(dim(X))) {
    seq_Q = .estimateQ(Y, X)
  } else {
    estimateQFixedY <- function(id){
      return(.estimateQ(Y, X[, id]))
    }
    seq_Q = parallel::mclapply(seq(1, p), estimateQFixedY, mc.cores = numCores)
    seq_Q = unlist(seq_Q)
  }
  
  Q[1] = max(seq_Q)
  if (Q[1] <= 0 & stop == TRUE) return(NULL)
  index_max = min(which(seq_Q == Q[1]))
  index_select[1] = index_max
  
  print(namesX[index_select[1]])
  print(Q[1])
  
  count = 1
  
  # select rest of the variables
  while (count < num_features) {
    seq_Q = rep(0, p - count)
    # indices that have not been selected yet
    index_left = setdiff(seq(1, p), index_select[1:count])
    
    # find the next best feature
    estimateQFixedYandSubX <- function(id){
      return(.estimateQ(Y, cbind(X[, c(index_select[1:count], id)])))
    }
    
    if (length(index_left) == 1) {
      seq_Q = estimateQFixedYandSubX(index_left[1])
    } else {
      seq_Q = parallel::mclapply(index_left, estimateQFixedYandSubX, mc.cores = numCores)
      seq_Q = unlist(seq_Q)
    }
    Q[count + 1] = max(seq_Q)
    index_max = min(which(seq_Q == Q[count + 1]))
    
    if (Q[count + 1] <= Q[count] & stop == TRUE) break
    index_select[count + 1] = index_left[index_max]
    
    print(namesX[index_select[count + 1]])
    print(Q[count + 1])
    
    count = count + 1
  }
  
  selectedVar = data.table::data.table(index = index_select[1:count], names = namesX[index_select[1:count]])
  stepT = Q / .estimateS(Y)
  result = list(selectedVar = selectedVar, stepT = stepT[1:count])
  class(result) = "foci"
  return(result)
}


foci <- function(Y, X, num_features = NULL, stop = TRUE, na.rm = TRUE,
                 standardize = "scale", numCores = parallel::detectCores(),
                 parPlat = 'none', printIntermed = TRUE) {
  core_time <- Sys.time()
  if (is.null(colnames(X))) {
    colnames(X) <- paste0('V',1:ncol(X))
    warning('X lacked column names, has been assigned V1, V2,...')
  }
  namesX <- colnames(X)
  
  # if inputs are not in proper format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y <- as.vector(unlist(Y))
  }
  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  if (is.null(num_features)) num_features <- dim(X)[2]
  
  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok <- complete.cases(Y,X)
    X <- as.matrix(X[ok,])
    Y <- Y[ok]
  }
  
  
  n <- length(Y)
  
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")
  
  p = ncol(X)
  
  if (num_features > p) stop("Number of features should not be larger than maximum number of original features.")
  if ((floor(num_features) != num_features) || (num_features <= 0)) stop("Number of features should be a positive integer.")
  
  if (!is.numeric(Y)) stop("currently FOCI does not handle factor Y")
  
  if (standardize == "scale") {
    for (i in 1:p) {
      if(length(unique(X[, i])) > 1) {
        X[,i] <- (X[,i] - mean(X[,i])) / sd(X[,i])
      }else{
        stop(paste0("Column ", i, " of X is constant."))
      }
    }
  }
  
  if (standardize == "bounded") {
    for (i in 1:p) {
      if(length(unique(X[, i])) > 1) {
        X[,i] <- (X[,i] - min(X[,i])) / (max(X[, i]) - min(X[, i]))
      }else{
        stop(paste0("Column ", i, " of X is constant."))
      }
    }
  }
  
  if (parPlat[1] == 'none') {
    res <- foci_main(Y, X, num_features = num_features,
                     stop = stop, numCores = numCores)
    feature <- data.frame(col_remain = res$selectedVar$index, stepT = res$stepT)
    
    core_time <- difftime(Sys.time(), core_time, units = 'sec')
    results <- list(feature = feature, 
                    time = core_time)
    return(results)
  }
  
  nr <- nrow(X)
  permRowNums <- sample(1:nr,nr,replace=FALSE)
  X <- X[permRowNums,]
  Y <- Y[permRowNums]
  
  rowNums <- parallel::splitIndices(length(Y), numCores)
  selectFromChunk <- function(nodeNum) {
    myRows <- rowNums[[nodeNum]]
    sel <- foci_main(Y[myRows], X[myRows,], stop = stop,
                     numCores = numCores)$selectedVar$index
  }
  
  if(inherits(parPlat,'cluster')) {
    cls <- parPlat
  }else if(parPlat == 'locThreads') {
    # set up the cluster (in multicore case, it's virtual)
    cls <- parallel::makeCluster(numCores)
  } else stop('invalid parPlat')
  
  # worker nodes load library
  parallel::clusterEvalQ(cls, library(FOCI))
  # ship data to workers
  parallel::clusterExport(cls, c('Y', 'X', 'rowNums', 'selectFromChunk'),
                          envir = environment())
  # drop-in replacement for mclapply
  slc <- parallel::parLapply(cls, seq(1, length(cls)), selectFromChunk)
  
  if (printIntermed) print(slc)
  
  slc <- Reduce(union, slc)
  
  res <- foci_main(Y, X[, slc], num_features, stop)
  # must translate indices in reduced system to those of original
  newIdxs <- res$selectedVar$index
  origIdxs <- slc[newIdxs]
  res$selectedVar$index <- origIdxs
  
  res$stepT = res$stepT[1:num_features]
  parallel::stopCluster(cls)
  
  feature <- data.frame(col_remain = res$selectedVar$index, stepT = res$stepT)
  
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = feature, 
                  time = core_time)
  return(results)
}


###############################################################################

## The Gram-Schimidt function used to reduce the dimensions of the explanatory variables

###############################################################################

Gram_Schmidt_delete <- function(D_GS, x, alpha = 0.1){
  N = nrow(x)
  p = ncol(x)
  index_x <- array(0, dim = p)
  u <- x
  
  for (i in 1:p){
    for (j in 1:ncol(D_GS)){
      d <- t(u[,i]) %*% D_GS[,j] / (t(D_GS[,j]) %*% D_GS[,j])
      u[,i] <- u[,i]-c(d)*D_GS[,j]
    }
  }
  # non-supervised G-S
  t_value<- rep(0, times = p)
  for (i in 1:p){
    if (var(u[,i])==0 | var(x[,i])==0){
      t_value[i] = 0
    }else{
      r = cor(u[,i], x[,i])
      # t_value[i] = r*sqrt(N-2)/sqrt(1-r^2)
      t_value[i] = var(u[,i])
    }
  }
  
  Gk <- colnames(u)[which(t_value < alpha)]
  # print(Gk)
  # print(t_value[t_value < alpha])
  
  # t_bound = qt(1-alpha, df = N-2)
  # Gk <- colnames(u[, which(abs(t_value) < t_bound)])
  # print(Gk)
  # print(t_value[abs(t_value) < t_bound])
  
  return(Gk)
}


###############################################################################

## Baseline: Mutual information Feature Selection (MIFS)
# Using Mutual Information for Selecting Features in Supervised Neural Net Learning
# Roberto Battiti
# IEEE Transactions on Neural Networks

###############################################################################
MI_main <- function(X, Y, beta = 0.5, num_features = NULL, method = "emp"){
  core_time <- Sys.time()
  library(infotheo)
  X <- infotheo::discretize(X)
  Y <- infotheo::discretize(Y)
  if (is.null(num_features)) num_features = ncol(X)
  feature_set <- colnames(X)
  selected_set <- NULL
  MI <- rep(0, times = length(feature_set))
  for (i in 1:length(feature_set)){
    MI[i] <- mutinformation(X[,feature_set[i]], Y$y, method=method)
  }
  MI.Cf <- MI
  max_idex <- order(MI, decreasing = TRUE)[1]
  IC.S <- MI.Cf[max_idex]
  feature <- feature_set[max_idex]
  MI.Cf <- MI.Cf[-max_idex]
  feature_set <- feature_set[-max_idex]
  selected_set <- feature
  s = 1
  print(s)
  print(feature)
  print(IC.S)
  print(MI[max_idex])
  output_col <- data.frame(col_remain = feature, IC.S = IC.S, MI.max_idex = MI[max_idex])
  while(s < num_features){
    f = length(feature_set)
    MI.fs <- matrix(rep(0, times = f*s), nrow = f)
    for (i in 1:f){
      for (j in 1:s){
        MI.fs[i, j] <- mutinformation(X[,feature_set[i]], 
                                      X[,selected_set[j]], method=method)
      }
    }
    MI.f <- rowSums(MI.fs)
    MI <- MI.Cf-beta*MI.f
    max_idex <- order(MI, decreasing = TRUE)[1]
    IC.S <- MI.Cf[max_idex]
    feature <- feature_set[max_idex]
    MI.Cf <- MI.Cf[-max_idex]
    feature_set <- feature_set[-max_idex]
    selected_set <- c(selected_set, feature)
    s = length(selected_set)
    
    print(s)
    print(feature)
    print(IC.S)
    print(MI[max_idex])
    output_col <- rbind(output_col, c(feature, IC.S, MI[max_idex]))
  }
  colnames(output_col) <- c("feature", "IC.S", "MIFS.max_idex")
  
  library(stringr)
  output_col$feature <- as.numeric(str_remove(string = output_col$feature, pattern = "X"))
  
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = output_col, 
                  time = core_time)
  return(results)
}


###############################################################################

## Baseline: Normalized Mutual information Feature Selection (NMIFS)
# Normalized Mutual Information Feature Selection
# Pablo A. Estévez, Michel Tesmer, Claudio A. Perez, and Jacek M. Zurada
# IEEE Transactions on Neural Networks

###############################################################################
NMIFS <- function(X, Y, num_features = NULL, method = "emp"){
  core_time <- Sys.time()
  library(infotheo)
  X <- infotheo::discretize(X)
  Y <- infotheo::discretize(Y)
  if (is.null(num_features)) num_features = ncol(X)
  feature_set <- colnames(X)
  selected_set <- NULL
  MI <- rep(0, times = length(feature_set))
  H <- rep(0, times = length(feature_set))
  for (i in 1:length(feature_set)){
    MI[i] <- mutinformation(X[,feature_set[i]], Y$y, method=method)
    H[i] <- entropy(X[,feature_set[i]], method = method)
  }
  MI.Cf <- MI
  max_idex <- order(MI, decreasing = TRUE)[1]
  IC.S <- MI.Cf[max_idex]
  feature <- feature_set[max_idex]
  HS <- H[max_idex]
  HF <- H[-max_idex]
  MI.Cf <- MI.Cf[-max_idex]
  feature_set <- feature_set[-max_idex]
  selected_set <- feature
  s = 1
  print(s)
  print(feature)
  print(IC.S)
  print(MI[max_idex])
  output_col <- data.frame(col_remain = feature, IC.S = IC.S, MI.max_idex = MI[max_idex])
  while(s < num_features){
    f = length(feature_set)
    MI.fs <- matrix(rep(0, times = f*s), nrow = f)
    for (i in 1:f){
      for (j in 1:s){
        denominator <- min(HF[i], HS[j])
        MI.fs[i, j] <- mutinformation(X[,feature_set[i]], 
                                      X[,selected_set[j]], 
                                      method=method) / denominator
      }
    }
    MI.f <- rowSums(MI.fs)
    beta <- 1/s
    MI <- MI.Cf-beta*MI.f
    max_idex <- order(MI, decreasing = TRUE)[1]
    IC.S <- MI.Cf[max_idex]
    feature <- feature_set[max_idex]
    HS <- c(HS, HF[max_idex])
    HF <- HF[-max_idex]
    MI.Cf <- MI.Cf[-max_idex]
    feature_set <- feature_set[-max_idex]
    selected_set <- c(selected_set, feature)
    s = length(selected_set)
    
    print(s)
    print(feature)
    print(IC.S)
    print(MI[max_idex])
    output_col <- rbind(output_col, c(feature, IC.S, MI[max_idex]))
  }
  colnames(output_col) <- c("feature", "IC.S", "NMIFS.max_idex")
  
  library(stringr)
  output_col$feature <- as.numeric(str_remove(string = output_col$feature, pattern = "X"))
  
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = output_col, 
                  time = core_time)
  return(results)
}


###############################################################################

## Baseline: minimal-redundancy-maximal-relevance criterion (mRMR)
# Feature Selection Based on Mutual Information: Criteria of Max-Dependency, 
# Max-Relevance, and Min-Redundancy
# Hanchuan Peng, Fuhui Long, and Chris Ding
# IEEE Transactions on Neural Networks

###############################################################################
mRMR <- function(X, Y, num_features = NULL, method = "emp"){
  core_time <- Sys.time()
  library(infotheo)
  X <- infotheo::discretize(X)
  Y <- infotheo::discretize(Y)
  if (is.null(num_features)) num_features = ncol(X)
  feature_set <- colnames(X)
  selected_set <- NULL
  MI <- rep(0, times = length(feature_set))
  for (i in 1:length(feature_set)){
    MI[i] <- mutinformation(X[,feature_set[i]], Y$y, method=method)
  }
  MI.Cf <- MI
  max_idex <- order(MI, decreasing = TRUE)[1]
  IC.S <- MI.Cf[max_idex]
  feature <- feature_set[max_idex]
  MI.Cf <- MI.Cf[-max_idex]
  feature_set <- feature_set[-max_idex]
  selected_set <- feature
  s = 1
  print(s)
  print(feature)
  print(IC.S)
  print(MI[max_idex])
  output_col <- data.frame(col_remain = feature, IC.S = IC.S, MI.max_idex = MI[max_idex])
  while(s < num_features){
    f = length(feature_set)
    MI.fs <- matrix(rep(0, times = f*s), nrow = f)
    for (i in 1:f){
      for (j in 1:s){
        MI.fs[i, j] <- mutinformation(X[,feature_set[i]], 
                                      X[,selected_set[j]], method=method)
      }
    }
    MI.f <- rowSums(MI.fs)
    beta <- 1/s
    MI <- MI.Cf-beta*MI.f
    max_idex <- order(MI, decreasing = TRUE)[1]
    IC.S <- MI.Cf[max_idex]
    feature <- feature_set[max_idex]
    MI.Cf <- MI.Cf[-max_idex]
    feature_set <- feature_set[-max_idex]
    selected_set <- c(selected_set, feature)
    s = length(selected_set)
    
    print(s)
    print(feature)
    print(IC.S)
    print(MI[max_idex])
    output_col <- rbind(output_col, c(feature, IC.S, MI[max_idex]))
  }
  colnames(output_col) <- c("feature", "IC.S", "mRMR.max_idex")
  
  library(stringr)
  output_col$feature <- as.numeric(str_remove(string = output_col$feature, pattern = "X"))
  
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = output_col, 
                  time = core_time)
  return(results)
}


###############################################################################

## Baseline: DC-SIS
# Feature Screening via Distance Correlation Learning
# Runze Li, Wei Zhong, and Liping Zhu
# Journal of the American Statistical Association

###############################################################################
DC_SIS <- function(Y, X, d = NULL) {
  core_time <- Sys.time()
  if (is.null(dim(Y))){
    if (nrow(X) != length(Y) ) {
      stop("X and Y should have same number of rows!")
    }
  }else
    if (nrow(X) != nrow(Y) ) {
      stop("X and Y should have same number of rows!")
    }
  results <- data.frame()
  
  # The distance correlation measure is from the energy package by Rizzo and Szekely. 
  dcory<-function(Xk){energy::dcor(Xk,Y)}  
  v<-abs(apply(X,2,dcory))
  # rank order of estimated association strengths;
  Dword <- order(v,decreasing = T) 
  rank <- v[Dword]
  # rank<-match(v[Dword], v)
  if (is.null(d)){
    # d <- min(round(nrow(X)/log(nrow(X))), round(ncol(X)/log(nrow(X))))
    d <- ncol(X)
  }
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = data.frame(col_remain = Dword[1:d], measurement = rank[1:d]), 
                  time = core_time)
  return(results)
}


###############################################################################

## Baseline: ReliefF when classification, RReliefF when regression

###############################################################################
reliefF <- function(X, Y, y_class = "regression", num_features = NULL){
  core_time <- Sys.time()
  if (y_class == "regression")
    estimator = "RReliefFexpRank"
  if (y_class == "classification")
    estimator = "ReliefFexpRank"
  library(CORElearn)
  data <- data.frame(y = Y, X)
  x_index <- colnames(data[,-1])
  f <- as.formula(paste("y ~", paste(x_index, collapse = " + ")))
  estReliefF <- attrEval(f, data,
                         estimator= estimator # "RReliefFequalK", "RReliefFexpRank", "RReliefdistance"
                         # , ReliefIterations=30
  )
  if (is.null(num_features))
    num_features = length(estReliefF)
  col_remain <- order(estReliefF, decreasing = TRUE)[1:num_features]
  measurement <- estReliefF[col_remain]
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = data.frame(col_remain = col_remain, measurement = measurement), 
                  time = core_time)
  return(results)
}


###############################################################################

## The Main function of the proposed CDNFS

###############################################################################
CDNFS_main <- function(data, class_col = 1, d = NULL){
  core_time <- Sys.time()
  # reduce the dimension with Gram-Schmidt transformation. Train the parameters 
  # using only the training set
  y <- data.frame(class = data[, class_col])
  x = data.frame(data[,-class_col])
  name_x <- colnames(x)
  y$class <- as.numeric(y$class)
  p = ncol(x)
  colnames(x) <- paste0("V", c(1:p))
  print(dim(x))
  
  # delete the variables without information (sd(x) < 1e-5)
  library(infotheo)
  non_zero <- which(apply(infotheo::discretize(x), 2, entropy) > 1e-2)
  #  non_zero <- which(apply(x, 2, function(x) sd(x) > 1e-2))
  zero <- colnames(x[,-non_zero])
  delete_variables <- zero
  name_x <- colnames(x)[non_zero]
  x <- data.frame(x[, non_zero])
  colnames(x) <- name_x
  # print(zero)
  print(dim(x))
  
  x <- scale(x)
  N <- nrow(x)
  p = ncol(x)
  z = NULL
  stepT = c()
  
  if (is.null(d)) d.stop <- p
  else d.stop <- d
  
  while (p>0 && length(z)<d.stop){
    cor_yx <- array(0, dim = p)
    cor_S <- codec_S(y[,1], z)
    if (cor_S <= 0) break
    for (i in 1:p){
      cor_yx[i] <- codec_Q(y[,1], x[,i], z)/cor_S
    }
    k = order(cor_yx, decreasing = TRUE)[1]
    t = max(cor_yx)
    stepT = c(stepT, t)
    D = data.frame(x[,k])
    colnames(D) <- colnames(x)[k]
    print(c(round(t, 3), colnames(D)))
    if (is.null(z)){
      z = data.frame(D)
      
      D_GS <- z
    }else{
      name_z <- colnames(z)
      z = cbind(z, D)
      colnames(z) <- c(name_z, colnames(D))
      
      D_GS <- cbind(D_GS, D)
      colnames(D_GS) <- colnames(z)
      for (i in 1:(ncol(D_GS)-1)){
        d_GS <- t(D) %*% D_GS[,i] / (t(D_GS[,i]) %*% D_GS[,i])
        D_GS[, ncol(D_GS)] <- D_GS[, ncol(D_GS)] - c(d_GS) * D_GS[, i]
      }
    }
    if (t<=0 && is.null(d)) break
    if(ncol(x)>2){
      x <- x[, -k]
    }else
      if (ncol(x)==2){
        xnames <- colnames(x)
        x <- data.frame(x = x[,-k])
        colnames(x) <- xnames[-k]
      }else break
    
    # delete the variables which is highly correlated with the selected x[,k]
    # using Gram-Schimidt Orthogonalization
    Gk <- Gram_Schmidt_delete(D_GS, x, alpha = 0.05)
    # print(length(Gk))
    if (length(Gk)>0){
      # x[, Gk] <- list(NULL)
      x <- x[, (!colnames(x) %in% Gk), drop = FALSE]
    }
    delete_variables <- c(delete_variables, Gk)
    p = ncol(x)
    
    # print(ncol(z))
    print(dim(x))
    
    # filename = paste0(simulation,"_variables_selected_CDNFS.csv")
    # write.csv(data.frame(col_remain = colnames(z), stepT = stepT), filename, row.names = FALSE)
  }
  
  # Save the parameters of the Gram-Schimidt Orthogonalization transformation.
  col_remain <- colnames(z)
  library(stringr)
  col_remain <- as.numeric(str_remove(string = col_remain, pattern = "V"))
  output_col <- data.frame(col_remain, stepT = stepT)
  
  # print(col_remain)
  # print(delete_variables)
  
  if (is.null(d))
    output_col <- output_col[-length(col_remain),]
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = output_col, 
                  time = core_time)
  return(results)
}


plot_neural_net <- function(net = model.sol, model = "full", iteration, min_rep) {
  png(paste0('net_', model, '_', iteration, '.png'))
  library(NeuralNetTools)
  plot(net, rep = min_rep)
  dev.off()
}
NN_regression <- function(train, test, f, active_f, num.hidden.node, scaled=NULL, model, simulation, iteration, i_fold=1){
  ##########################################################################
  ## MLP
  library(neuralnet)
  library(NeuralNetTools)
  is.regre <- TRUE
  err.fct <- "sse"
  
  model.sol <- neuralnet(f, data = train, 
                         hidden = num.hidden.node, 
                         stepmax = 1e+05,
                         act.fct = active_f, # "tanh", "logistic", sigmoid, softmax
                         algorithm = "rprop+", # "rprop-", "sag", "slr", OR "backprop", learningrate=0.01,
                         threshold=0.05,
                         #                        learningrate.factor = list(minus = 0.1, plus = 1.0),
                         rep = 3,                   # the training repetitions for the NN
                         err.fct = err.fct,
                         lifesign = "full",      # "minimal", "full" # verbose output during training
                         lifesign.step = 5000,
                         linear.output = is.regre, 
                         likelihood = TRUE,
                         still.times = 1e+04)
  
  zero.list <- NULL
  rep <- model.sol$call$rep
  for (i in 1:rep){
    if (is.null(model.sol$weights[[i]])){
      zero.list <- c(zero.list, i)
    }
  }
  model.sol$weights[zero.list] <- NULL
  model.sol$generalized.weights[zero.list] <- NULL
  model.sol$net.result[zero.list] <- NULL
  model.sol$startweights[zero.list] <- NULL
  model.sol$call$rep <- rep-length(zero.list)
  
  min_rep = order(as.numeric(model.sol$result.matrix[
    rownames(model.sol$result.matrix)=="error",]))[1]
  AIC_x = model.sol$result.matrix[rownames(model.sol$result.matrix)=="aic",min_rep]
  BIC_x = model.sol$result.matrix[rownames(model.sol$result.matrix)=="bic",min_rep]
  ERR_x = model.sol$result.matrix[rownames(model.sol$result.matrix)=="error",min_rep]
  print(ERR_x)
  
  #  plot_neural_net(net = model.sol, model, iteration, min_rep)
  
  parameter <- data.frame(model.sol$weights[[1]][[1]])
  rownames(parameter) <- c("intercept", x_index) 
  # write.csv(parameter, paste0(simulation, "_", model, "_weight_NN2_layer1_", iteration, "_", i_fold, ".csv"))
  for (i_layer in 1:length(num.hidden.node)){
    parameter <- data.frame(model.sol$weights[[1]][[i_layer+1]])
    rownames(parameter) <- c("intercept", 1:num.hidden.node[i_layer]) 
    write.csv(parameter, paste0(simulation, "_", model, "_weight_NN2_layer",(i_layer+1) ,"_", iteration, "_", i_fold, ".csv"))
  }
  # write.csv(model.sol$result.matrix, 
            # file = paste0(simulation, "_", model, "_result_NN2_", iteration, "_", i_fold, ".csv"))
  
  
  insample_result <- predict(object = model.sol, newdata = train[,-1], rep = min_rep)
  if (is.null(scaled)){
    insample <- data.frame(A_y = train$y, E_y = insample_result)
  }else{
    maxs = scaled[1]
    mins = scaled[2]
    pr.nn <- insample_result*(maxs-mins)+mins
    train.r <- (train$y)*(maxs-mins)+mins
    insample <- data.frame(A_y = train.r, E_y = pr.nn)
  }
  #  print(insample)
  # write.csv(insample, paste0(simulation, "_insample_", model, "_NN2_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  MSPE = mean((insample[,1] - insample[,2])^2)
  cat("training MSE = ", MSPE,"\n")
  output_insample = data.frame(MSPE=MSPE)
  
  # out-of-sample prediction
  outsample_result <- predict(object = model.sol, newdata = test[,-1], rep = min_rep)
  if (is.null(scaled)){
    outsample <- data.frame(A_y = test$y, E_y = outsample_result)
  }else{
    maxs = scaled[1]
    mins = scaled[2]
    pr.nn <- outsample_result*(maxs-mins)+mins
    test.r <- (test$y)*(maxs-mins)+mins
    outsample <- data.frame(A_y = test.r, E_y = pr.nn)
  }
  #  print(outsample)
  # write.csv(outsample, paste0(simulation, "_outsample_", model, "_NN2_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  MSPE = mean((outsample[,1] - outsample[,2])^2)
  cat("test MSE = ", MSPE,"\n")
  output_outsample = data.frame(MSPE=MSPE)
  output <- c(output_insample, output_outsample)
  return(output)
}

class_matrix_to_vector <- function(modeled_output) {
  size <- nrow(modeled_output)
  classvalues <- c(rep(0, size))
  for (i in 1:size) {
    class <- which.max(modeled_output[i, ])
    if (modeled_output[i, class] > 0) {
      classvalues[i] <- class
    } else {
      classvalues[i] = NA
    }
  }
  return(classvalues)
}

accuracy <- function(confusion){
  frequent <- 0
  for (i in 1:ncol(confusion)){
    if (colnames(confusion)[i] %in% rownames(confusion)){
      frequent = frequent + confusion[which(rownames(confusion)==colnames(confusion)[i]), i]
    }
  }
  return(frequent / sum(confusion))
}

NN_classification <- function(train, test, col_remain, active_f, hidden, model, simulation, iteration, i_fold = 1){
  #          Neural Network
  #          --------------
  library(neuralnet)
  library(NeuralNetTools)
  
  nclasses <- length(unique(train$y)) # Number of classes
  training.size <- nrow(train) # Number of training samples
  test.size <- nrow(test) # Number of test samples
  # hidden <- c(10, 5) # hidden nodes
  learningfunc <- "rprop+"
  act.fct <-"logistic"
  iterations <- 1e+05
  is.regre <- TRUE
  err.fct <- "sse"
  still.times <- 1e+04
  
  data.size <- training.size + test.size
  
  cat("\n\ntraining size = ", training.size, " samples\n")
  cat("number of classes = ", nclasses, "\n")
  cat("hidden layer(s) = ", hidden, "\n")
  cat("learning function = ", learningfunc, "\n")
  
  classcol <- function (col, dim) {
    # converts a category to a column vector of dimension dim
    m <- matrix(0.1, dim, 1)
    m[col, 1] <- 0.9
    m
  }
  
  x_index <- paste0("X", col_remain)
  training.data <- data.frame(train[, x_index])
  colnames(training.data) <- x_index
  test.data <- data.frame(test[, x_index])
  colnames(test.data) <- x_index
  training.class <- as.numeric(train$y)
  test.class <- as.numeric(test$y)
  
  # convert each training.class to a column vector of length nclasses
  # where all entries are zero except for the one at column
  # training.class. Join all the column vectors (for each sample)
  # into a matrix of size training.size by nclasses.
  class_matrix <- sapply(training.class, classcol, nclasses)
  cdict <- paste0("c",1:nclasses)
  rownames(class_matrix) <- cdict
  # Assign variables c1,c2, and etc. to the nrows of the class_matrix.
  
  
  f <- as.formula(paste(paste(cdict, collapse = " + "), "~", paste(x_index, collapse = " + ")))
  
  model.sol <- neuralnet(f, data = cbind(t(class_matrix), training.data), 
                         hidden =  hidden, 
                         stepmax = iterations,
                         act.fct = act.fct, # "tanh", "logistic", sigmoid, softmax
                         algorithm = learningfunc, # "rprop-", "sag", "slr", OR "backprop", learningrate=0.01,
                         threshold=0.05,
                         #                        learningrate.factor = list(minus = 0.1, plus = 1.0),
                         rep = 3,                   # the training repetitions for the NN
                         err.fct = err.fct,
                         lifesign = "full",      # "minimal", "full" # verbose output during training
                         lifesign.step = 5000,
                         linear.output = is.regre, 
                         likelihood = TRUE,
                         still.times = still.times
  )
  
  zero.list <- NULL
  rep <- model.sol$call$rep
  for (i in 1:rep){
    if (is.null(model.sol$weights[[i]])){
      zero.list <- c(zero.list, i)
    }
  }
  model.sol$weights[zero.list] <- NULL
  model.sol$generalized.weights[zero.list] <- NULL
  model.sol$net.result[zero.list] <- NULL
  model.sol$startweights[zero.list] <- NULL
  model.sol$call$rep <- rep-length(zero.list)
  
  min_rep = order(model.sol$result.matrix[rownames(model.sol$result.matrix)=="error",])[1]
  ERR_x = model.sol$result.matrix[rownames(model.sol$result.matrix)=="error",min_rep]
  print(ERR_x)
  
  #  plot_neural_net(net = model.sol, model, iteration, min_rep)
  
  parameter <- data.frame(model.sol$weights[[1]][[1]])
  rownames(parameter) <- c("intercept", x_index) 
  # write.csv(parameter, paste0(simulation, "_", model, "_weight_NN2_layer1_", iteration, "_", i_fold, ".csv"))
  for (i_layer in 1:length(num.hidden.node)){
    parameter <- data.frame(model.sol$weights[[1]][[i_layer+1]])
    rownames(parameter) <- c("intercept", 1:num.hidden.node[i_layer]) 
    write.csv(parameter, paste0(simulation, "_", model, "_weight_NN2_layer",(i_layer+1) ,"_", iteration, "_", i_fold, ".csv"))
  }
  # write.csv(model.sol$result.matrix, 
            # file = paste0(simulation, "_", model, "_result_NN2_", iteration, "_", i_fold, ".csv"))
  
  insample_result <- predict(object = model.sol, newdata = training.data, rep = min_rep)
  classvalues <- class_matrix_to_vector(insample_result)
  
  insample <- table(training.class, classvalues)
  # write.csv(insample, paste0(simulation, "_insample_", model, "_NN2_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  training.accuracy <- accuracy(insample)
  cat("training data\n confusion matrix")
  print(insample)
  cat("training accuracy = ", training.accuracy,"\n")
  output_insample = data.frame(acc=training.accuracy)
  
  
  # out-of-sample prediction
  outsample_result <- predict(object = model.sol, newdata = test.data, rep = min_rep)
  classvalues <- class_matrix_to_vector(outsample_result)
  outsample <- table(test.class, classvalues)
  # write.csv(outsample, paste0(simulation, "_outsample_", model, "_NN2_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  cat("\n\ntest data\n confusion matrix")
  print(outsample)
  test.accuracy <- accuracy(outsample)
  cat("\n test accuracy = ", test.accuracy,"\n")
  output_outsample = data.frame(acc=test.accuracy)
  
  output <- c(output_insample, output_outsample)
  
  return(output)
}

SVM_classification <- function(train, test, f, model, simulation, iteration, i_fold = 1){
  ##########################################################################
  ## SVM
  core_time <- Sys.time()
  library(e1071)
  tuned <- tune.svm(f, data = train, gamma = 10^(-6:-1), cost = 10^(1:2),
                    type = "C-classification",  kernel = "radial")
  summary(tuned)
  model.tuned <- svm(f, data = train, gamma=tuned$best.parameters$gamma,
                     cost=tuned$best.parameters$cost,
                     type = "C-classification",  kernel = "radial")
  summary(model.tuned)
  # save(model.tuned, file = paste0(simulation, "_", model, "_model_SVM_", iteration, "_", i_fold, ".Rdata"))
  
  insample_result <- predict(model.tuned, train)
  insample <- table(train$y, insample_result)
  
  # write.csv(insample, paste0(simulation, "_insample_", model, "_SVM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  training.accuracy <- accuracy(insample)
  cat("training data\n confusion matrix")
  print(insample)
  cat("training accuracy = ", training.accuracy,"\n")
  output_insample = data.frame(acc=training.accuracy)
  
  
  # out-of-sample prediction
  outsample_result <- predict(object = model.tuned, newdata = test)
  outsample <- table(test$y, outsample_result)
  # write.csv(outsample, paste0(simulation, "_outsample_", model, "_SVM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  cat("\n\ntest data\n confusion matrix")
  print(outsample)
  test.accuracy <- accuracy(outsample)
  cat("\n test accuracy = ", test.accuracy,"\n")
  output_outsample = data.frame(acc=test.accuracy)
  
  output <- c(output_insample, output_outsample)
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  # write.csv(core_time, paste0(simulation, "_time_", model, "_SVM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  return(output)
}


SVM_regression <- function(train, test, f, scaled=NULL, model, simulation, iteration, i_fold=1){
  #          SVM
  #          --------------
  core_time <- Sys.time()
  library(e1071)
  tuned <- tune.svm(f, data = train, gamma = 10^(-6:-1), cost = 10^(1:2),
                    type = "eps-regression", kernel = "radial")
  print(summary(tuned))
  model.sol <- svm(f, data = train, gamma=tuned$best.parameters$gamma,
                   cost=tuned$best.parameters$cost)
  print(summary(model.sol))
  
  # save(model.sol, file = paste0(simulation, "_", model, "_model_SVM_", iteration, "_", i_fold, ".Rdata"))
  
  
  insample_result <- predict(object = model.sol, newdata = train)
  if (is.null(scaled)){
    insample <- data.frame(A_y = train$y, E_y = insample_result)
  }else{
    maxs = scaled[1]
    mins = scaled[2]
    pr.nn <- insample_result*(maxs-mins)+mins
    train.r <- (train$y)*(maxs-mins)+mins
    insample <- data.frame(A_y = train.r, E_y = pr.nn)
  }
  #  print(insample)
  # write.csv(insample, paste0(simulation, "_insample_", model, "_SVM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  MSPE = mean((insample[,1] - insample[,2])^2)
  cat("training MSE = ", MSPE,"\n")
  output_insample = data.frame(MSPE=MSPE)
  
  
  outsample_result <- predict(model.sol, test)
  if (is.null(scaled)){
    outsample <- data.frame(A_y = test$y, E_y = outsample_result)
  }else{
    maxs = scaled[1]
    mins = scaled[2]
    pr.nn <- outsample_result*(maxs-mins)+mins
    test.r <- (test$y)*(maxs-mins)+mins
    outsample <- data.frame(A_y = test.r, E_y = pr.nn)
  }
  #  print(outsample)
  # write.csv(outsample, paste0(simulation, "_outsample_", model, "_SVM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  MSPE = mean((outsample[,1] - outsample[,2])^2)
  cat("test MSE = ", MSPE,"\n")
  output_outsample = data.frame(MSPE=MSPE)
  print(output_outsample)
  output <- c(output_insample, output_outsample)
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  # write.csv(core_time, paste0(simulation, "_time_", model, "_SVM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  return(output)
}


CART_classification <- function(train, test, f, model, simulation, iteration, i_fold=1){
  core_time <- Sys.time()
  
  train$y <- as.numeric(as.factor(train$y))
  test$y <- as.numeric(as.factor(test$y))
  
  library(rpart) 
  library(maptree)
  model.sol <- rpart(f, method="class", data=train)
  #  draw.tree(model.sol)
  #  plotcp(model.sol)
  
  # save(model.sol, file = paste0(simulation, "_", model, "_model_CART_", iteration, "_", i_fold, ".Rdata"))
  
  insample_result <- predict(object = model.sol, newdata = train)
  insample_value <- class_matrix_to_vector(insample_result)
  insample <- table(train$y, insample_value)
  # write.csv(insample, paste0(simulation, "_insample_", model, "_CART_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  training.accuracy <- accuracy(insample)
  cat("training data\n confusion matrix")
  print(insample)
  cat("training accuracy = ", training.accuracy,"\n")
  output_insample = data.frame(acc=training.accuracy)
  
  
  # out-of-sample prediction
  outsample_result <- predict(object = model.sol, newdata = test)
  outsample_value <- class_matrix_to_vector(outsample_result)
  outsample <- table(test$y, outsample_value)
  # write.csv(outsample, paste0(simulation, "_outsample_", model, "_CART_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  cat("\n\ntest data\n confusion matrix")
  print(outsample)
  test.accuracy <- accuracy(outsample)
  cat("\n test accuracy = ", test.accuracy,"\n")
  output_outsample = data.frame(acc=test.accuracy)
  
  output <- c(output_insample, output_outsample)
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  # write.csv(core_time, paste0(simulation, "_time_", model, "_CART_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  return(output)
  
}


CART_regression <- function(train, test, f, scaled=NULL, model, simulation, iteration, i_fold=1){
  core_time <- Sys.time()
  library(rpart) 
  library(maptree)
  model.sol <- rpart(f, method="anova", data=train)
  #  draw.tree(model.sol)
  #  plotcp(model.sol)
  
  # save(model.sol, file = paste0(simulation, "_", model, "_model_CART_", iteration, "_", i_fold, ".Rdata"))
  
  
  insample_result <- predict(object = model.sol, newdata = train)
  if (is.null(scaled)){
    insample <- data.frame(A_y = train$y, E_y = insample_result)
  }else{
    maxs = scaled[1]
    mins = scaled[2]
    pr.nn <- insample_result*(maxs-mins)+mins
    train.r <- (train$y)*(maxs-mins)+mins
    insample <- data.frame(A_y = train.r, E_y = pr.nn)
  }
  #  print(insample)
  # write.csv(insample, paste0(simulation, "_insample_", model, "_CART_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  MSPE = mean((insample[,1] - insample[,2])^2)
  cat("training MSE = ", MSPE,"\n")
  output_insample = data.frame(MSPE=MSPE)
  
  
  outsample_result <- predict(model.sol, test)
  if (is.null(scaled)){
    outsample <- data.frame(A_y = test$y, E_y = outsample_result)
  }else{
    maxs = scaled[1]
    mins = scaled[2]
    pr.nn <- outsample_result*(maxs-mins)+mins
    test.r <- (test$y)*(maxs-mins)+mins
    outsample <- data.frame(A_y = test.r, E_y = pr.nn)
  }
  #  print(outsample)
  # write.csv(outsample, paste0(simulation, "_outsample_", model, "_CART_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  MSPE = mean((outsample[,1] - outsample[,2])^2)
  cat("test MSE = ", MSPE,"\n")
  output_outsample = data.frame(MSPE=MSPE)
  output <- c(output_insample, output_outsample)
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  # write.csv(core_time, paste0(simulation, "_time_", model, "_CART_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  return(output)
}


lightGBM_regression <- function(train, test, scaled=NULL, model, simulation, iteration, i_fold=1){
  core_time <- Sys.time()
  library(lightgbm)
  index_valid <- sample(x = 1:nrow(train), size = 0.2*nrow(train), replace = FALSE)
  train_x <- as.matrix(train[-index_valid,-1])
  train_y <- train[-index_valid,1]
  valid_x <- as.matrix(train[index_valid,-1])
  valid_y <- train[index_valid, 1]
  dtrain <- lgb.Dataset(data = train_x, label = train_y, is_sparse = FALSE,
                        colnames = x_index)
  dtest <- lgb.Dataset.create.valid(dataset=dtrain,
                                    data=valid_x,
                                    label=valid_y)
  valids <- list(test=dtest)
  
  params <- list(objective = "regression", metric = "l2") # L2 not twelve
  model.sol <- lgb.train(params=params,
                         data=dtrain,
                         valids = valids,
                         min_data = 1, # min data in a group
                         learning_rate = 0.1, # smaller,slower,maybe more accurate
                         nrounds = 1000,
                         early_stopping_rounds = 100 #if not better than last 20 rounds,stop
  )
  
  # save(model.sol, file = paste0(simulation, "_", model, "_model_lightGBM_", iteration, "_", i_fold, ".Rdata"))
  
  
  insample_result <- predict(model.sol, as.matrix(train[-1]))
  
  if (is.null(scaled)){
    insample <- data.frame(A_y = train$y, E_y = insample_result)
  }else{
    maxs = scaled[1]
    mins = scaled[2]
    pr.nn <- insample_result*(maxs-mins)+mins
    train.r <- (train$y)*(maxs-mins)+mins
    insample <- data.frame(A_y = train.r, E_y = pr.nn)
  }
  #  print(insample)
  # write.csv(insample, paste0(simulation, "_insample_", model, "_lightGBM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  MSPE = mean((insample[,1] - insample[,2])^2)
  cat("training MSE = ", MSPE,"\n")
  output_insample = data.frame(MSPE=MSPE)
  
  test_x <- as.matrix(test[,-1])
  test_y <- test[,1]
  outsample_result = predict(model.sol, test_x)
  if (is.null(scaled)){
    outsample <- data.frame(A_y = test$y, E_y = outsample_result)
  }else{
    maxs = scaled[1]
    mins = scaled[2]
    pr.nn <- outsample_result*(maxs-mins)+mins
    test.r <- (test$y)*(maxs-mins)+mins
    outsample <- data.frame(A_y = test.r, E_y = pr.nn)
  }
  #  print(outsample)
  # write.csv(outsample, paste0(simulation, "_outsample_", model, "_lighGBM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  MSPE = mean((outsample[,1] - outsample[,2])^2)
  cat("test MSE = ", MSPE,"\n")
  output_outsample = data.frame(MSPE=MSPE)
  output <- c(output_insample, output_outsample)
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  # write.csv(core_time, paste0(simulation, "_time_", model, "_lightGBM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  return(output)
}


lightGBM_classification <- function(train, test, model, simulation, iteration, i_fold=1){
  core_time <- Sys.time()
  library(lightgbm)
  train$y <- as.numeric(as.factor(train$y)) - 1
  test$y <- as.numeric(as.factor(test$y)) - 1
  index_valid <- sample(x = 1:nrow(train), size = 0.2*nrow(train), replace = FALSE)
  train_x <- as.matrix(train[-index_valid,-1])
  train_y <- train[-index_valid,1]
  valid_x <- as.matrix(train[index_valid,-1])
  valid_y <- train[index_valid, 1]
  dtrain <- lgb.Dataset(data = train_x, label = train_y, is_sparse = FALSE,
                        colnames = x_index)
  dtest <- lgb.Dataset.create.valid(dataset=dtrain, data=valid_x, label=valid_y)
  valids <- list(test=dtest)
  
  num_class <- length(unique(c(train$y, test$y)))
  params <- list(objective = "multiclass", metric = "multi_error", num_class = num_class)
  model.sol <- lgb.train(params=params,
                         data=dtrain,
                         valids = valids,
                         min_data = 1, # min data in a group
                         learning_rate = 0.1,
                         nrounds = 1000,
                         early_stopping_rounds = 100 #if not better than last 100 rounds,stop
  )
  
  # save(model.sol, file = paste0(simulation, "_", model, "_model_lightGBM_", iteration, "_", i_fold, ".Rdata"))
  
  insample_result <- predict(model.sol, as.matrix(train[-1]), reshape = TRUE)
  insample_value <- class_matrix_to_vector(insample_result)
  insample <- table(train$y, insample_value)
  colnames(insample) <- as.numeric(colnames(insample)) - 1
  
  # write.csv(insample, paste0(simulation, "_insample_", model, "_lightGBM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  training.accuracy <- accuracy(insample)
  cat("training data\n confusion matrix")
  print(insample)
  cat("training accuracy = ", training.accuracy,"\n")
  output_insample = data.frame(acc=training.accuracy)
  
  
  # out-of-sample prediction
  outsample_result <- predict(model.sol, as.matrix(test[,-1]), reshape = TRUE)
  outsample_value <- class_matrix_to_vector(outsample_result)
  outsample <- table(test$y, outsample_value)
  colnames(outsample) <- as.numeric(colnames(outsample)) - 1
  
  # write.csv(outsample, paste0( simulation, "_outsample_", model, "_lightGBM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  cat("\n\ntest data\n confusion matrix")
  print(outsample)
  test.accuracy <- accuracy(outsample)
  cat("\n test accuracy = ", test.accuracy,"\n")
  output_outsample = data.frame(acc=test.accuracy)
  
  output <- c(output_insample, output_outsample)
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  # write.csv(core_time, paste0(simulation, "_time_", model, "_lightGBM_", iteration, "_", i_fold, ".csv"), row.names = TRUE)
  return(output)
  
}

###############################################################################

## The Gram-Schimidt Orthogonalization in CDNFS for iterpretation

###############################################################################

Gram_Schmidt_delete_visual <- function(D_GS, x, alpha = 0.05){
  N = nrow(x)
  p = ncol(x)
  index_x <- array(0, dim = p)
  u <- x
  
  for (i in 1:p){
    for (j in 1:ncol(D_GS)){
      d <- t(u[,i]) %*% D_GS[,j] / (t(D_GS[,j]) %*% D_GS[,j])
      u[,i] <- u[,i]-c(d)*D_GS[,j]
    }
  }
  # non-supervised G-S
  t_value<- rep(0, times = p)
  for (i in 1:p){
    if (var(u[,i])==0 | var(x[,i])==0){
      t_value[i] = 0
    }else{
      r = cor(u[,i], x[,i])
      # t_value[i] = r*sqrt(N-2)/sqrt(1-r^2)
      t_value[i] = var(u[,i])
    }
  }
  
  Gk <- colnames(u)[which(t_value < alpha)]
  # print(round(t_value, 3))
  Gk_data <- data.frame(Gk, alpha = t_value[t_value < alpha])
  print(Gk)
  
  # t_bound = qt(1-alpha, df = N-2)
  # Gk <- colnames(u[, which(abs(t_value) < t_bound)])
  # print(Gk)
  # print(t_value[abs(t_value) < t_bound])
  
  return(Gk_data)
}

###############################################################################

## The Main function of the CDNFS used in the interpretation

###############################################################################
CDNFS_visual <- function(data, class_col = 1, d = NULL){
  core_time <- Sys.time()
  # reduce the dimension with Gram-Schmidt transformation. Train the parameters 
  # using only the training set
  y <- data.frame(class = data[, class_col])
  x = data.frame(data[,-class_col])
  name_x <- colnames(x)
  y$class <- as.numeric(y$class)
  p = ncol(x)
  colnames(x) <- paste0("V", c(1:p))
  print(dim(x))
  
  # delete the variables without information (sd(x) < 1e-5)
  library(infotheo)
  entropy_x <- apply(infotheo::discretize(x), 2, entropy)
  non_zero <- which( entropy_x > 1e-2)
  #  non_zero <- which(apply(x, 2, function(x) sd(x) > 1e-2))
  
  zero <- colnames(x[,-non_zero])
  zero_info <- entropy_x[entropy_x <= 1e-2]
  zero_data <- data.frame(zero, entropy=zero_info)
  filename = paste0(simulation,"_delete_lack_info_", iteration, ".csv")
  write.csv(zero_data, filename, row.names = FALSE)
  
  name_x <- colnames(x)[non_zero]
  x <- data.frame(x[, non_zero])
  colnames(x) <- name_x
  # print(zero)
  print(dim(x))
  
  x <- scale(x)
  N <- nrow(x)
  p = ncol(x)
  z = NULL
  stepT = c()
  delete_variables <- NULL
  k_GS <- 1
  
  if (is.null(d)) d.stop <- p
  else d.stop <- d
  
  while (p>0 && length(z)<d.stop){
    cor_yx <- array(0, dim = p)
    cor_S <- codec_S(y[,1], z)
    if (cor_S <= 0) break
    for (i in 1:p){
      cor_yx[i] <- codec_Q(y[,1], x[,i], z)/cor_S
    }
    k = order(cor_yx, decreasing = TRUE)[1]
    t = max(cor_yx)
    stepT = c(stepT, t)
    D = data.frame(x[,k])
    colnames(D) <- colnames(x)[k]
    print(c(round(t, 3), colnames(D)))
    if (is.null(z)){
      z = data.frame(D)
      
      D_GS <- z
    }else{
      name_z <- colnames(z)
      z = cbind(z, D)
      colnames(z) <- c(name_z, colnames(D))
      
      D_GS <- cbind(D_GS, D)
      colnames(D_GS) <- colnames(z)
      for (i in 1:(ncol(D_GS)-1)){
        d_GS <- t(D) %*% D_GS[,i] / (t(D_GS[,i]) %*% D_GS[,i])
        D_GS[, ncol(D_GS)] <- D_GS[, ncol(D_GS)] - c(d_GS) * D_GS[, i]
      }
    }
    if (t<=0 && is.null(d)) break
    if(ncol(x)>2){
      x <- x[, -k]
    }else
      if (ncol(x)==2){
        xnames <- colnames(x)
        x <- data.frame(x = x[,-k])
        colnames(x) <- xnames[-k]
      }else break
    
    # delete the variables which is highly correlated with the selected x[,k]
    # using unsupervised Gram-Schimidt
    Gk_data <- Gram_Schmidt_delete_visual(D_GS, x, alpha = 0.05)
    if (length(Gk_data)>0){
      Gk <- Gk_data$Gk
      # print(length(Gk))
      # x[, Gk] <- list(NULL)
      x <- x[, (!colnames(x) %in% Gk), drop = FALSE]
      delete_variables <- rbind(delete_variables, Gk_data)
    }
    p = ncol(x)
    
    # print(ncol(z))
    print(dim(x))
    
    filename = paste0(simulation,"_variables_selected_codec_GS_", iteration, ".csv")
    write.csv(data.frame(col_remain = colnames(z), stepT = stepT), filename, row.names = FALSE)
    
    filename = paste0(simulation,"_variables_deleted_GS_", iteration, "_", k_GS, ".csv")
    write.csv(Gk_data, filename, row.names = FALSE)
    k_GS <- k_GS + 1
  }
  
  
  # Save the parameters of the Gram-Schimidt transformation.
  col_remain <- colnames(z)
  library(stringr)
  col_remain <- as.numeric(str_remove(string = col_remain, pattern = "V"))
  output_col <- data.frame(col_remain, stepT = stepT)
  
  # print(col_remain)
  # print(delete_variables)
  if (is.null(d))
    output_col <- output_col[-length(col_remain),]
  core_time <- difftime(Sys.time(), core_time, units = 'sec')
  results <- list(feature = output_col, 
                  time = core_time)
  return(results)
}
