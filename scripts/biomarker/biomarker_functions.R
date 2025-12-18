
#' Fits a simple linear model by regressing y over each column of X.
#'
#' @param X : matrix of n by m, it can have NA's in it, columns should be named.
#' @param Y : matrix of n by p, rows are in the same order of with the rows of X.
#' @param v.th : Minimum variance for the columns of X, columns with smaller variances will be dropped.
#' @param n.min : Minimum number of finite pairs between columns of X and Y, column pairs not satisfying this condition will be dropped.
#' @param rank.max : For each column of Y, only the top rank.max elements will be returned.
#' @param q.val.max : Associations with a q_val > q.val.max will be dropped.
#' @param regression_coef : Should the regression coefficients be calculated
#' @param stability_score : Should the stability scores be calculated (valid only if regression_coef = TRUE) 
#' @param ns.min : If the stability scores are calculated, associations with score less than ns.min will be dropped
#' 
#'
#' @return A data frame with: x (corresponding column of X), y (corresponding column of Y), correlation_coeff (Pearson correlation), 
#'         regression_coeff (regression coefficient), p_val / q_val (homoskedastic p-value / q-value),
#'         n (number of non-na samples), ns (a stability score to represent how many points of y needs to be replaced with
#'         the mean value to flip the sign of regression coefficient), rank (rank of the significance for each column of Y), 
#
#' @export
#'
#' @examples
#'
linear_model <- function(X, Y, v.X.min = 0.0025, n.min = 25, rank.max = 250, q.val.max = 0.1, regression_coef = TRUE, stability_score = FALSE, ns.min = 3) {
  require(tidyverse)
  require(matrixStats)
  require(WGCNA)
  
  X <- X[, (matrixStats::colVars(X, na.rm = T) >= v.X.min) & (colSums(!is.na(X)) >= n.min), drop = FALSE]
  if (ncol(X)==0){
    print(paste("returning empty frame since data frame has ncol of X=",ncol(X)))
    return (tibble(x = character(), y = character(), feat = character())) 
    #' empty data frame if no columns left. 
    #' for eg, this can occur when we are finding biomarkers within a lineage and we try to run biomarkers for lineages
    #' then we have perfect agreement and so v.X will be 0.
  }
  print("Calculating correlation coefficients and q-values...")
  res <- WGCNA::corAndPvalue(X,Y, use = "p")[c(1,2,5)] %>%
    reshape2::melt() %>% 
    tidyr::pivot_wider(names_from = L1, values_from = value) %>% 
    dplyr::filter(nObs >= n.min) %>% 
    dplyr::group_by(Var2) %>%  
    dplyr::mutate(q = p.adjust(p, method = "BH")) %>% 
    dplyr::arrange(q) %>% dplyr::mutate(rank = 1:n()) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(q <= q.val.max) %>%
    dplyr::group_by(Var2, sign(cor)) %>% 
    dplyr::arrange(desc(abs(cor))) %>% dplyr::mutate(rank = pmin(rank, 1:n())) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank <= rank.max) %>%
    dplyr::rename(x = Var1, y = Var2, rho = cor, p.val = p, q.val = q, n = nObs) %>%
    dplyr::mutate(x = as.character(x), y = as.character(y)) %>%
    dplyr::distinct(x, y, rho, p.val, q.val, n, rank) %>%
    dplyr::rename(correlation_coef = rho, p_val = p.val, q_val = q.val)
  
  # Do we need the regression coefficients?
  if(regression_coef & (nrow(res) > 0)){
    print("Calculating regression coefficients...")
    X_ <- X[, res$x]; Y_ <- Y[, res$y]
    mask <- !is.finite(X_ * Y_)
    X_[mask] <- NA; Y_[mask] <- NA
    Y_ <- scale(Y_, center = TRUE, scale = FALSE)
    X_ <- scale(X_, center = TRUE, scale = FALSE)
    
    res$var.x <- colMeans(X_^2, na.rm = T)
    res$var.y <- colMeans(Y_^2, na.rm = T) 
    
    res <- res %>% 
      dplyr::mutate(regression_coef = correlation_coef * sqrt(var.y/var.x)) %>% 
      dplyr::select(-var.y)
    
    #' Robustness score that determines the number of samples, that when modified
    #' flips the sign of corr/reg coefficient
    if(stability_score){
      print("Calculating stability scores...")
      S <- X_ * Y_
      n <- colSums(!is.na(X_))
      ns <- apply(S,2, function(s) sum(cumsum(sort(s/sum(s, na.rm = T), decreasing = T)) < 1, na.rm = T)) + 1
      res$ns <- pmin(ns, n - apply(X_, 2, function(x) max(table(x))))
      res <- dplyr::filter(res, ns >= ns.min)
    }
  }
  return(res)
}


#' Calculates and returns univariate analysis results by correlating each column of Y with features sets from depmap.org.
#'
#' @param Y : Matrix n x m, make sure rownames are depmap_id's and columns are named. There can be NA's.
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param features : You can give a subset of available feature sets, but if left as NULL, it computes for everything.
#' @param homoskedastic : Should homoskedastic or robust q-values be used for ranking and filtering. Default is homoskedastic.
#' @param n.X.min : Results with less thana given sample size are dropped, default is 100
#' @param n_stable.min : Results that can be flipped with less than n_stable.min data points are dropped, the default is 3
#' @param q_val.max : Results with q-values less than q_val.max are dropped, default is 0.2
#' @param parallel : Should multiple cores to be used to decrease the compute time. Default is FALSE.
#'
#' @return Returns a data-table with each row corresponds to a particular (feature, feature_set, y) triplet. See linear_model for the other columns.
#' @export
#'
#' @examples
#'
univariate_biomarker_table <- function(Y, file,
                                       features = NULL, n.X.min = 100,
                                       v.X.min = 0.0025,
                                       n_stable.min = 3, q_val_max = .1,
                                       regression_coef = TRUE,
                                       stability_score = TRUE, rank.max = 250){
  require(tidyverse)
  require(magrittr)
  require(rhdf5)
  require(WGCNA)

  
  if(!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  
  if(is.null(features)){
    features <-  substr(setdiff(h5ls(file)$group, c("/",  "/Lineage_Table")),2,100)
  }else{
    features <- intersect(features, substr(setdiff(h5ls(file)$group, c("/",  "/Lineage_Table")),2,100))
  }
  print(features)
  
  RESULTS <- list()
  for(feat in features){
    
    
    X <- read_dataset(file , feat)
    cl = intersect(rownames(X), rownames(Y))
    X <- X[cl, , drop = FALSE] ## drop=false to handle single row situation
    
    if((nrow(X) >= n.X.min) & (ncol(X) > 0)){ ## use nrow and ncol in case X is null
      
      print(paste0(feat, " - ", dim(X)[1] ,'x', dim(X)[2]))
      
      RESULTS[[feat]] <- linear_model(X = X, Y = Y[cl, , drop = FALSE],
                                                        v.X.min = v.X.min, n.min = n.X.min, 
                                                        rank.max = rank.max, q.val.max = q_val_max, 
                                                        regression_coef = regression_coef, stability_score = stability_score, ns.min = n_stable.min) %>%
        dplyr::rename(feature = x) %>%
        dplyr::mutate(feature_set = feat) 
      print(feat)
    }else{
      print(paste("Skipping feature set", feat, "because it has too few samples or columns for the cell lines.",
                  nrow(X),ncol(X) ))
      }
  }
  
  RESULTS <- dplyr::bind_rows(RESULTS)
  return(RESULTS)
}





#' Fits a random forest from X to y using k-fold cross validation. Returns
#'
#' @param X : Feature matrix.
#' @param y : Response vector.
#' @param k : Number of cross-validation folds
#' @param vc : Variance threshold, columns of X with lower variance than this threshold are dropped.
#' @param lm : Parameter for gausscov package for feature selection within each fold.
#' @param p0 : Parameter for gausscov package for feature selection within each fold.
#' @param folds : Fold structure for the cross-validation. If NULL, stratified folds are created
#' @param X.test : If provided the the predictions are provided as part of the outputs.
#'
#' @return model_table: Model performance parameters and feature importances.
#'         predictions: Predictions.
#'         folds = Cross-validation folds.
#'         y.test.hat: Test predictions if X.test is provided.
#'
#' @export
#'
#' @examples
random_forest <- function (X, y, k = 5, vc = 0.01, lm = 25, p0 = 0.01, folds = NULL, X.test = NULL) {
  y.clean <- y[is.finite(y)]
  cl <- sample(dplyr::intersect(rownames(X), names(y.clean)))
  X.clean <- X[cl, apply(X[cl, ], 2, function(x) all(is.finite(x)))]
  if(is.null(folds)){
    N = floor(length(cl)/k)
    folds <- sapply(1:k, function(kx) cl[((kx - 1) * N + 1):(kx * N)])
  }

  X.clean <- X.clean[cl, ]
  y.clean <- y.clean[cl]
  colnames(X.clean) %<>% make.names()
  yhat_rf <- rep(NA, length(y.clean))
  names(yhat_rf) <- names(y.clean)
  SS = tibble()

  if(!is.null(X.test)){
    y.test.hat = list()
    px = 1
  }

  for (kx in 1:k) {
    test <- intersect(folds[,kx], cl)
    train <- dplyr::setdiff(cl, test)
    X_train <- X.clean[train, , drop = F]
    X_test <- X.clean[test, , drop = F]
    X_train <- X_train[, apply(X_train, 2, var) > vc, drop = F]
    b <- tryCatch(gausscov::f2st(y = y.clean[train], x = X_train, lm = lm, p0 = p0, sub = FALSE), error = function(e) NULL)
    if(is.null(b)){
      cc = cor(X_train, y.clean[train])
      features <- union(rownames(cc[rank(-cc) <= 100, ,drop = FALSE]), rownames(cc[rank(cc) <= 100, ,drop = FALSE]))
    }else{
      features <- colnames(X_train)[b[[1]][, 2]]
    }

    if (length(features) > 0) {
      X_train <- X_train[, features, drop = F]
      rf <- ranger::ranger(y = y.clean[train], x = X_train,
                           importance = "impurity")
      yhat_rf[test] <- predict(rf, data = as.data.frame(X_test[, colnames(X_train), drop = F]))$predictions

      if(!is.null(X.test)){
        X.test.temp <- X.test[, colnames(X_train), drop = F]
        X.test.temp <- X.test.temp[apply(is.finite(X.test.temp), 1, all), , drop = FALSE]
        y.test.hat[[px]] = tibble(y.hat = predict(rf, data = as.data.frame(X.test.temp))$predictions,
                                  sample_id = rownames(X.test))
        px = px + 1
      }

      ss <- tibble::tibble(feature = names(rf$variable.importance),
                           RF.imp = rf$variable.importance/sum(rf$variable.importance),
                           fold = kx)
      SS %<>% dplyr::bind_rows(ss)
    }
  }
  if (nrow(SS) > 0) {
    RF.importances <- SS %>%
      dplyr::distinct(feature, RF.imp, fold) %>%
      reshape2::acast(feature ~ fold, value.var = "RF.imp", fill = 0)

    RF.table <- tibble::tibble(feature = rownames(RF.importances),
                               RF.imp.mean = rowMeans(RF.importances), 
                               RF.imp.sd = apply(RF.importances, 1, function(x) sd(x, na.rm = T)),
                               RF.imp.stability = apply(RF.importances, 1, function(x) sum(x > 0)/k)) %>%
      dplyr::filter(feature != "(Intercept)") %>%
      dplyr::arrange(desc(RF.imp.mean)) 
    
    mse <- mean((yhat_rf - y.clean)^2, na.rm = T)
    mse.se <- sqrt(var((yhat_rf - y.clean)^2, na.rm = T))/sqrt(length(y.clean))
    r2 <- 1 - (mse/var(y.clean, na.rm = T))
    ps <- cor(y.clean, yhat_rf, use = "pairwise.complete.obs")
    RF.table %<>% dplyr::mutate(MSE = mse, MSE.se = mse.se,
                                R2 = r2, PearsonScore = c(ps))

    if(is.null(X.test)){
      return(list(model_table = RF.table, predictions = yhat_rf, folds = folds))
    }else{
      y.test.hat <- dplyr::bind_rows(y.test.hat)
      return(list(model_table = RF.table, predictions = yhat_rf, folds = folds, y.test.hat = y.test.hat))
    }
  }
  else {
    if(is.null(X.test)){
      return(list(model_table = tibble(), predictions = tibble(), folds = folds))
    }
    else{
      return(list(model_table = tibble(), predictions = tibble(), folds = folds, y.test.hat = tibble()))
    }
  }
}


#' Exports individual datasets from depmap_datasets.h5 file. You can specify the row names either as ccle_names or depmap_ids.
#'
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param dataset : The dataset you want to export. You can check the available datasets with substr(setdiff(rhdf5::h5ls(file)$group, "/"),2,100)
#' @param rownames_depmap_ids : Default TRUE, you can get the rownames as ccle_names by switching to FALSE.
#'
#' @return The requested data matrix.
#' @export
#'
#' @examples
read_dataset <- function(file = "/cmap/obelix/pod/biomarker/depmap_public_24q4.h5", dataset,
                         rownames_depmap_ids = TRUE) {
  require(rhdf5)
  if(word(file, sep = fixed("://")) %in% c("s3", "http", "https")){
    s3 = TRUE
    print(paste0("Reading ", file, " from S3"))
  } else{
    s3 = FALSE
    print(paste0("Reading ", file, " from local"))
  }
  X <- h5read(file, name = paste0(dataset, "/mat"), s3 = s3)
  row_meta <- h5read(file, name = paste0(dataset, "/row_meta"), s3 = s3)
  column_meta <- h5read(file, name = paste0(dataset, "/column_meta"), s3 = s3)
  colnames(X) <- column_meta$column_name
  if(rownames_depmap_ids){
    rownames(X) <- row_meta$ModelID
  }else{
    rownames(X) <- row_meta$CCLEName
  }

  X <- X[rownames(X) != "NA", colnames(X) != "NA", drop = FALSE]
  X <- X[!duplicated(rownames(X)), !duplicated(colnames(X)), drop = FALSE]
  return(X)
}


#' Exports individual columns of a selected dataset from depmap_datasets.h5 file. You can specify the row names either as ccle_names or depmap_ids.
#'
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param dataset : The dataset you want to export. You can check the available datasets with substr(setdiff(rhdf5::h5ls(file)$group, "/"),2,100)
#' @param feature_names : Particular column names to pull from the given dataset.
#' @param rownames_depmap_ids : Default TRUE, you can get the rownames as ccle_names by switching to FALSE.
#'
#' @return
#' @export
#'
#' @examples
read_features <- function(file = "/cmap/obelix/pod/biomarker/depmap_public_24q4.h5",
                          dataset, feature_names = NULL,
                          rownames_depmap_ids = TRUE) {
  require(rhdf5)
  if(word(file, sep = fixed("://")) %in% c("s3", "http", "https")){
    s3 = TRUE
  } else{
    s3 = FALSE
  }

  row_meta <- h5read(file, name = paste0(dataset, "/row_meta"), s3 = s3)
  column_meta <- h5read(file, name = paste0(dataset, "/column_meta"), s3 = s3)

  if(!is.null(feature_names)){
    column_meta <- dplyr::filter(column_meta, column_name %in% feature_names)
  }

  if(nrow(column_meta) > 0){
    X <- h5read(file, name = paste0(dataset, "/mat"),
                index = list(NULL, column_meta$ix), s3 = s3)
    colnames(X) <- column_meta$column_name
    if(rownames_depmap_ids){
      rownames(X) <- row_meta$ModelID
    }else{
      rownames(X) <- row_meta$CCLEName
    }
    X <- X[rownames(X) != "NA", colnames(X) != "NA", drop = FALSE]
    X <- X[!duplicated(rownames(X)), !duplicated(colnames(X)), drop = FALSE]
    return(X)
  } else{
    return(NULL)
  }
}




#' Creates and returns two concatenated depmap feature matrices for given response and confounder matrices.
#' Returned matrices shares the same row names (depmap_id's). The rownames are limited to the intersection of the rownames of Y and W.
#' Additionally, any columns or rows that has more than 10% NA's are dropped.
#'
#' @param Y : Response matrix, rows are indexed by depmap_id's.
#' @param W : Confounder matrix, rows are indexed by depmap_id's.
#' @param file : File path for the compressed depmap datasets (.h5 file).
#'
#' @return Returns a list of three matrices X.DNA, X.RNA, and X.CRISPR. X.DNA contains CopyNumber, Mutation, Fusion, Lineage, and OmicsSignatures features as columns.
#'         X.RNA additionally includes the RNA expression data as well, lastly X.CRISPR adds the CRISPR profiles for the overlapping rows.
#' @export
#'
#' @examples
RF_feature_sets <- function(Y, W = NULL, file) {
  require(magrittr)

  if(!is.null(W)){
    cl = intersect(rownames(Y), rownames(W))
    Y <- Y[cl, , drop = FALSE]
    W <- W[cl, , drop = FALSE]
    colnames(W) %<>% paste0("CONF_", .)
  }


  CN <- read_dataset(file = file, dataset = 'CopyNumber')
  colnames(CN) %<>% paste0("CN_", .)
  MUT <- read_dataset(file = file, dataset = 'Mutation')
  colnames(MUT) %<>% paste0("MUT_", .)
  FUS <- read_dataset(file = file, dataset = 'Fusion')
  colnames(FUS) %<>% paste0("FUS_", .)
  LIN <- read_dataset(file = file, dataset = 'Lineage')
  colnames(LIN) %<>% paste0("LIN_", .)

  cl = intersect(rownames(MUT), rownames(CN)) %>%
    intersect(rownames(FUS)) %>%
    intersect(rownames(LIN)) %>%
    # intersect(rownames(SIG)) %>%
    intersect(rownames(Y))

  X.DNA <- cbind(CN[cl, ], MUT[cl, ]) %>%
    cbind(FUS[cl, ]) %>%
    cbind(LIN[cl, ])

  if(!is.null(W)){
    X.DNA <- cbind(X.DNA, W[cl, , drop = FALSE])
  }

  X.DNA <- X.DNA[rowMeans(is.finite(X.DNA)) > 0.9, ,drop = FALSE]

  rm(CN, MUT, FUS, LIN)
  EXP <- read_dataset(file = file, dataset = 'Expression')
  EXP <- EXP[rowMeans(is.finite(EXP)) > 0.9, ,drop = FALSE]
  colnames(EXP) %<>% paste0("EXP_", .)
  
  cl = intersect(rownames(X.DNA), rownames(EXP))

  X.DNA <- X.DNA[cl, ]
  X.RNA <- EXP[cl, ]
  rm(EXP)

  CRISPR <- read_dataset(file = file, dataset = "CRISPR")
  CRISPR <- CRISPR[rowMeans(is.finite(CRISPR)) > 0.9, ,drop = FALSE]
  colnames(CRISPR) %<>% paste0("CRISPR_", .)
  
  cl = intersect(cl, rownames(CRISPR))
  X.CRISPR <- CRISPR[cl, ]
  rm(CRISPR)

  colnames(X.DNA) %<>% make.names()
  colnames(X.RNA) %<>% make.names()
  colnames(X.CRISPR) %<>% make.names()

  return(list(X.DNA = X.DNA,  X.RNA = X.RNA, X.CRISPR = X.CRISPR))
}


#' Fit two random forest models (DNA and DNA+RNA) for each column of Y and return the
#' summary results: Cross validated R2, PearsonScore, feature importances (mean, sd, stability)
#'
#' @param Y : Response matrix, rows are depmap_id's. 2 models are fit for each column. Columns should be named.
#' @param W : Confounder matrix. Should not contain and NA's, rows are indexed by depmap_ids, and aligned with Y. Default is NULL.
#' @param file : Path to the h5 file containing depmap datasets.
#' @param k : Number of cross validation folds. Default is 10.
#'
#' @return Returns a data frame containing the feature importances and model properties for each column of Y (y) and model (DNA or DNA+RNA)
#' @export
#'
#' @examples
multivariate_biomarker_table <- function(Y, W = NULL, file, k = 10) {
  X <- RF_feature_sets(Y, W = W, file = file)
  cl <- intersect(rownames(X$X.DNA), rownames(X$X.RNA))

  rf.DNA <- list(); rf.RNA <- list(); rf.CRISPR <- list()
  for(ix in 1:dim(Y)[2]){

    y = Y[,ix]; y = y[is.finite(y)]
    cl_ = intersect(cl, names(y))

    rf_DNA <- random_forest(X$X.DNA[cl_, ], y[cl_], k = k)
    rf_RNA <- random_forest(cbind(X$X.DNA[cl_, ], X$X.RNA[cl_, ]), y[cl_], folds = rf_DNA$folds, k = k)
    
    cl_ = intersect(cl_, rownames(X$X.CRISPR))
    
    
    rf_CRISPR <- random_forest(cbind(X$X.DNA[cl_, ], X$X.RNA[cl_, ], X$X.CRISPR[cl_, ]), y[cl_], k = k)

    rf.DNA[[ix]] <- rf_DNA$model_table %>%
      dplyr::mutate(y = colnames(Y)[ix],
                    model = "dna") %>%
      dplyr::filter(RF.imp.stability >= 0.5) %>% 
      dplyr::mutate(rank = 1:n())
    

    rf.RNA[[ix]] <- rf_RNA$model_table %>%
      dplyr::mutate(y = colnames(Y)[ix],
                    model = "dna_rna") %>%
      dplyr::filter(RF.imp.stability >= 0.5) %>% 
      dplyr::mutate(rank = 1:n())
    
    
    rf.CRISPR[[ix]] <- rf_CRISPR$model_table %>%
      dplyr::mutate(y = colnames(Y)[ix],
                    model = "dna_rna_crispr") %>%
      dplyr::filter(RF.imp.stability >= 0.5) %>% 
      dplyr::mutate(rank = 1:n())

    print(paste0(ix, " out of ", dim(Y)[2]))
  }

  RF <- dplyr::bind_rows(dplyr::bind_rows(rf.DNA), dplyr::bind_rows(rf.RNA), dplyr::bind_rows(rf.CRISPR))

  return(RF)
}








#' Function that creates the multivariate biomarker results from the input file.
#'
#' @param in_path : Path for the input file.
#' @param out_path : Path to write the output file.
#' @param input_file_name : Name of the input file.
#' @param output_file_name : Name of the output file.
#' @param cell_line_column : Index for of cell lines. Input file should have this column. The strong default is depmap_id.
#' @param depmap_file : File path for the .h5 files for the depmap datasets.
#' @param treatment_columns : Set of columns of the input file to identify each profile.
#' @param response_column : Column of the input file that has the respose file.
#' @param aggregate_function : How to aggregate if the treatment columns and cell line column doesn't uniquely identifies the response value.
#' @param transform_function : Should the response variable be transformed before the analysis. Default is identity, a common choice is log2.
#'
#' @return
#' @export
#'
#' @examples
#' create_multivariate_biomarker_table( in_path = "data/MTS_SEQ003_KF/",
#'                                      treatment_columns = c("x_project_id", "pert_plate", "pert_name", "pert_dose"),
#'                                      depmap_file = "~/Downloads/mts_sequencing_analysis/data/depmap_datasets_internal_24q4.h5")
#'
#' create_multivariate_biomarker_table( in_path = "data/MTS_SEQ003_KF/", input_file_name = "drc.csv",
#'                                      output_file_name = "l2auc_multivariate_biomarkers",
#'                                      treatment_columns = c("x_project_id", "pert_plate", "pert_name"),
#'                                      response_column = "auc", transform_function = log2,
#'                                      depmap_file = "~/Downloads/mts_sequencing_analysis/data/depmap_datasets_internal_24q4.h5")
#'
create_multivariate_biomarker_table <- function(in_path, out_path = NULL,
                                                output_file_name = "l2fc_multivariate_biomarkers",
                                                depmap_file,
                                                treatment_columns = c("pert_id", "x_project_id", "pert_name", "pert_plate", "pert_dose"),
                                                response_column = "median_l2fc", aggregate_function = median, transform_function = function(x){x}) {
  require(data.table)
  require(tidyverse)
  require(rlang)

  # check if the input path exists
  print(paste0("DEBUG inpath: ", in_path))
  if(!file.exists(in_path)){
    stop("Input path does not exist!")
  }

  # load the input file
  input_file <- read_data_table(in_path)

  # check if the input table populated
  if(nrow(input_file) == 0){
    stop(paste0("Input file for multivariate biomarker is empty!"))
  }

  # check for missing critical columns -  soft warning
  necessary_columns <- unique(c(response_column, treatment_columns, "depmap_id"))

  if(any(!necessary_columns %in% colnames(input_file))){
    stop(paste0(paste0(setdiff(necessary_columns, colnames(input_file)), collapse = ", "),
                " are missing from multivariate biomarker input file!"))
  }

  # cast the responses into a matrix
  M <- input_file %>%
    dplyr::filter(is.finite(.data[[response_column]])) %>%
    tidyr::unite(cn, all_of(treatment_columns), sep = "::") %>%
    reshape2::acast(depmap_id ~ cn, value.var = response_column, fun.aggregate = aggregate_function) %>%
    transform_function()

  # generate the biomarker table
  multivariate_biomarker_table <- multivariate_biomarker_table(Y = M, file = depmap_file, k = 10)

  # Export the biomarker table -----
  print(paste0("Writing the multivariate output file to ", paste0(out_path, "/", output_file_name)))
  input_file %>%
    dplyr::select(all_of(treatment_columns)) %>%
    dplyr::distinct() %>% 
    tidyr::unite(y, all_of(treatment_columns), sep = "::", remove = FALSE) %>% 
    dplyr::inner_join(multivariate_biomarker_table) %>% 
    dplyr::select(-y) %>%
    write_out_table(path = file.path(out_path, output_file_name))
}



#' Function that creates the univariate biomarker results from the input file.
#'
#' @param in_path : Path for the input file.
#' @param out_path : Path to write the output file.
#' @param input_file_name : Name of the input file.
#' @param output_file_name : Name of the output file.
#' @param cell_line_column : Index for of cell lines. Input file should have this column. The strong default is depmap_id.
#' @param depmap_file : File path for the .h5 files for the depmap datasets.
#' @param treatment_columns : Set of columns of the input file to identify each profile.
#' @param response_column : Column of the input file that has the respose file.
#' @param aggregate_function : How to aggregate if the treatment columns and cell line column doesn't uniquely identifies the response value.
#' @param transform_function : Should the response variable be transformed before the analysis. Default is identity, a common choice is log2.
#' @param features : Set of depmap datasets to be used in the analysis. The default NULL makes sure all the datasets in the dpemap_file are used.
#' @param min_sample_size : Minimum sample size for a feature to be included in the analysis. The default is 100.
#' See the description of univariate_biomarker_table() for the remaining parameters.
#'
#' @return
#' @export
#'
#' @examples
create_univariate_biomarker_table <- function(in_path, out_path,
                                              output_file_name = "l2fc_univariate_biomarkers", depmap_file,
                                                treatment_columns = c("pert_id", "x_project_id", "pert_name", "pert_plate", "pert_dose"),
                                                response_column = "median_l2fc", aggregate_function = median, transform_function = function(x){x},
                                              features = NULL, min_sample_size = 100,
                                              regression_coef = TRUE, stability_score = TRUE, rank.max = 250,
                                              min_x_variance = 0.0025, n_stable_min = 3, q_val_max = .1
                                              ) {
  require(data.table)
  require(tidyverse)
  require(rlang)


  
  # check if the input path exists
  if(!file.exists(in_path)){
    stop("Input path does not exist!")
  }

  # load the input file
  input_file <- read_data_table(in_path)

  # check if the input table populated
  if(nrow(input_file) == 0){
    stop(paste0("Input file for univariate biomarker is empty!"))
  }

  # check for missing critical columns
  necessary_columns <- unique(c(response_column, treatment_columns, "depmap_id"))
  
  # soft warning
  if(any(!necessary_columns %in% colnames(input_file))){
    stop(paste0(paste0(setdiff(necessary_columns, colnames(input_file)), collapse = ", "),
                " are missing in univariate biomarker input file!"))
  }

  # cast the responses into a matrix
  M <- input_file %>%
    dplyr::filter(is.finite(.data[[response_column]])) %>%
    tidyr::unite(cn, tidyselect::all_of(treatment_columns), sep = "::") %>%
    reshape2::acast(depmap_id ~ cn, value.var = response_column, fun.aggregate = aggregate_function) %>%
    transform_function()

  #TODO: Decide on agg function

  # generate the biomarker table
  univariate_biomarker_table <- univariate_biomarker_table(Y = M, features = features,file = depmap_file, n.X.min = min_sample_size,
                                         v.X.min = min_x_variance, n_stable.min = n_stable_min, q_val_max = q_val_max, 
                                         regression_coef = regression_coef, stability_score = stability_score, rank.max = rank.max)

  # Export the biomarker table -----
  print(paste0("Writing the univariate output file to ", paste0(out_path, "/", output_file_name)))
  input_file %>%
    dplyr::select(tidyselect::all_of(treatment_columns)) %>%
    dplyr::distinct() %>%
    tidyr::unite(y, all_of(treatment_columns), sep = "::", remove = FALSE) %>%
    dplyr::inner_join(univariate_biomarker_table) %>%
    dplyr::select(-y) %>%
    write_out_table(path = file.path(out_path, output_file_name))
}
