

#' Fits Random Forest models for biomarker prediction using univariate target recovery.
#'
#' @param X : Feature matrix (n samples by p features) with row names as DepMap IDs.
#' @param Y : Response matrix (n samples by m compounds) with row names as DepMap IDs.
#' @param biomarker_file : Path to the downloaded depmap_datasets.h5 file.
#' @param CompoundList : Dataframe containing compound annotations (must include cn, CompoundName, and GeneSymbolOfTargets).
#' @param test_samples : Optional vector of sample IDs (DepMap IDs) to hold out for testing. Default is NULL.
#' @param bm_th : Minimum r-squared correlation threshold for feature selection, default is 0.05.
#' @param bm_R : Rank threshold for initial feature selection, default is 10.
#' @param bm_R2 : Final rank threshold after q-value sorting, default is 50.
#' @param features : Character vector of feature sets to include.
#'
#' @return A list containing four dataframes: model performances, predictions, variable importances, and the baseline univariate biomarkers table.
#' @export
biomarker_suite_rf <- function(X, Y, biomarker_file, CompoundList, test_samples = NULL, bm_th = 0.05, bm_R = 10,  bm_R2 = 50, features = c("CRISPR", "RNAi", "Expression", "Mutation", "CopyNumber", "Fusion", "Lineage")) {
  require(tidyverse)
  require(ranger)
  
  
  # Fast set operations to define train/test splits
  train <- intersect(setdiff(rownames(Y), test_samples), rownames(X))
  test <- intersect(intersect(rownames(Y), test_samples), rownames(X))
  
  # Compute target recovery biomarkers on the training set
  bm.auc <- target_recovery(Y[train, , drop = FALSE], biomarker_file, CompoundList, features = features)
  
  CompoundList <- CompoundList %>%
    dplyr::rename(cn = SampleID) 
  
  # Clean up targets mapping
  targets <- CompoundList %>% 
    dplyr::distinct(cn, CompoundName, GeneSymbolOfTargets) %>% 
    tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>% 
    dplyr::mutate(GeneSymbolOfTargets = trimws(GeneSymbolOfTargets)) %>% 
    dplyr::distinct() %>%
    tidyr::drop_na()
  
  # Map X column names to actual gene symbols
  feat_map <- tibble::tibble(cn = colnames(X)) %>%  
    dplyr::mutate(
      t1 = stringr::word(stringr::word(cn, 2, sep = stringr::fixed("_")), sep = stringr::fixed(".")),
      t2 = stringr::word(stringr::word(cn, 2, sep = stringr::fixed("_")), sep = stringr::fixed("--")),
      t3 = stringr::word(stringr::word(cn, 2, sep = stringr::fixed("_")), -1, sep = stringr::fixed("--"))
    ) %>% 
    tidyr::pivot_longer(cols = c(t1, t2, t3), values_to = "GeneSymbolOfTargets", names_to = "d") %>% 
    dplyr::filter(GeneSymbolOfTargets != "X", !is.na(GeneSymbolOfTargets)) %>%
    dplyr::select(-d) %>% 
    dplyr::distinct() %>% 
    dplyr::rename(cn.feat = cn) %>% 
    dplyr::inner_join(targets, by = "GeneSymbolOfTargets")
  
  # Filter for selected features based on thresholds
  selected_features <- bm.auc %>% 
    dplyr::distinct(cn, CompoundName, feature_set, feature, rank, correlation_coef, q_val, status) %>%
    dplyr::filter((status == "Other") | (feature_set == "Lineage")) %>% 
    dplyr::filter(((rank <= bm_R) & (correlation_coef^2 >= bm_th)) | (rank <= 1)) 
  
  if (nrow(selected_features) > 0) {
    selected_features <- selected_features %>%
      dplyr::group_by(cn) %>% 
      dplyr::mutate(rank_ = dplyr::min_rank(q_val)) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(rank_ <= bm_R2) 
  }
  
  # Named vector dictionary to replace the nested ifelse statement
  prefix_map <- c(
    "CRISPR"     = "CRISPR_", 
    "RNAi"       = "RNAi_", 
    "Expression" = "EXP_", 
    "CopyNumber" = "CN_", 
    "Mutation"   = "MUT_", 
    "Fusion"     = "FUS_", 
    "Lineage"    = "LIN_"
  )
  
  selected_columns <- selected_features %>%
    dplyr::distinct(cn, CompoundName, feature_set, feature) %>% 
    dplyr::mutate(
      tar = prefix_map[feature_set],
      y = paste0(tar, feature)
    ) %>% 
    dplyr::filter(!is.na(tar)) 
  
  # Inner model fitting function
  fit <- function(x, y) {
    cl <- intersect(train, names(y))
    cl.test <- intersect(test, names(y))
    all_cl <- union(cl, cl.test)
    
    rf <- ranger::ranger(x = x[cl, , drop = FALSE], y = y[cl], importance = "impurity")
    pr <- predict(rf, data = x[all_cl, , drop = FALSE])
    
    y.hat <- tibble::tibble(
      depmap_id = all_cl, 
      y.hat = pr$predictions, 
      y = y[all_cl],
      type = ifelse(depmap_id %in% cl, "train", "test")
    ) %>% 
      dplyr::left_join(tibble::tibble(depmap_id = cl, y.hat.oob = rf$predictions), by = "depmap_id")
    
    imp <- tibble::tibble(var = names(rf$variable.importance), imp = rf$variable.importance) %>%
      dplyr::arrange(desc(imp))
    
    res <- tibble::tibble(
      mse.oob = mean((rf$predictions - y[cl])^2, na.rm = TRUE),
      var.y.train = var(y[cl], na.rm = TRUE), 
      r2.oob = mse.oob / var.y.train,
      r.oob = cor(rf$predictions, y[cl], use = "p")[, 1]
    )
    
    if (length(cl.test) > 0) {
      test_res <- y.hat %>% 
        dplyr::filter(type == "test") %>% 
        dplyr::summarise(
          var.y.test = var(y, na.rm = TRUE),
          mse = mean((y - y.hat)^2, na.rm = TRUE),
          r2 = 1 - mse / var.y.test,
          r = cor(y, y.hat, use = "p")[, 1]
        )
      res <- dplyr::bind_cols(test_res, res)
    }
    
    return(list(res, y.hat, imp))
  }
  
  
  # Initialize lists for accumulating RF models
  biomarker_table <- list()
  prediction_table <- list()
  importance_table <- list()
  jx <- 1
  
  for (cmp in colnames(Y)) {
    tars <- feat_map %>% dplyr::filter(cn == cmp) %>% dplyr::pull(GeneSymbolOfTargets) %>% unique()
    extras <- selected_columns %>% dplyr::filter(cn == cmp) %>% dplyr::pull(y) %>% intersect(colnames(X))
    
    y <- Y[, cmp]
    y <- y[is.finite(y)]
    
    res_list <- list()
    pred_list <- list()
    imp_list <- list()
    ix <- 1 
    
    if (length(tars) > 0) {
      # Fit a model for each individual target
      for (tar in tars) {
        tar_feats <- feat_map %>% dplyr::filter(cn == cmp, GeneSymbolOfTargets == tar) %>% dplyr::pull(cn.feat) %>% unique()
        temp <- fit(x = X[names(y), tar_feats, drop = FALSE], y = y)
        
        res_list[[ix]]  <- temp[[1]] %>% dplyr::mutate(model = tar)
        pred_list[[ix]] <- temp[[2]] %>% dplyr::mutate(model = tar)
        imp_list[[ix]]  <- temp[[3]] %>% dplyr::mutate(model = tar)
        ix <- ix + 1
      }
      
      # Fit a model for ALL targets together
      all_tar_feats <- feat_map %>% dplyr::filter(cn == cmp, GeneSymbolOfTargets %in% tars) %>% dplyr::pull(cn.feat) %>% unique()
      temp <- fit(x = X[names(y), all_tar_feats, drop = FALSE], y = y)
      
      res_list[[ix]]  <- temp[[1]] %>% dplyr::mutate(model = "targets")
      pred_list[[ix]] <- temp[[2]] %>% dplyr::mutate(model = "targets")
      imp_list[[ix]]  <- temp[[3]] %>% dplyr::mutate(model = "targets")
      ix <- ix + 1
    }
    
    cols <- unique(union(feat_map %>% dplyr::filter(cn == cmp, GeneSymbolOfTargets %in% tars) %>% dplyr::pull(cn.feat), extras))
    if (length(cols) > 0) {
      # Fit the extended model (Targets + Extracted Predictors)
      temp <- fit(x = X[names(y), cols, drop = FALSE], y = y)
      
      res_list[[ix]]  <- temp[[1]] %>% dplyr::mutate(model = "extended")
      pred_list[[ix]] <- temp[[2]] %>% dplyr::mutate(model = "extended")
      imp_list[[ix]]  <- temp[[3]] %>% dplyr::mutate(model = "extended")
    }
    
    # Bundle compound results together
    biomarker_table[[jx]]  <- dplyr::bind_rows(res_list) %>% dplyr::mutate(cn = cmp)
    prediction_table[[jx]] <- dplyr::bind_rows(pred_list) %>% dplyr::mutate(cn = cmp)
    importance_table[[jx]] <- dplyr::bind_rows(imp_list) %>% dplyr::mutate(cn = cmp)
    
    message(paste0(cmp, " - ", jx))
    jx <- jx + 1
  }
  
  cmp_annotations <- CompoundList %>% dplyr::distinct(cn, CompoundName)
  
  return(list(
    dplyr::bind_rows(biomarker_table) %>% dplyr::left_join(cmp_annotations, by = "cn"), 
    dplyr::bind_rows(prediction_table) %>% dplyr::left_join(cmp_annotations, by = "cn"), 
    dplyr::bind_rows(importance_table) %>% dplyr::left_join(cmp_annotations, by = "cn"), 
    bm.auc
  ))
}


#' Performs K-fold cross-validation for the Random Forest biomarker suite.
#'
#' @param X : Feature matrix (n samples by p features) with row names as DepMap IDs.
#' @param Y : Response matrix (n samples by m compounds) with row names as DepMap IDs.
#' @param biomarker_file : Path to the downloaded depmap_datasets.h5 file.
#' @param CompoundList : Dataframe containing compound annotations.
#' @param bm_th : Minimum r-squared correlation threshold for feature selection, default is 0.05.
#' @param bm_R : Rank threshold for initial feature selection, default is 10.
#' @param bm_R2 : Final rank threshold after q-value sorting, default is 50.
#' @param K : Number of cross-validation folds, default is 10.
#' @param seed : Random seed for reproducible fold generation, default is NULL.
#'
#' @return A list containing combined model_performances, predictions, variable_importances (across all K folds + full model where K=0), and univariate_biomarkers.
#' @export
biomarker_suite_rf_cv <- function(X, Y, biomarker_file, CompoundList, bm_th = 0.05, bm_R = 10, bm_R2 = 50, K = 10, seed = NULL) {
  require(tidyverse)
  
  if (!is.null(seed)) set.seed(seed)
  
  cl <- intersect(rownames(X), rownames(Y)) %>% sample()
  
  RES <- list()
  PRED <- list()
  IMP <- list()
  
  # Cross Validation Folds
  for (k in 1:K) {
    message("Running fold: ", k)
    
    # Calculate fold indices safely
    fold_samples <- cl[seq.int(k, length(cl), by = K)]
    
    temp <- biomarker_suite_rf(
      X, Y, 
      biomarker_file = biomarker_file, 
      CompoundList = CompoundList, 
      test_samples = fold_samples, 
      bm_th = bm_th, bm_R = bm_R, bm_R2 = bm_R2
    )
    
    RES[[k]]  <- temp[[1]] %>% dplyr::mutate(K = k)
    PRED[[k]] <- temp[[2]] %>% dplyr::mutate(K = k)
    IMP[[k]]  <- temp[[3]] %>% dplyr::mutate(K = k)
  }
  
  message("Running Full Model (K = 0)")
  
  # Full Model (No test samples)
  temp_full <- biomarker_suite_rf(
    X, Y, 
    biomarker_file = biomarker_file, 
    CompoundList = CompoundList, 
    test_samples = NULL, 
    bm_th = bm_th, bm_R = bm_R, bm_R2 = bm_R2
  )
  
  RES[[K + 1]]  <- temp_full[[1]] %>% dplyr::mutate(K = 0)
  PRED[[K + 1]] <- temp_full[[2]] %>% dplyr::mutate(K = 0)
  IMP[[K + 1]]  <- temp_full[[3]] %>% dplyr::mutate(K = 0)
  
  return(list(
    model_performances = dplyr::bind_rows(RES), 
    predictions = dplyr::bind_rows(PRED), 
    variable_importances = dplyr::bind_rows(IMP), 
    univariate_biomarkers = temp_full[[4]]
  ))
}

#' Calculates and returns univariate analysis results by correlating each column of Y with features sets from depmap.org.
#'
#' @param Y Matrix n x m, make sure rownames are depmap_id's and columns are named. There can be NA's.
#' @param file Please point out to the downloaded depmap_datasets.h5 file.
#' @param features You can give a subset of available feature sets, but if left as NULL, it computes for everything.
#' @param n.X.min Results with less than a given sample size are dropped, default is 250.
#' @param v.X.min Minimum variance for the columns of X, default is 0.0025.
#' @param q_val_max Results with q-values less than q_val_max are dropped, default is 0.1.
#' @param rank.max Results with ranks (by q-value) greater than rank.max are dropped, default is 250.
#'
#' @return Returns a data-table with each row corresponds to a particular (feature, feature_set, y) triplet.
#' @export
univariate_biomarker_table <- function(Y, file, features = NULL, n.X.min = 250, v.X.min = 0.0025, q_val_max = 0.1, rank.max = 250, LM = NULL) {
  require(tidyverse)
  require(rhdf5)
  require(WGCNA)
  
  
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  
  # Fetch available features once to avoid redundant h5ls calls
  h5_groups <- rhdf5::h5ls(file)
  available_features <- substr(unique(h5_groups$group[h5_groups$name == "mat"]), 2, 100)
  
  if (is.null(features)) {
    features <- available_features
  } else {
    features <- intersect(features, available_features)
  }
  
  message("Features to process: ", paste(features, collapse = ", "))
  
  RESULTS <- list()
  for (feat in features) {
    
    if(feat == "Lineage" & !is.null(LM)){
      X <- LM
    }else{
      X <- read_dataset(file, feat)
    }
    cl <- intersect(rownames(X), rownames(Y))
    X <- X[cl, , drop = FALSE]
    
    if (nrow(X) >= n.X.min && ncol(X) > 0) {
      message(paste0("Processing ", feat, " - ", nrow(X), 'x', ncol(X)))
      
      RESULTS[[feat]] <- linear_model(
        X = X, 
        Y = Y[cl, , drop = FALSE],
        v.X.min = v.X.min, 
        n.min = n.X.min
      ) %>% 
        dplyr::filter(rank <= rank.max, q_val <= q_val_max) %>% 
        dplyr::rename(feature = x) %>%
        dplyr::mutate(feature_set = feat) 
    }
  }
  
  return(dplyr::bind_rows(RESULTS))
}


#' Fits a simple linear model by regressing y over each column of X.
#'
#' @param X Matrix of n by m, it can have NA's in it, columns should be named.
#' @param Y Matrix of n by p, rows are in the same order of with the rows of X.
#' @param v.X.min Minimum variance for the columns of X, columns with smaller variances will be dropped.
#' @param n.min Minimum number of finite pairs between columns of X and Y, column pairs not satisfying this condition will be dropped.
#'
#' @return A dataframe containing the correlation coefficient, p-value, q-value, rank, and variable variances for each valid column pair.
#' @export
linear_model <- function(X, Y, v.X.min = 0.0025, n.min = 100) {
  require(tidyverse)
  require(matrixStats)
  require(WGCNA)
  
  # Extremely fast matrix-to-dataframe conversion (Replaces reshape2::melt)
  cor_res <- WGCNA::corAndPvalue(X, Y, use = "p")
  
  cor.table <- tibble::tibble(
    x = rep(rownames(cor_res$cor), times = ncol(cor_res$cor)),
    y = rep(colnames(cor_res$cor), each = nrow(cor_res$cor)),
    correlation_coef = as.vector(cor_res$cor),
    p_val = as.vector(cor_res$p),
    n = as.vector(cor_res$nObs)
  ) %>% 
    dplyr::filter(n >= n.min)
  
  
  # Variance map for Y
  masks_Y <- is.finite(Y)
  masks_Y <- masks_Y[, !duplicated(t(masks_Y)), drop = FALSE]
  colnames(masks_Y) <- paste0("m", seq_len(ncol(masks_Y)))
  
  map_Y <- apply(masks_Y, 2, function(m) apply(is.finite(Y) == m, 2, all)) %>% apply(1, which.max)
  map_Y_df <- tibble::tibble(y = names(map_Y), m = colnames(masks_Y)[map_Y])
  
  masks_Y[!masks_Y] <- NA
  vX <- apply(masks_Y, 2, function(m) matrixStats::colVars(X * m, na.rm = TRUE))
  vX <- vX[, map_Y_df$m, drop = FALSE]
  colnames(vX) <- map_Y_df$y
  
  # Variance map for X
  masks_X <- is.finite(X)
  masks_X <- masks_X[, !duplicated(t(masks_X)), drop = FALSE]
  colnames(masks_X) <- paste0("m", seq_len(ncol(masks_X)))
  
  map_X <- apply(masks_X, 2, function(m) apply(is.finite(X) == m, 2, all)) %>% apply(1, which.max)
  map_X_df <- tibble::tibble(x = names(map_X), m = colnames(masks_X)[map_X])
  
  masks_X[!masks_X] <- NA
  vY <- apply(masks_X, 2, function(m) matrixStats::colVars(Y * m, na.rm = TRUE))
  vY <- vY[, map_X_df$m, drop = FALSE]
  colnames(vY) <- map_X_df$x
  
  # Helper to melt variance matrices cleanly
  melt_var <- function(mat, val_name, col1, col2) {
    as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
      stats::setNames(c(col1, col2, val_name))
  }
  
  cor.table <- cor.table %>% 
    dplyr::left_join(melt_var(vX, "var.x", "x", "y"), by = c("x", "y")) %>% 
    dplyr::left_join(melt_var(vY, "var.y", "y", "x"), by = c("x", "y")) %>% 
    dplyr::mutate(regression_coef = correlation_coef * sqrt(var.y) / sqrt(var.x)) %>% 
    dplyr::filter(n >= n.min, var.x >= v.X.min)
  
  if (nrow(cor.table) > 0) {
    cor.table <- cor.table %>% 
      dplyr::group_by(y) %>% 
      dplyr::mutate(
        q_val = p.adjust(p_val, method = "BH"),
        rank = dplyr::min_rank(q_val) 
      ) %>% 
      dplyr::ungroup() %>%
      dplyr::select(x, y, correlation_coef, regression_coef, q_val, rank, p_val, n, var.x, var.y)
  }
  
  return(cor.table) 
}


#' Exports individual datasets from depmap_datasets.h5 file.
#'
#' @param file Please point out to the downloaded depmap_datasets.h5 file.
#' @param dataset The dataset you want to export.
#' @param rownames_depmap_ids Default TRUE, you can get the rownames as ccle_names by switching to FALSE.
#'
#' @return The requested data matrix with cleaned row and column names.
#' @export
read_dataset <- function(file = '/data/biomarker/current/depmap_datasets_public.h5', dataset, rownames_depmap_ids = TRUE) {
  require(rhdf5)
  
  # Cleaner check for S3 vs Local
  is_s3 <- grepl("^(s3|http|https)://", file)
  message("Reading ", file, if(is_s3) " from S3" else " from local")
  
  X <- h5read(file, name = paste0(dataset, "/mat"), s3 = is_s3)
  row_meta <- h5read(file, name = paste0(dataset, "/row_meta"), s3 = is_s3)
  column_meta <- h5read(file, name = paste0(dataset, "/column_meta"), s3 = is_s3)
  
  colnames(X) <- column_meta$column_name
  rownames(X) <- if(rownames_depmap_ids) row_meta$ModelID else row_meta$CCLEName
  
  # Filter missing dimensions
  X <- X[rownames(X) != "NA", colnames(X) != "NA", drop = FALSE]
  X <- X[!duplicated(rownames(X)), !duplicated(colnames(X)), drop = FALSE]
  
  return(X)
}


#' Extracts target recovery biomarkers and classifies their correlations.
#'
#' @param Y Matrix n x m, make sure rownames are depmap_id's and columns are named.
#' @param file Path to the downloaded depmap_datasets.h5 file.
#' @param compound_annotations Dataframe containing prior knowledge target mappings.
#' @param features Character vector of feature sets to evaluate.
#' @param rank.max Maximum rank limit for univariate feature inclusion.
#' @param q.max Maximum q-value limit for univariate feature inclusion.
#' @param n.min Minimum number of overlapping samples.
#'
#' @return A dataframe with features classified into Target, Target-Correlate, Lineage-Correlate, or Other.
#' @export
target_recovery <- function(Y, file, compound_annotations, features = c("CRISPR", "RNAi", "Expression", "Mutation", "CopyNumber", "Fusion", "Lineage", "Repurposing.Primary"), 
                            rank.max = 100, q.max = 0.1, n.min = 250, LM = NULL) {
  require(tidyverse)
  
  compound_annotations <- compound_annotations %>%
    dplyr::rename(cn = SampleID) 
  
  
  # ---------------------------------------------------------
  # HELPER: Calculate Partial Correlation
  # ---------------------------------------------------------
  calc_partial_cor <- function(cor_xy, cor_xz, cor_yz) {
    numerator <- cor_xy - (cor_xz * cor_yz)
    denominator <- sqrt(1 - cor_yz^2) * sqrt(1 - cor_xz^2)
    return(numerator / denominator)
  }
  
  # ---------------------------------------------------------
  # PHASE 1: Compute Initial Univariate Biomarkers
  # ---------------------------------------------------------
  bm <- univariate_biomarker_table(
    Y, file, features = features, 
    rank.max = rank.max, q_val_max = q.max, n.X.min = n.min, LM = LM
  )
  
  if (nrow(bm) > 0) {
    bm <- bm %>% 
      dplyr::mutate(cn = stringr::word(y, 1, sep = stringr::fixed("::")))
  }
  
  # ---------------------------------------------------------
  # PHASE 2: Map Features to Gene Symbols
  # ---------------------------------------------------------
  
  
  feature_gene_map <- bm %>%
    dplyr::distinct(feature) %>%
    dplyr::mutate(
      feature2 = stringr::word(feature, sep = stringr::fixed(".")),
      feature3 = stringr::word(feature2, sep = stringr::fixed("--")),
      feature4 = stringr::word(feature2, 2, sep = stringr::fixed("--"))
    ) %>%
    tidyr::pivot_longer(
      cols = c(feature2, feature3, feature4),
      names_to = "dummy",
      values_to = "GeneSymbolOfTargets"
    ) %>%
    dplyr::select(-dummy) %>%
    dplyr::mutate(GeneSymbolOfTargets = trimws(GeneSymbolOfTargets)) %>% 
    tidyr::drop_na() %>%
    dplyr::filter(GeneSymbolOfTargets != "", GeneSymbolOfTargets != "NA", GeneSymbolOfTargets != "X") %>% # <-- Strict filter
    dplyr::distinct()
  
  
  # ---------------------------------------------------------
  # PHASE 3: Annotate Known Drug Targets
  # ---------------------------------------------------------
  bm <- compound_annotations %>%
    dplyr::distinct(cn, GeneSymbolOfTargets) %>%
    tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>%
    dplyr::mutate(GeneSymbolOfTargets = trimws(GeneSymbolOfTargets)) %>% 
    tidyr::drop_na() %>% 
    dplyr::filter(GeneSymbolOfTargets != "", GeneSymbolOfTargets != "NA", GeneSymbolOfTargets != "X") %>% # <-- Strict filter
    dplyr::inner_join(feature_gene_map, by = "GeneSymbolOfTargets") %>% 
    dplyr::select(-GeneSymbolOfTargets) %>%
    dplyr::distinct() %>% 
    dplyr::mutate(is.target = TRUE) %>%
    dplyr::right_join(bm, by = c("cn", "feature"))
  
  
  
  # Distinct list of features to process
  feats_to_process <- bm %>% dplyr::distinct(feature, feature_set, is.target, cn)
  
  # ---------------------------------------------------------
  # PHASE 4: Load Lineage Matrix
  # ---------------------------------------------------------
  if(is.null(LM)){
    mat_lineage <- read_dataset(file, "Lineage")
  }else{
    mat_lineage <- LM
  }

  cl_common <- intersect(rownames(mat_lineage), rownames(Y))
  mat_lineage <- mat_lineage[cl_common, , drop = FALSE]
  
  
  # ---------------------------------------------------------
  # PHASE 5: Loop Over Each Feature Set to Calculate Partial Correlations
  # ---------------------------------------------------------
  results_list <- list()
  
  for (current_feat_set in unique(feats_to_process$feature_set)) {
    
    # Load current feature matrix
    if(current_feat_set == "Lineage" & !is.null(LM)){
      mat_current <- LM
    }else{
      mat_current <- read_dataset(file, current_feat_set)
    }
    cl_current <- intersect(rownames(mat_current), rownames(Y))
    mat_current <- mat_current[cl_current, , drop = FALSE]
    
    # Separate known targets from current features
    target_feats_df <- feats_to_process %>% dplyr::filter(is.target) %>% dplyr::distinct(feature_set, feature)
    current_feats_vec <- feats_to_process %>% dplyr::filter(feature_set == current_feat_set) %>% dplyr::pull(feature) %>% unique() %>% as.character()
    
    # Filter base biomarker table for the current loop
    current_results <- bm %>% dplyr::filter(feature_set == current_feat_set)
    
    # --- 5A: Account for Target Confounding ---
    if (nrow(target_feats_df) > 0) {
      target_cors_list <- list()
      
      # Calculate correlation between target features and current features
      for (tgt_feat_set in unique(target_feats_df$feature_set)) {
        mat_target <- read_dataset(file, tgt_feat_set)
        cl_tgt <- intersect(cl_current, rownames(mat_target))
        
        tgt_feats_vec <- target_feats_df %>% dplyr::filter(feature_set == tgt_feat_set) %>% dplyr::pull(feature) %>% unique() %>% as.character()
        
        cor_mat <- WGCNA::cor(mat_target[cl_tgt, tgt_feats_vec, drop = FALSE], 
                              mat_current[cl_tgt, current_feats_vec, drop = FALSE], use = "p")
        
        target_cors_list[[tgt_feat_set]] <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE) %>%
          stats::setNames(c("target", "feature", "feature_target.cor")) %>%
          dplyr::mutate(
            feature_set_target = tgt_feat_set,
            target = as.character(target),
            feature = as.character(feature)
          ) %>% 
          tibble::as_tibble()
      }
      
      
      # Merge target correlations into the results
      tars2 <- bm %>%
        dplyr::filter(is.target) %>%
        dplyr::rename(target = feature, target.cor = correlation_coef, feature_set_target = feature_set) %>%
        dplyr::left_join(dplyr::bind_rows(target_cors_list), by = c("target", "feature_set_target"), relationship = "many-to-many") %>%
        dplyr::distinct(cn, target, feature, target.cor, feature_target.cor, feature_set_target) %>%
        dplyr::mutate(feature_set = current_feat_set)
      
      current_results <- current_results %>%
        dplyr::left_join(tars2, by = c("cn", "feature", "feature_set")) %>%
        dplyr::mutate(
          target.cor = tidyr::replace_na(target.cor, 0),
          target_part_cor_coef = calc_partial_cor(correlation_coef, feature_target.cor, target.cor)
        ) %>%
        dplyr::group_by(cn, y, feature, feature_set) %>%
        dplyr::arrange(target_part_cor_coef^2 * (1 - target.cor^2), target) %>%
        dplyr::slice(1) %>% 
        dplyr::ungroup()
    }
    
    # --- 5B: Account for Lineage Confounding ---
    cl_lineage <- intersect(rownames(mat_current), rownames(mat_lineage))
    l_cor_mat <- WGCNA::cor(mat_lineage[cl_lineage, , drop = FALSE], 
                            mat_current[cl_lineage, current_feats_vec, drop = FALSE], use = "p")
    
    lineage_cor_df <- as.data.frame(as.table(l_cor_mat), stringsAsFactors = FALSE) %>%
      stats::setNames(c("lineage", "feature", "feature_lineage.cor")) %>%
      dplyr::mutate(lineage = as.character(lineage), feature = as.character(feature)) %>% 
      tibble::as_tibble()
    
    tars3 <- bm %>%
      dplyr::filter(feature_set == "Lineage") %>%
      dplyr::rename(lineage = feature, lineage.cor = correlation_coef) %>%
      dplyr::left_join(lineage_cor_df, by = "lineage", relationship = "many-to-many") %>%
      dplyr::distinct(cn, lineage, feature, lineage.cor, feature_lineage.cor) %>%
      dplyr::mutate(feature_set = current_feat_set)
    
    # Calculate lineage partial correlation and finalize loop results
    current_results <- current_results %>%
      dplyr::left_join(tars3, by = c("cn", "feature", "feature_set")) %>%
      dplyr::mutate(
        lineage_part_cor_coef = calc_partial_cor(correlation_coef, feature_lineage.cor, lineage.cor)
      ) %>%
      dplyr::group_by(cn, y, feature, feature_set) %>%
      dplyr::arrange(lineage_part_cor_coef^2 * (1 - lineage.cor^2), lineage) %>%
      dplyr::slice(1) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(is.target = !is.na(is.target))
    
    results_list[[current_feat_set]] <- current_results
  }
  
  
  # ---------------------------------------------------------
  # PHASE 6: Final Status Classification
  # ---------------------------------------------------------

  
  final_res <- dplyr::bind_rows(results_list) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(
      target_part_cor_coef  = tidyr::replace_na(target_part_cor_coef, correlation_coef),
      target.cor            = tidyr::replace_na(target.cor, 0),
      lineage_part_cor_coef = tidyr::replace_na(lineage_part_cor_coef, correlation_coef),
      lineage.cor           = tidyr::replace_na(lineage.cor, 0),
      
      # Calculate variance fractions
      target_var_fraction   = (1 - target.cor^2) * target_part_cor_coef^2 / correlation_coef^2,
      lineage_var_fraction  = (1 - lineage.cor^2) * lineage_part_cor_coef^2 / correlation_coef^2
    ) %>%
    dplyr::mutate(
      status = dplyr::case_when(
        is.target ~ "Target",
        target_var_fraction < 0.5 ~ "Target-Correlate",
        lineage_var_fraction < 0.5 ~ "Lineage-Correlate",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::left_join(compound_annotations, by = "cn")
  
  return(final_res)
}



