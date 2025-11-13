##################################################
## R script for ExpressAnalyst
## Description: functions for DE analysis
##
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#' Select Metadata Factors for DE Analysis
#'
#' This function is used to select metadata factors for a Differential Expression (DE) analysis.
#'
#' @param dataName A character string specifying the name of the dataset.
#' @param meta0 The primary metadata factor
#' @param meta1 The secondary metadata factor
#' @param block1 The blocking factor for
#'
#' @author Guangyan Zhou, \email{guangyan.zhou@mail.mcgill.ca}
#' @details Additional details about the function, if needed.
#'
#' @export
#'
SetSelectedMetaInfo <- function(dataName="", meta0, meta1, block1){
  dataSet <- readDataset(dataName);
  if(meta0 == "NA"){
    RegisterData(dataSet, 0);
  }else{
    rmidx <- which(dataSet$meta.info.[, meta0]=="NA")
    if(meta1 != "NA"){
        rmidx <- c(rmidx,which(dataSet$meta.info[, meta1]=="NA"))
    }
    if(length(rmidx)>0){
       meta<- dataSet$meta.info[-rmidx,]
     #print(meta);
       for(col in 1:ncol(meta)){
        meta[,col]<- droplevels(meta[,col])
       }
       dataSet$rmidx <- rmidx
    }else{
        meta<- dataSet$meta.info
    }
    cls <- meta[, meta0];
    dataSet$fst.cls <- cls; # for PCA plotting
    block <- NULL;
    dataSet$sec.cls <- "NA";
    if(meta1 != "NA"){
      if(block1){
        block <- meta[, meta1];
      }else{ # two factor
        cls <- interaction(meta[, c(meta0, meta1)], sep = "_", lex.order = TRUE);
      }
      dataSet$sec.cls <- meta[, meta1]; # for pca coloring
    }
    dataSet$analysisVar <- meta0 
    dataSet$secondVar <- meta1
    dataSet$cls <- cls; # record main cls;
    dataSet$block <- block;
    RegisterData(dataSet, levels(cls)[levels(cls)!="NA"]);
  }
}

#' Perform Differential Analysis
#'
#' This function performs differential analysis based on different types of comparisons.
#'
#' @param dataName A character string specifying the name of the dataset.
#' @param anal.type The type of analysis to perform. Options: "default", "custom", "time", "reference", "nested".
#' @param par1 Parameter 1, depending on the analysis type.
#' @param par2 Parameter 2, depending on the analysis type.
#' @param nested.opt Option for nested analysis. Options: "intonly" (default), "all".
#' @param robustTrend Logical. If TRUE, use robust trend test.
#'
#' @return Results of the differential analysis.
#' @details default: all pair-wise comparison (A-B) + (B-C) + (A-C), custom: only compare two groups (A-C), time: compare consecutive groups (B-A) + (C-B), reference: all others against common reference (A-C) + (B-C), nested: (A-B)+(C-D) 
#' @author Guangyan Zhou, \email{guangyan.zhou@mail.mcgill.ca}
#' @export
#' @license MIT
#'
PerformDEAnal<-function (dataName="", anal.type = "default", par1 = NULL, par2 = NULL, nested.opt = "intonly", robustTrend=F){

  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  if (dataSet$de.method == "deseq2") {
    dataSet <- prepareContrast(dataSet, anal.type, par1, par2, nested.opt);
    .prepare.deseq(dataSet, anal.type, par1, par2 , nested.opt);
    .perform.computing();
    dataSet <- .save.deseq.res(dataSet);
  }else if (dataSet$de.method == "limma"){
    dataSet <- prepareContrast(dataSet, anal.type, par1, par2, nested.opt);
    dataSet <- .perform_limma_edger(dataSet, robustTrend);
  }else if (dataSet$de.method == "edger"){
    dataSet <- prepareEdgeRContrast(dataSet, anal.type, par1, par2, nested.opt);
    dataSet <- .perform_limma_edger(dataSet, robustTrend);
  }else{ #dataSet$de.method == "wtt"
    dataSet <- prepareContrast(dataSet, anal.type, par1, par2, nested.opt);
    dataSet <- .perform_williams_trend(dataSet, robustTrend);
  }
  return(RegisterData(dataSet));
}

.prepare.deseq <- function(dataSet, anal.type, par1, par2, nested.opt) {

  my.fun <- function() {
    require(DESeq2)

    # Helper: prefix numeric-looking group labels
    formatLevel <- function(x) {
      if (grepl("^[0-9]", x)) paste0(dataSet$analysisVar, "_", x) else x
    }

    # Helper: parse contrast string "A vs. B"
    parse_contrast_groups <- function(cstr) {
      comps <- strsplit(cstr, " vs\\.?\\s*")[[1]]
      if (length(comps) != 2) stop(paste("Invalid contrast format:", cstr))
      return(comps)
    }

    # Extract count data
    data.anot <- .get.annotated.data();
    if (length(dataSet$rmidx) > 0) {
      data.anot <- data.anot[, -dataSet$rmidx]
    } 

    # Format class labels
    if (any(grepl("(^[0-9]+).*", dataSet$fst.cls))) {
      fst.cls <- paste0(dataSet$analysisVar, "_", dataSet$fst.cls)
    } else {
      fst.cls <- dataSet$fst.cls
    }
    fst.cls <- as.character(fst.cls)                   # <<< NEW

    all_conditions <- unique(fst.cls)
    contrast_list <- list()

    # ---- Single-factor designs ----
    colData <- data.frame(condition = factor(fst.cls,     # <<< REPLACED
                                             levels = all_conditions))

    if (!is.null(dataSet$block)) {
      colData$block <- factor(dataSet$block)
      design <- ~ block + condition
    } else {
      design <- ~ condition
    }

    if (anal.type == "default") {
      for (i in 1:(length(all_conditions) - 1))
        for (j in (i + 1):length(all_conditions)) {
          contrast_name <- paste0(all_conditions[i], " vs ", all_conditions[j])
          contrast_list[[contrast_name]] <-
            c("condition",
              all_conditions[j],   # NUMERATOR = second term
              all_conditions[i])   # DENOMINATOR = first term
        }

    } else if (anal.type == "reference") {
      ref <- formatLevel(par1)
      if (!(ref %in% all_conditions))
        stop("Reference level not found: ", ref)

      for (cond in setdiff(all_conditions, ref)) {
        contrast_name <- paste0(ref, " vs ", cond)
        contrast_list[[contrast_name]] <-
          c("condition", cond, ref)   # numerator = second term
      }

    } else if (anal.type == "custom") {
      comps <- parse_contrast_groups(par1)
      comps <- vapply(comps, formatLevel, "")
      if (!all(comps %in% all_conditions))
        stop("Invalid custom contrast: ", par1)

      contrast_name <- paste0(comps[1], " vs ", comps[2])
      contrast_list[[contrast_name]] <-
        c("condition", comps[2], comps[1]) # numerator = second term
    }

    # ---- Run DESeq2 ----
    dds <- DESeqDataSetFromMatrix(countData = round(data.anot),
                                  colData   = colData,
                                  design    = design)
    dds <- DESeq(dds, betaPrior = FALSE)
    qs::qsave(dds, "deseq.res.obj.rds")

    # ---- Extract contrast results ----
    results_list <- list()
    if (length(contrast_list) > 0) {
      for (contrast_name in names(contrast_list)) {
        res <- results(dds,
                       contrast            = contrast_list[[contrast_name]],
                       independentFiltering = FALSE,
                       cooksCutoff          = Inf)

        topFeatures <- data.frame(res@listData)
        rownames(topFeatures) <- rownames(res)
        colnames(topFeatures) <- sub("padj", "adj.P.Val",  colnames(topFeatures))
        colnames(topFeatures) <- sub("pvalue", "P.Value",  colnames(topFeatures))
        colnames(topFeatures) <- sub("log2FoldChange","logFC",colnames(topFeatures))
        topFeatures <- topFeatures[c("logFC","baseMean","lfcSE",
                                     "stat","P.Value","adj.P.Val")]
        topFeatures <- topFeatures[order(topFeatures$P.Value), ]

        results_list[[contrast_name]] <- topFeatures
      }
    } else {
      results_list[[1]] <- .get.interaction.results()
    }

    return(results_list)
  }

  dat.in <- list(data = dataSet, my.fun = my.fun)
  qs::qsave(dat.in, file = "dat.in.qs")
  return(1)
}




.save.deseq.res <- function(dataSet){
  dat.in <- qs::qread("dat.in.qs"); 
  my.res <- dat.in$my.res;
  dataSet$comp.res.list <- my.res;
  dataSet$comp.res <- my.res[[1]];
  qs::qsave(my.res, file="dat.comp.res.qs");
  return(dataSet);
}


prepareContrast <-function(dataSet, anal.type = "reference", par1 = NULL, par2 = NULL, nested.opt = "intonly"){
  
  msgSet <- readSet(msgSet, "msgSet");
  cat(anal.type, par1, par2, nested.opt, "\n")
  set.seed(1337);
  myargs <- list();
  cls <- dataSet$cls;
  dataSet$comp.type <- anal.type;
  grp.nms <- levels(cls);
  analysisVar <- dataSet$analysisVar
  if(dataSet$cont.inx[analysisVar] |  any(grepl("(^[0-9]+).*", grp.nms))){
    if(grepl( "vs",par1)){
      par1 <- strsplit(par1, " vs. ")[[1]]
      par1 <- paste0(analysisVar,"_",par1[1]," vs. ",analysisVar,"_",par1[2])
    }else{
      par1<- paste0(analysisVar,"_",par1)
    }
    if(par2 != "NA"){
      if(grepl( "vs",par2)){
        par2 <- strsplit(par2, " vs. ")[[1]]
        par2 <- paste0(analysisVar,"_",par2[1]," vs. ",analysisVar,"_",par2[2])
      }else{
        par2<- paste0(analysisVar,"_",par2)
      }
    }

    if(any(grepl("(^[0-9]+).*",  colnames(dataSet$design)))){
      colnames(dataSet$design) = as.character(sapply( colnames(dataSet$design),function(x) paste0(analysisVar,"_",x)))
    }
    grp.nms <- paste0(analysisVar,"_",grp.nms)
    
  }

  dataSet$par1 <- par1;
  
  if (anal.type == "default") {
    inx <- 0
    for (m in 1:(length(grp.nms) - 1)) {
      for (n in (m + 1):length(grp.nms)) {
        inx <- inx + 1
        myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep = "");
      }
    }
    filename <- "sig_genes_pairwise";
  } else if (anal.type == "time") {
    for (i in 2:length(grp.nms)) {
      myargs[[i - 1]] <- paste(grp.nms[i], "-", grp.nms[i-1], sep = "")
    }
    filename <- "sig_genes_time_series";
  } else if (anal.type == "custom") {
    grp.nms <- strsplit(par1, " vs. ")[[1]]
    myargs[[1]] <- paste(grp.nms, collapse = "-")
    dataSet$grp.nms <- grp.nms;
    filename <- paste("sig_genes_", paste(grp.nms, collapse = "_vs_"), sep = "")
    dataSet$contrast <- paste(grp.nms, collapse = "_vs_");
  } else if (anal.type == "reference") {
    ref <- par1;
    cntr.cls <- grp.nms[grp.nms != ref]
    myargs <- as.list(paste(cntr.cls, "-", ref, sep = ""));
    dataSet$ref <- ref; 
    filename <- paste("sig_genes_reference_", ref, sep = "");
  } else if (anal.type == "nested") {
    grp.nms1 <- strsplit(par1, " vs. ")[[1]]
    grp.nms2 <- strsplit(par2, " vs. ")[[1]]
    if (all(grp.nms1 == grp.nms2)) {
      msgSet$current.msg <-"The two nested groups are the same. Please choose two different groups."
      saveSet(msgSet, "msgSet");      
      return(0)
    }
    grp.nms <- unique(c(grp.nms1, grp.nms2))
    if (nested.opt == "intonly") {
      dataSet$nested.int.opt <- "True";
      myargs[[1]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    } else {
      dataSet$nested.int.opt <- "False";
      myargs[[1]] <- paste(grp.nms1, collapse = "-")
      myargs[[2]] <- paste(grp.nms2, collapse = "-")
      myargs[[3]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    }
    dataSet$contrast <- paste(paste(paste(grp.nms1, collapse = "_vs_"), "_", paste(grp.nms2, collapse = "_vs_"), sep = ""), sep = "")
    filename <- paste("sig_genes_nested_", paste(paste(grp.nms1, collapse = "_vs_"), "_", paste(grp.nms2, collapse = "_vs_"), sep = ""), sep = "")
  } else {
    print(paste("Not supported: ", anal.type))
  }
    
  dataSet$filename <- paste0(filename, "_", dataSet$de.method);
  require(limma);
  design <- dataSet$design;
  myargs[["levels"]] <- design;
  dataSet$contrast.type <- anal.type;
  contrast.matrix <- do.call(makeContrasts, myargs);
  dataSet$contrast.matrix <- contrast.matrix;
  RegisterData(dataSet);
  return(dataSet);
}

.perform_limma_edger <- function(dataSet, robustTrend = FALSE) {
  ## ------------------------------------------------------------------ ##
  ## 1 · Input checks & dependencies                                    ##
  ## ------------------------------------------------------------------ ##
  require(limma)
  require(edgeR)

  if (is.null(dataSet$design) || is.null(dataSet$contrast.matrix)) {
    stop("design and/or contrast.matrix missing in dataSet. Run prepareEdgeRContrast() first.")
  }

  design           <- dataSet$design
  contrast.matrix  <- dataSet$contrast.matrix
  contrast.names   <- colnames(contrast.matrix)

  paramSet <- readSet(paramSet, "paramSet")
  msgSet   <- readSet(msgSet,   "msgSet")

  data.norm <- if (length(dataSet$rmidx) > 0) {
    dataSet$data.norm[, -dataSet$rmidx, drop = FALSE]
  } else {
    dataSet$data.norm
  }

  if (dataSet$de.method == "limma") {
    if (is.null(dataSet$block)) {
      fit <- lmFit(data.norm, design)
    } else {
      corfit <- duplicateCorrelation(data.norm, design, block = dataSet$block)
      fit <- lmFit(data.norm, design, block = dataSet$block, correlation = corfit$consensus)
    }
    
    if (!is.fullrank(design)) {
      msgSet$current.msg <- "This metadata combination is not full rank! Please use other combination.";
      saveSet(msgSet, "msgSet");  
      return(0)
    }
    
    df.residual <- fit$df.residual
    if (all(df.residual == 0)) {
      msgSet$current.msg <- "All residuals equal 0. There is not enough replicates in each group (no residual degrees of freedom)!";
      saveSet(msgSet, "msgSet");  
      return(0);
    }
    fit2 <- contrasts.fit(fit, contrast.matrix);
    fit2 <- eBayes(fit2, trend=robustTrend, robust=robustTrend);

result.list <- list()
for (nm in colnames(contrast.matrix)) {
  tbl <- topTable(fit2, coef = nm, number = Inf, adjust.method = "fdr")
  colnames(tbl)[colnames(tbl) == "FDR"] <- "adj.P.Val"
  result.list[[nm]] <- tbl
}

dataSet$comp.res.list      <- result.list     
  dataSet$comp.res <- result.list[[1]]

  } else if (dataSet$de.method == "edger") {

    set.seed(1)

    # Retrieve the raw (un‑normalised) count matrix
    cnt.mat <- .get.annotated.data()
    if (length(dataSet$rmidx) > 0)
      cnt.mat <- cnt.mat[, -dataSet$rmidx, drop = FALSE]

    grp.fac <- factor(dataSet$cls)

    if (!is.null(dataSet$block)) {
      blk.fac <- factor(dataSet$block)
      design  <- model.matrix(~ grp.fac + blk.fac)
    } else {
      # Use the stored design if created with ~0+grp.fac; otherwise rebuild
      if (is.null(attr(design, "assign")))
        design <- model.matrix(~ 0 + grp.fac)
    }

    y <- DGEList(counts = cnt.mat, group = grp.fac)
    y <- calcNormFactors(y)

    ## Dispersions
    y <- estimateGLMCommonDisp(y, design, verbose = FALSE)
    y <- tryCatch(
      estimateGLMTrendedDisp(y, design),
      error   = function(e) { msgSet$current.msg <- e$message ; saveSet(msgSet, "msgSet"); return(0) },
      warning = function(w) { msgSet$current.msg <- c(msgSet$current.msg, w$message); saveSet(msgSet, "msgSet"); }
    )
    y <- estimateGLMTagwiseDisp(y, design)

    fit <- glmFit(y, design)

    result.list <- vector("list", length(contrast.names))
    names(result.list) <- contrast.names
    for (nm in contrast.names) {
      lrt <- glmLRT(fit, contrast = contrast.matrix[, nm])
      tbl <- topTags(lrt, n = Inf)$table
      colnames(tbl)[colnames(tbl) == "FDR"]    <- "adj.P.Val"
      colnames(tbl)[colnames(tbl) == "PValue"] <- "P.Value"
      result.list[[nm]] <- tbl
    }
  dataSet$comp.res.list <- result.list
  dataSet$comp.res <- result.list[[1]]

}
  return(dataSet);

}
.perform_williams_trend <- function(dataSet,
                                    robustTrend = FALSE,
                                    verbose     = TRUE)
{

  expr <- if (length(dataSet$rmidx) > 0)
    dataSet$data.norm[, -dataSet$rmidx, drop = FALSE] else
      dataSet$data.norm

  cls_vals <- if (length(dataSet$rmidx) > 0)
    dataSet$cls[-dataSet$rmidx] else
      dataSet$cls

  unique_cls <- unique(cls_vals)
  dose_order <- suppressWarnings(as.numeric(as.character(unique_cls)))

  if (all(!is.na(dose_order))) {
    ord_levels <- unique_cls[order(dose_order)]
  } else {
    ord_levels <- sort(unique_cls)
  }

  grp <- factor(cls_vals, levels = ord_levels, ordered = TRUE)
  grp <- droplevels(grp)

  if (nlevels(grp) < 3L)
    stop("Williams trend test requires ≥ 3 ordered doses.")

  grp_levels <- levels(grp)
  grp_index  <- lapply(grp_levels, function(lv) which(grp == lv))
  n_per_grp  <- vapply(grp_index, length, integer(1))

  total_n  <- sum(n_per_grp)
  df_resid <- total_n - length(grp_levels)
  if (df_resid <= 0)
    stop("Williams trend test requires replication (df <= 0).")

  gene_ids   <- rownames(expr)
  gene_count <- nrow(expr)

  group_means <- matrix(NA_real_, nrow = gene_count,
                        ncol = length(grp_levels),
                        dimnames = list(gene_ids, grp_levels))
  SSE <- numeric(gene_count)

  for (i in seq_along(grp_levels)) {
    idx <- grp_index[[i]]
    mat <- expr[, idx, drop = FALSE]
    gm  <- rowMeans(mat)
    group_means[, i] <- gm

    if (length(idx) > 1L) {
      centered <- sweep(mat, 1, gm, "-", check.margin = FALSE)
      SSE <- SSE + rowSums(centered^2)
    }
  }

  MSE <- SSE / df_resid

  control_mean <- group_means[, 1, drop = TRUE]
  dose_means   <- group_means[, -1, drop = FALSE]
  pair_logFC   <- sweep(dose_means, 1, control_mean, "-")
  colnames(pair_logFC) <- paste0("Dose_", grp_levels[1],
                                 ".Dose_", grp_levels[-1])

  ave_expr <- rowMeans(expr)

  combo_idx <- lapply(seq_len(length(grp_levels) - 1L),
                      function(k) seq(from = length(grp_levels) - k + 1L,
                                      to   = length(grp_levels)))
  combo_n  <- vapply(combo_idx, function(idx) sum(n_per_grp[idx]), numeric(1))

  diff_mat <- matrix(NA_real_, nrow = gene_count,
                     ncol = length(combo_idx))
  for (j in seq_along(combo_idx)) {
    idx      <- combo_idx[[j]]
    weights  <- n_per_grp[idx] / sum(n_per_grp[idx])
    submeans <- group_means[, idx, drop = FALSE]
    weighted <- sweep(submeans, 2, weights, `*`)
    diff_mat[, j] <- rowSums(weighted) - control_mean
  }

  sqrt_terms <- sqrt(1 / n_per_grp[1] + 1 / combo_n)
  denom_mat  <- outer(sqrt(MSE), sqrt_terms, `*`)
  t_mat      <- diff_mat / denom_mat

  zero_mask <- denom_mat == 0
  if (any(zero_mask, na.rm = TRUE)) {
    t_mat[zero_mask & diff_mat > 0]  <- Inf
    t_mat[zero_mask & diff_mat < 0]  <- -Inf
    t_mat[zero_mask & diff_mat == 0] <- 0
  }

  t_pos <- apply(t_mat, 1, max, na.rm = TRUE)
  t_neg <- apply(t_mat, 1, min, na.rm = TRUE)
  t_pos[is.infinite(t_pos) & t_pos < 0] <- NA_real_
  t_neg[is.infinite(t_neg) & t_neg > 0] <- NA_real_

  both_na <- is.na(t_pos) & is.na(t_neg)
  abs_pos <- abs(t_pos); abs_pos[is.na(abs_pos)] <- -Inf
  abs_neg <- abs(t_neg); abs_neg[is.na(abs_neg)] <- -Inf
  use_pos <- abs_pos >= abs_neg
  t_stat  <- ifelse(use_pos, t_pos, t_neg)
  t_stat[both_na] <- NA_real_

  P.Value   <- 2 * pt(-abs(t_stat), df = df_resid)
  adj.P.Val <- p.adjust(P.Value, "fdr")

  topFeatures <- data.frame(
    pair_logFC,
    AveExpr   = ave_expr,
    t         = t_stat,
    P.Value   = P.Value,
    adj.P.Val = adj.P.Val,
    check.names = FALSE)

  ord <- order(topFeatures$adj.P.Val, na.last = TRUE)
  topFeatures <- topFeatures[ord, , drop = FALSE]

  dataSet$comp.res <- topFeatures
  dataSet$comp.res.list <- make_comp_res_list(
    topFeatures,
    stat.cols = c("AveExpr", "t", "P.Value", "adj.P.Val"))
  dataSet$de.method <- "wtt"

  if (verbose) {
    cat("Williams trend test (BMDExpress-style):",
        nrow(topFeatures), "features processed;",
        "df =", df_resid, "\n")
  }

  return(dataSet)
}



SetupDesignMatrix<-function(dataName="", deMethod){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  cls <- dataSet$cls; 
  design <- model.matrix(~ 0 + cls) # no intercept
  colnames(design) <- levels(cls);
  dataSet$design <- design;
  dataSet$de.method <- deMethod;
  dataSet$pval <- 0.05;
  dataSet$fc.val <- 1;

  saveSet(paramSet, "paramSet");
  return(RegisterData(dataSet));
}


# perform limma on given two groups selected 
# used by integarative analysis

#'Perform differential analysis using Limma method (meta-analysis)
#'@description Detect number of DE genes in individual data set for quality checking purposes 
#'@param dataName File name of data set.
#'@param grps Selected two groups for comparison
#'@param p.lvl P-value threshold
#'@param fc.lvl Fold-change threshold
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
PerformLimmaDE<-function(dataName="", grps, p.lvl, fc.lvl=NULL){
  #save.image("limma.RDAta");
  dataSet <- readDataset(dataName);
  dataSet$pval <- p.lvl;
    dataSet$fc.lvl <- 0;

  if(length(levels(dataSet$cls))>2){ 
    grp.nms <- strsplit(grps, " vs. ")[[1]];
    sel.inx <- as.character(dataSet$cls) %in% grp.nms;
  }else{
    sel.inx <- rep(T, ncol(dataSet$data.norm));
  }
  
  group <- factor(dataSet$cls[sel.inx]); # note regenerate factor to drop levels 
  data <- dataSet$data.norm[, sel.inx];
  
  res.limma <- PerformLimma(data, group);
  res.all <- GetLimmaResTable(res.limma$fit.obj);
  
  if(!is.null(fc.lvl)){
    hit.inx <- abs(res.all$logFC)>= fc.lvl & res.all$adj.P.Val <= p.lvl
    dataSet$fc.lvl <- fc.lvl;

  }else{
    hit.inx <- res.all$adj.P.Val <= p.lvl
    dataSet$fc.lvl <- 0;

  }
  if(sum(hit.inx) == 0){
    return (c(1, 0, nrow(res.all)));
  }
  # note, hit.inx can contain NA, not T/F
  hit.inx <- which(hit.inx);
  res <- res.all[hit.inx, , drop=F];
  
  # rm .txt suffix for new names
  shortNm <- substring(dataName, 0, nchar(dataName)-4);
  fast.write(signif(res[,-1],5), file=paste("SigGenes_", shortNm, ".csv",sep=""));
  
  sig.count <- nrow(res);
  de.genes <- rownames(res);
  non.sig.count <- nrow(data)-sig.count;
  rm(res.all);
  
  gc();
  
  # record the sig gene vec
  output <- c(1, sig.count, non.sig.count);
  return(RegisterData(dataSet, output));
}

# perfor differential analysis for array/RNA seq data
# for two groups only (used for meta-analysis)
PerformLimma<-function(data, group){
  #print(identical(colnames(data), group))
  #print(colnames(data))
  #print(group)

  #print("limma");
  require(limma);
  design <- model.matrix(~-1 + group);
  fit <- lmFit(data, design);
  
  grps.cmp <- paste("group", levels(group)[2], " - ", "group", levels(group)[1], sep="");
  myargs <- list(grps.cmp, levels = design);
  contrast.matrix <- do.call(makeContrasts, myargs);
  fit <- contrasts.fit(fit, contrast.matrix);
  fit <- eBayes(fit);
  gc();
  return (list(fit.obj=fit));
}

# get result table from eBayes fit object
GetLimmaResTable<-function(fit.obj){
  require(limma);
  resTable <- topTable(fit.obj, number=Inf, adjust.method="BH");
  if(!is.null(resTable$ID)){ # for older version
    rownames(resTable) <- resTable$ID;
    resTable$ID <- NULL;
  }
  return (resTable);
}

MultiCovariateRegression <- function(fileName,
                                     analysis.var, # metadata variable name
                                     ref = NULL, # reference class from analysis.var metadata (only if categorical)
                                     contrast = "anova",  # comparison class from analysis.var (only if categorical)
                                     blocking.factor = NULL, 
                                     robustTrend = F, 
                                     internal=F,
                                     useBatchCorrected=TRUE){ # whether returns 0/1 or dataset object

  dataSet <- readDataset(fileName);
  if(!exists('adj.vec')){ ## covariate to adjust for
    adj.vec <- "NA";
  }else{
    if(length(adj.vec) > 0){

    }else{
      adj.vec <- "NA"
    }
  }
  interim <- .multiCovariateRegression(dataSet, analysis.var, ref, contrast, blocking.factor, adj.vec, robustTrend, F, useBatchCorrected)
  if(is.list(interim)){
    res <- 1;
  }else{
    res <- interim;
  }
  return(res);
}

.multiCovariateRegression <- function(dataSet,
                                     analysis.var, # metadata variable name
                                     ref = NULL, # reference class from analysis.var metadata (only if categorical)
                                     contrast = "anova",  # comparison class from analysis.var (only if categorical)
                                     # fixed.effects = NULL,  
                                     blocking.factor = NULL, 
                                     adj.factors="NA",# metadata variables to adjust for
                                     robustTrend = F, 
                                     internal=F, # whether returns 0/1 or dataset object, T for metaanal covariate
                                     useBatchCorrected=TRUE){
  # load libraries
  library(limma);
  library(dplyr);
  
  # need a line for read dataSet
  msgSet <- readSet(msgSet, "msgSet");
  dataSet$rmidx <- NULL;

  # for embedded inside tools (ExpressAnalyst etc)
  paramSet <- readSet(paramSet, "paramSet");
  useMeta <- !is.null(paramSet$performedBatch) && isTRUE(paramSet$performedBatch);
  if(internal){
    if(useMeta){
    print("batchadjustedcov");
    inmex.meta <- qs::qread("inmex_meta.qs");
    }else{
    print("notbatchadjusted");

    inmex.meta <- qs::qread("inmex.meta.orig.qs");

    }
    feature_table <- inmex.meta$data[, colnames(inmex.meta$data) %in% colnames(dataSet$data.norm)];
  }else{
    feature_table <- dataSet$data.norm;
  }
  covariates <- dataSet$meta.info;

  matched_indices <- match(colnames(feature_table), rownames(covariates));
  covariates <- covariates[matched_indices, ,drop=F ];
  dataSet$meta.info <- covariates;
  #fixed.effects <- adj.vec
  # process covariates
  var.types <- lapply(covariates, class) %>% unlist();
  covariates[,c(var.types == "character")] <- lapply(covariates[,c(var.types == "character")], factor);

  # aggregate vars
  all.vars <- c(analysis.var);
  vars <- c(analysis.var);

  if(all(adj.factors != "NA")){
    vars = c(vars, adj.factors);
  }
  
  if(!is.null(blocking.factor) && !is.na(blocking.factor) && blocking.factor!="NA" && blocking.factor!="" ){
    all.vars = c(all.vars, blocking.factor);
  }

  all.vars<- unique(all.vars);
    
  covariates <- covariates[,unique(c(vars, all.vars)),drop=F];
  rmidx <-which(apply(covariates, 1, function(x) "NA" %in% x))
  
  if(length(rmidx)>0){
    covariates <- covariates[-rmidx,,drop=F];
    dataSet$rmidx <- rmidx;
  }
  feature_table <- feature_table[,colnames(feature_table) %in% rownames(covariates)];
  
  if(!identical(colnames(feature_table), rownames(covariates))){
    msgSet$current.msg <- "Error - order of samples got mixed up between tables";
    saveSet(msgSet, "msgSet");
    return(0)
  }
  
  
  # get analysis type
  analysis.type = ifelse(dataSet$disc.inx[analysis.var],"disc","cont")
  if(is.na(analysis.type)){
    msgSet$current.msg <- "Analysis var not found in our database!";
    saveSet(msgSet, "msgSet");
    return(0)
  }
  
  if(analysis.type == "disc"){
    # build design and contrast matrix
    #  covariates[, analysis.var] <- covariates[, analysis.var] %>% make.names() %>% factor();
    #str(covariates)
    grp.nms <- levels(covariates[, analysis.var]);
    
    if(any(grepl("(^[0-9]+).*", grp.nms))){
      grp.nms <- paste0(analysis.var,"_",grp.nms);
      if(!(is.null(ref))& ref!="NA"){
        ref <- paste0(analysis.var,"_", ref)
      }
      if(contrast!="anova"){
        contrast <- paste0(analysis.var,"_", contrast)
      }     
    }
    
    for(col in 1:ncol(covariates)){
      if(dataSet$cont.inx[colnames(covariates)[col]]){
        covariates[,col] <- as.numeric(as.character(covariates[,col]))
      }
    }
    
    design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    colnames(design)[1:length(grp.nms)] <- grp.nms;
    myargs <- list();
    if(contrast == "anova"){ 
      contrasts <- grp.nms[grp.nms != ref];
      myargs <- as.list(paste(contrasts, "-", ref, sep = "")); 
    } else {
      myargs <- as.list(paste(contrast, "-", ref, sep = ""));
    }
    myargs[["levels"]] <- design;
    contrast.matrix <- do.call(makeContrasts, myargs);
    
    # handle blocking factor
    if (is.null(blocking.factor) | is.na(blocking.factor) | blocking.factor=="NA" | blocking.factor == "") {
      fit <- lmFit(feature_table, design);
    } else {
      block.vec <- covariates[,blocking.factor];
      corfit <- duplicateCorrelation(feature_table, design, block = block.vec);
      fit <- lmFit(feature_table, design, block = block.vec, correlation = corfit$consensus);
    }

    # get results
    fit <- contrasts.fit(fit, contrast.matrix);
    fit <- eBayes(fit, trend=robustTrend, robust=robustTrend);
    rest <- topTable(fit, number = Inf);
    
    if(contrast != "anova"){    
      colnames(rest)[1] <- myargs[[1]];
      grp.nms<-c(ref,contrast)
      
    }
    #for meta-anal
    if(internal){
        ### get results with no adjustment
        design <- model.matrix(formula(paste0("~ 0", paste0(" + ", analysis.var, collapse = ""))), data = covariates);
        colnames(design)[1:length(grp.nms)] <- grp.nms;
        myargs[["levels"]] <- design;
        contrast.matrix <- do.call(makeContrasts, myargs);
        fit <- lmFit(feature_table, design)
        fit <- contrasts.fit(fit, contrast.matrix);
        fit <- eBayes(fit, trend=robustTrend, robust=robustTrend);
        res.noadj <- topTable(fit, number = Inf);
        dataSet$res.noadj <- res.noadj;
    }

    dataSet$contrast.matrix <- contrast.matrix;
    dataSet$par1 <-  myargs[[1]];
    dataSet$grp.nms <- ifelse(any(grepl("(^[0-9]+).*", grp.nms)), paste0(analysis.var,"_",grp.nms),grp.nms);
  } else { 
    
    # build design matrix
    types <- dataSet$cont.inx[vars];
    contIdx <- as.numeric(which(types))
    covariates[,contIdx] <- unlist(lapply(covariates[,contIdx], function(x) as.numeric(as.character((x)))));
    
    if (all(types)) {
      design <- model.matrix(formula(paste0("~", paste0(" + ", vars, collapse = ""))), data = covariates);
    } else {
      design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    }
    
    # handle blocking factor
    if (is.null(blocking.factor) | is.na(blocking.factor) | blocking.factor=="NA") {
      fit <- lmFit(feature_table, design);
    } else {
      block.vec <- covariates[, blocking.factor];
      corfit <- duplicateCorrelation(feature_table, design, block = block.vec);
      fit <- lmFit(feature_table, design, block = block.vec, correlation = corfit$consensus);
    }
    
    # get results
    fit <- eBayes(fit, trend=robustTrend, robust=robustTrend);
    rest <- topTable(fit, number = Inf, coef = analysis.var);
    colnames(rest)[1] <- dataSet$par1 <- analysis.var;
    
    ### get results with no adjustment
    if(internal){
    design <- model.matrix(formula(paste0("~", analysis.var)), data = covariates);
    fit <- eBayes(lmFit(feature_table, design), trend=robustTrend, robust=robustTrend);
    res.noadj <- topTable(fit, number = Inf);
    dataSet$res.noadj <- res.noadj;
    }
  }
  
  dataSet$design <- design;
  dataSet$contrast.type <- analysis.type;
  dataSet$comp.res <- rest;
  dataSet$comp.res.list <- make_comp_res_list(rest)
  dataSet$de.method <- "limma"
  dataSet$comp.type <- "default"
  dataSet$fit.obj <- fit;

  dataSet$pval <- 0.05;
  dataSet$fc.val <- 1;
  dataSet$analysis.var <- analysis.var;
  dataSet$de.adj <- adj.factors;

  if(all(adj.factors != "NA")){
    dataSet$de.adj <- "NA"
  }else{
    dataSet$de.adj <- adj.factors;
  }

  RegisterData(dataSet);
  return(dataSet);
}
make_comp_res_list <- function(resTab,
                               stat.cols = c("AveExpr", "F", "t",
                                             "P.Value", "adj.P.Val",
                                             "B", "FDR"))
{
  ## detect logFC columns automatically -------------------------------
  lfc.cand <- setdiff(colnames(resTab), stat.cols)

  ## strip common suffix/prefix for the test
  strip_logfc <- function(x)
      sub("\\.logFC$","",
      sub("^logFC\\.","", x, ignore.case = TRUE), ignore.case = TRUE)

  uniq.core <- unique(strip_logfc(lfc.cand))

  if (length(uniq.core) == 0L)
      stop("Could not detect any log-fold-change columns automatically. ",
           "Pass lfc.cols explicitly.")

  ## build list --------------------------------------------------------
  out <- lapply(uniq.core, function(core) {

            ## match any of the naming patterns for this contrast
            pat <- paste0("^(", core,
                          "|logFC\\.", core,
                          "|", core, "\\.logFC)$")
            this.lfc <- grep(pat, colnames(resTab), value = TRUE)

            if (length(this.lfc) == 0L)
                stop("Unexpected: no column matched for contrast ", core)

            df <- resTab[ , c(this.lfc[1], stat.cols[stat.cols %in% colnames(resTab)]),
                          drop = FALSE]
            colnames(df)[1] <- "logFC"
            df
         })
  names(out) <- uniq.core
  out
}

parse_contrast_groups <- function(contrast_str) {
  comps <- strsplit(contrast_str, " vs\\.?\\s*")[[1]]
  if (length(comps) != 2) stop(paste("Invalid contrast format:", contrast_str))
  return(comps)
}

.get.interaction.results <- function(dds.path = "deseq.res.obj.rds") {
  dds <- qs::qread(dds.path)
  cat("Available result names:\n")
  
  # Automatically detect the interaction term
  interaction_name <- grep("factorA.*factorB.*", resultsNames(dds), value = TRUE)
  if (length(interaction_name) == 0) {
    stop("No interaction term found in model.")
  }

  cat("Extracting interaction term:", interaction_name, "\n")
  res <- results(dds, name = interaction_name[1], independentFiltering = FALSE, cooksCutoff = Inf)

  # Format results
  topFeatures <- data.frame(res@listData)
  rownames(topFeatures) <- rownames(res)
  colnames(topFeatures) <- sub("padj", "adj.P.Val", colnames(topFeatures))
  colnames(topFeatures) <- sub("pvalue", "P.Value", colnames(topFeatures))
  colnames(topFeatures) <- sub("log2FoldChange", "logFC", colnames(topFeatures))
  topFeatures <- topFeatures[c("logFC", "baseMean", "lfcSE", "stat", "P.Value", "adj.P.Val")]
  topFeatures <- topFeatures[order(topFeatures$P.Value), ]

  return(topFeatures)
}

prepareEdgeRContrast <- function(dataSet,
                                 anal.type  = "reference",
                                 par1       = NULL,
                                 par2       = NULL,
                                 nested.opt = "intonly") {
  msgSet <- readSet(msgSet, "msgSet")
  set.seed(1337)
  require(limma)

  cls_raw      <- factor(dataSet$cls)
  raw_levels   <- levels(cls_raw)
  syn_levels   <- make.names(raw_levels)

  cls_syn      <- factor(cls_raw, levels = raw_levels, labels = syn_levels)

  dataSet$cls_raw   <- cls_raw
  dataSet$cls       <- cls_syn               # <- modeling uses sanitized
  dataSet$grp.nms   <- syn_levels
  dataSet$grp.nms_raw <- raw_levels
  dataSet$comp.type <- anal.type
  dataSet$par1      <- par1

  design <- model.matrix(~ 0 + cls_syn)
  colnames(design) <- syn_levels

  # helper: map a UI label (raw or syn) into sanitized
  to_syn <- function(x) {
    if (is.null(x)) return(NULL)
    x <- trimws(as.character(x))
    if (!nzchar(x)) return(NULL)
    if (x %in% raw_levels) return(make.names(x))
    x
  }

  # parse "A vs. B" (accept raw or syn on either side)
  parse_vs <- function(x) {
    if (is.null(x)) return(c(NA_character_, NA_character_))
    parts <- strsplit(x, "\\s*vs\\.?\\s*", perl = TRUE)[[1]]
    if (length(parts) != 2) return(c(NA_character_, NA_character_))
    c(to_syn(parts[1]), to_syn(parts[2]))
  }

  ## build contrasts using sanitized tokens
  if (anal.type == "reference") {
    ref_syn <- to_syn(par1)
    if (is.null(ref_syn) || !(ref_syn %in% syn_levels)) {
      stop("`par1` must specify a valid reference level. You gave '",
           if (is.null(par1)) "NULL" else par1,
           "'. Valid (raw): ", paste(raw_levels, collapse = ", "))
    }
    others <- setdiff(syn_levels, ref_syn)
    conts  <- setNames(lapply(others, \(g) paste0(g, " - ", ref_syn)),
                       paste0(others, "_vs_", ref_syn))

  } else if (anal.type == "default") {
    combs <- combn(syn_levels, 2, simplify = FALSE)
    conts <- setNames(lapply(combs, \(x) paste0(x[1], " - ", x[2])),
                      sapply(combs, \(x) paste0(x[1], "_vs_", x[2])))

  } else if (anal.type == "time") {
    tm <- syn_levels
    conts <- setNames(lapply(seq_len(length(tm) - 1L),
                             \(i) paste0(tm[i + 1L], " - ", tm[i])),
                      paste0(tm[-1], "_vs_", tm[-length(tm)]))

  } else if (anal.type == "custom") {
    grp <- parse_vs(par1)
    if (anyNA(grp) || !all(grp %in% syn_levels)) {
      stop("`par1` must be 'A vs. B'. Valid (raw): ",
           paste(raw_levels, collapse = ", "))
    }
    conts <- setNames(list(paste0(grp[2], " - ", grp[1])),
                      paste0(grp[2], "_vs_", grp[1]))

  } else if (anal.type == "nested") {
    g1 <- parse_vs(par1); g2 <- parse_vs(par2)
    if (anyNA(g1) || anyNA(g2) || !all(c(g1, g2) %in% syn_levels)) {
      stop("`par1` and `par2` must be 'A vs. B'. Valid (raw): ",
           paste(raw_levels, collapse = ", "))
    }
    if (identical(nested.opt, "intonly")) {
      expr  <- paste0("(", g1[1], " - ", g1[2], ") - (", g2[1], " - ", g2[2], ")")
      nm    <- paste0(g1[1], g1[2], "_vs_", g2[1], g2[2], "_interaction")
      conts <- setNames(list(expr), nm)
    } else {
      expr1 <- paste0(g1[2], " - ", g1[1])
      expr2 <- paste0(g2[2], " - ", g2[1])
      expr3 <- paste0("(", g1[2], " - ", g1[1], ") - (", g2[2], " - ", g2[1], ")")
      conts <- c(setNames(list(expr1), paste0(g1[2], "_vs_", g1[1])),
                 setNames(list(expr2), paste0(g2[2], "_vs_", g2[1])),
                 setNames(list(expr3), paste0("int_", g1[2], g1[1], "_vs_", g2[2], g2[1])))
    }

  } else {
    stop("Unsupported `anal.type`: ", anal.type)
  }

  contrast.matrix <- do.call(makeContrasts, c(conts, list(levels = design)))

  dataSet$design          <- design
  dataSet$contrast.matrix <- contrast.matrix
  dataSet$contrast.names  <- colnames(contrast.matrix)
  dataSet$contrast.type   <- anal.type
  dataSet$filename        <- paste0("edgeR_", anal.type, "_", dataSet$de.method)

  RegisterData(dataSet)
  dataSet
}
