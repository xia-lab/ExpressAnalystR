##################################################
## R script for ExpressAnalyst
## Description: Compute upset diagram
## Authors: 
## G. Zhou, guangyan.zhou@mail.mcgill.ca
###################################################


#'Prepare data for Upset diagram
#'@param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#'@param fileNm file name of the json file output 
#'@export
PrepareUpsetData <- function(fileNm){
  print(paste("[DEBUG PrepareUpsetData] Starting with fileNm:", fileNm));
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  mdata.all <- paramSet$mdata.all;
  anal.type <- paramSet$anal.type;
  print(paste("[DEBUG PrepareUpsetData] anal.type:", anal.type));

  newDat <- list();

  # selected dataset or comparisons for onedata (single gene expression matrix)
  if(anal.type == "metadata"){
  print("[DEBUG PrepareUpsetData] Processing metadata type");
  hit.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[hit.inx];
  }else if(anal.type == "onedata"){
    # For onedata, use comparison-based upset
    print("[DEBUG PrepareUpsetData] Delegating to PrepareUpsetDataFromComparisons");
    return(PrepareUpsetDataFromComparisons(fileNm));
  }else{
  print("[DEBUG PrepareUpsetData] Processing genelist type");
  sel.nms <- names(mdata.all)
  }
  sel.dats <- list();


    if(!exists("analSet$inmex.de")){
      analSet$inmex.de <- list();
    }

  # populate gene lists for upset plot based on selected names
  for(nm in sel.nms){
    if(anal.type == "metadata"){
      dataSet <- readDataset(nm);
      sel.dats[[nm]] <- dataSet$sig.mat$ids;
    }else if(anal.type == "onedata"){
      # This branch should not be reached due to early return above
    }else{
      dataSet <- readDataset(nm);
      gene.mat <- dataSet$prot.mat;

      # convert to entrez
      expr.val <- gene.mat[,1];
      en.ids <- rownames(gene.mat);

      names(expr.val) <- en.ids;
      analSet$inmex.de[[nm]] <- en.ids;
      sel.dats[[nm]] <- en.ids;
    }

  }

  if(anal.type == "metadata" & paramSet$meta.selected){
    sel.dats[["meta_dat"]] <- as.character(rownames(analSet$meta.mat));
  }
  
  if(length(sel.dats) == 0){
    AddErrMsg("No signficant features for any dataset!");
    return(0);
  }
  
  sums <- unlist(lapply(sel.dats, length));
  names <- unlist(lapply(sel.dats, paste, collapse = ", "));
  if(anal.type == "metadata"){
    metasum <- length(analSet$meta.stat$idd);
    metaname <- paste(analSet$meta.stat$idd, collapse = ", ");
    allsums <- c(sums, metasum);
    allnames <- c(names, metaname);
  }else{
    allsums <- c(sums);
    allnames <- c(names);
  }


  require(reshape2)
  df <- reshape::melt(sel.dats, value.name="id")
  colnames(df) <- c("name", 'set')
  uniq.nms <- unique(df$name)
  new.df <- dcast(df, name ~ set, value.var='set', fill=0)
  rownames(new.df) <- new.df[,1]
  new.df <- new.df[,-1, drop=F]
  
  gene.map <-  queryGeneDB("entrez", paramSet$data.org);
  gene.map[] <- lapply(gene.map, as.character)
  
  json.list <- list()
  for(i in 1:nrow(new.df)){
    json.list[[i]] <- list()
    json.list[[i]][["sets"]] <- new.df[i,][new.df[i,] != 0]
    entrez.vec <- rownames(new.df)[i];

    # Handle case where gene.map might not be a proper matrix/dataframe
    symbols <- entrez.vec;  # Default to using original ID
    tryCatch({
      if(!is.null(gene.map) && is.matrix(gene.map) || is.data.frame(gene.map)) {
        if(ncol(gene.map) >= 2 && "gene_id" %in% colnames(gene.map) && "symbol" %in% colnames(gene.map)) {
          hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
          symbols <- gene.map[hit.inx, "symbol"];
          # if not gene symbol, use id by itself
          na.inx <- is.na(symbols);
          symbols[na.inx] <- entrez.vec[na.inx];
        }
      }
    }, error = function(e) {
      # If conversion fails, just use original ID
      symbols <<- entrez.vec;
    })

    json.list[[i]][["name"]] <- symbols;
    json.list[[i]][["entrez"]] <- entrez.vec;
  }
  
  col.vec <-gg_color_hue(length(sel.dats));
  
  jsonNm <- paste0(fileNm, ".json")
  json.mat <- RJSONIO::toJSON(list(json.list, col.vec));
  sink(jsonNm);
  cat(json.mat);
  sink();
  
  return(1); 
}

#'Prepare data for Upset diagram using DE comparisons (onedata)
#'@param fileNm file name of the json file output
#'@param dataName optional dataset name (defaults to paramSet$dataName)
#'@export
PrepareUpsetDataFromComparisons <- function(fileNm, dataName = ""){
  print(paste("[DEBUG PrepareUpsetDataFromComparisons] Starting with fileNm:", fileNm, "dataName:", dataName));
  paramSet <- readSet(paramSet, "paramSet");

  # If dataName is empty, use the current dataName from paramSet
  if(dataName == "" || is.null(dataName)){
    dataName <- paramSet$dataName;
    print(paste("[DEBUG PrepareUpsetDataFromComparisons] Using paramSet$dataName:", dataName));
  }

  dataSet <- readDataset(dataName);
  print(paste("[DEBUG PrepareUpsetDataFromComparisons] dataSet loaded, class:", class(dataSet)));

  if (is.null(dataSet) || is.null(dataSet$comp.res.list) ||
      length(dataSet$comp.res.list) == 0) {
    print("[DEBUG PrepareUpsetDataFromComparisons] ERROR: No comp.res.list found");
    AddErrMsg("No DE comparison results found.");
    return(0);
  }

  comp_list <- dataSet$comp.res.list
  comp_names <- names(comp_list)
  print(paste("[DEBUG PrepareUpsetDataFromComparisons] Found", length(comp_list), "comparisons:", paste(comp_names, collapse=", ")));
  if (is.null(comp_names) || length(comp_names) == 0) {
    comp_names <- paste0("Comparison_", seq_along(comp_list))
  }
  print(paste("[DEBUG PrepareUpsetDataFromComparisons] upset.comp.nms from paramSet:", paste(paramSet$upset.comp.nms, collapse=", ")));
  if (!is.null(paramSet$upset.comp.nms) && length(paramSet$upset.comp.nms) > 0) {
    keep <- comp_names %in% paramSet$upset.comp.nms
    print(paste("[DEBUG PrepareUpsetDataFromComparisons] Filtering to selected comparisons, keep:", paste(which(keep), collapse=", ")));
    comp_list <- comp_list[keep]
    comp_names <- names(comp_list)
    print(paste("[DEBUG PrepareUpsetDataFromComparisons] After filtering:", length(comp_list), "comparisons:", paste(comp_names, collapse=", ")));
  }
  if (is.null(comp_names) || length(comp_names) == 0) {
    print("[DEBUG PrepareUpsetDataFromComparisons] ERROR: No comparisons after filtering");
    AddErrMsg("No comparisons selected.");
    return(0);
  }

  p.lvl <- paramSet$BHth
  if (is.null(p.lvl) || !is.finite(p.lvl)) {
    p.lvl <- 0.05
  }
  fc.lvl <- dataSet$fc.lvl
  if (is.null(fc.lvl) || !is.finite(fc.lvl)) {
    fc.lvl <- paramSet$fc.thresh
  }
  if (is.null(fc.lvl) || !is.finite(fc.lvl)) {
    fc.lvl <- 0
  }
  print(paste("[DEBUG PrepareUpsetDataFromComparisons] Thresholds: p.lvl =", p.lvl, "fc.lvl =", fc.lvl));

  use_fdr <- if (is.null(paramSet$use.fdr)) TRUE else isTRUE(paramSet$use.fdr)
  print(paste("[DEBUG PrepareUpsetDataFromComparisons] use_fdr:", use_fdr));

  sel.dats <- list()
  for (i in seq_along(comp_list)) {
    rt <- comp_list[[i]]
    if (is.null(rt) || nrow(rt) == 0) {
      print(paste("[DEBUG PrepareUpsetDataFromComparisons] Comparison", comp_names[[i]], "has no data"));
      sel.dats[[comp_names[[i]]]] <- character(0)
      next
    }
    print(paste("[DEBUG PrepareUpsetDataFromComparisons] Processing comparison", comp_names[[i]], "with", nrow(rt), "features"));

    pcol <- NULL
    if (use_fdr && "adj.P.Val" %in% names(rt)) {
      pcol <- "adj.P.Val"
    } else if ("P.Value" %in% names(rt)) {
      pcol <- "P.Value"
    } else if ("padj" %in% names(rt)) {
      pcol <- "padj"
    } else if ("pvalue" %in% names(rt)) {
      pcol <- "pvalue"
    }

    lfc_col <- if ("logFC" %in% names(rt)) {
      "logFC"
    } else if ("log2FoldChange" %in% names(rt)) {
      "log2FoldChange"
    } else {
      NULL
    }

    deg_pass <- if (!is.null(pcol)) {
      is.finite(rt[[pcol]]) & (rt[[pcol]] <= p.lvl)
    } else {
      rep(TRUE, nrow(rt))
    }
    lfc_pass <- if (!is.null(lfc_col)) {
      is.finite(rt[[lfc_col]]) & (abs(rt[[lfc_col]]) >= fc.lvl)
    } else {
      rep(TRUE, nrow(rt))
    }

    sig_genes <- as.character(rownames(rt)[deg_pass & lfc_pass])
    print(paste("[DEBUG PrepareUpsetDataFromComparisons] Comparison", comp_names[[i]], "has", length(sig_genes), "significant features"));
    sel.dats[[comp_names[[i]]]] <- sig_genes
  }

  sel.dats <- sel.dats[sapply(sel.dats, length) > 0]
  print(paste("[DEBUG PrepareUpsetDataFromComparisons] After filtering empty sets:", length(sel.dats), "comparisons with data"));
  if (length(sel.dats) == 0) {
    print("[DEBUG PrepareUpsetDataFromComparisons] ERROR: No significant features");
    AddErrMsg("No significant features for any comparison.");
    return(0);
  }
  print(paste("[DEBUG PrepareUpsetDataFromComparisons] Final comparison names:", paste(names(sel.dats), collapse=", ")));

  require(reshape2)
  df <- reshape::melt(sel.dats, value.name="id")
  colnames(df) <- c("name", 'set')
  new.df <- dcast(df, name ~ set, value.var='set', fill=0)
  rownames(new.df) <- new.df[,1]
  new.df <- new.df[,-1, drop=FALSE]

  gene.map <- queryGeneDB("entrez", paramSet$data.org);
  gene.map[] <- lapply(gene.map, as.character)

  json.list <- list()
  for(i in 1:nrow(new.df)){
    json.list[[i]] <- list()
    json.list[[i]][["sets"]] <- new.df[i,][new.df[i,] != 0]
    entrez.vec <- rownames(new.df)[i]

    # Handle case where gene.map might not be a proper matrix/dataframe
    symbols <- entrez.vec;  # Default to using original ID
    tryCatch({
      if(!is.null(gene.map) && (is.matrix(gene.map) || is.data.frame(gene.map))) {
        if(ncol(gene.map) >= 2 && "gene_id" %in% colnames(gene.map) && "symbol" %in% colnames(gene.map)) {
          hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
          symbols <- gene.map[hit.inx, "symbol"];
          # if not gene symbol, use id by itself
          na.inx <- is.na(symbols);
          symbols[na.inx] <- entrez.vec[na.inx];
        }
      }
    }, error = function(e) {
      # If conversion fails, just use original ID
      symbols <<- entrez.vec;
    })

    json.list[[i]][["name"]] <- symbols
    json.list[[i]][["entrez"]] <- entrez.vec
  }

  col.vec <- gg_color_hue(length(sel.dats))

  jsonNm <- paste0(fileNm, ".json")
  json.mat <- RJSONIO::toJSON(list(json.list, col.vec))
  sink(jsonNm)
  cat(json.mat)
  sink()

  return(1)
}

#'Get significant counts for each comparison (for display purposes)
#'@param dataName optional dataset name (defaults to paramSet$dataName)
#'@export
GetComparisonSigCounts <- function(dataName = "") {
  paramSet <- readSet(paramSet, "paramSet");
  dataSet <- readDataset(dataName);

  if (is.null(dataSet) || is.null(dataSet$comp.res.list) ||
      length(dataSet$comp.res.list) == 0) {
    return(integer(0))
  }

  comp_list <- dataSet$comp.res.list
  comp_names <- names(comp_list)
  if (is.null(comp_names) || length(comp_names) == 0) {
    comp_names <- paste0("Comparison_", seq_along(comp_list))
  }

  p.lvl <- paramSet$BHth
  if (is.null(p.lvl) || !is.finite(p.lvl)) {
    p.lvl <- 0.05
  }
  fc.lvl <- dataSet$fc.lvl
  if (is.null(fc.lvl) || !is.finite(fc.lvl)) {
    fc.lvl <- paramSet$fc.thresh
  }
  if (is.null(fc.lvl) || !is.finite(fc.lvl)) {
    fc.lvl <- 0
  }

  use_fdr <- if (is.null(paramSet$use.fdr)) TRUE else isTRUE(paramSet$use.fdr)

  counts <- integer(length(comp_list))
  for (i in seq_along(comp_list)) {
    rt <- comp_list[[i]]
    if (is.null(rt) || nrow(rt) == 0) {
      counts[[i]] <- 0
      next
    }

    pcol <- NULL
    if (use_fdr && "adj.P.Val" %in% names(rt)) {
      pcol <- "adj.P.Val"
    } else if ("P.Value" %in% names(rt)) {
      pcol <- "P.Value"
    } else if ("padj" %in% names(rt)) {
      pcol <- "padj"
    } else if ("pvalue" %in% names(rt)) {
      pcol <- "pvalue"
    }

    lfc_col <- if ("logFC" %in% names(rt)) {
      "logFC"
    } else if ("log2FoldChange" %in% names(rt)) {
      "log2FoldChange"
    } else {
      NULL
    }

    deg_pass <- if (!is.null(pcol)) {
      is.finite(rt[[pcol]]) & (rt[[pcol]] <= p.lvl)
    } else {
      rep(TRUE, nrow(rt))
    }
    lfc_pass <- if (!is.null(lfc_col)) {
      is.finite(rt[[lfc_col]]) & (abs(rt[[lfc_col]]) >= fc.lvl)
    } else {
      rep(TRUE, nrow(rt))
    }

    counts[[i]] <- sum(deg_pass & lfc_pass, na.rm = TRUE)
  }

  names(counts) <- comp_names
  counts
}

#'Set selected comparisons for upset plot (onedata mode)
#'@export
SetUpsetComparisons <- function(){
  print("[DEBUG SetUpsetComparisons] Starting");
  paramSet <- readSet(paramSet, "paramSet");
  if (!exists("comp.vec")) {
    print("[DEBUG SetUpsetComparisons] ERROR: comp.vec does not exist");
    paramSet$upset.comp.nms <- NULL
    saveSet(paramSet, "paramSet")
    return(0)
  }
  print(paste("[DEBUG SetUpsetComparisons] comp.vec:", paste(comp.vec, collapse=", ")));
  paramSet$upset.comp.nms <- comp.vec
  print(paste("[DEBUG SetUpsetComparisons] Set paramSet$upset.comp.nms to:", paste(paramSet$upset.comp.nms, collapse=", ")));
  rm("comp.vec", envir = .GlobalEnv)
  saveSet(paramSet, "paramSet")
  print("[DEBUG SetUpsetComparisons] Saved and returning 1");
  return(1)
}

#Record upset intersection mode for report
SetUpsetMode <- function(mode){
      paramSet <- readSet(paramSet, "paramSet");
      paramSet$upsetMode <- mode;
  saveSet(paramSet, "paramSet");
}