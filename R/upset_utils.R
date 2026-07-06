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
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  mdata.all <- paramSet$mdata.all;
  anal.type <- paramSet$anal.type;

  newDat <- list();

  # selected dataset or comparisons for onedata (single gene expression matrix)
  if(anal.type == "metadata"){
  hit.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[hit.inx];
  }else if(anal.type == "onedata"){
    # For onedata, use comparison-based upset
    return(PrepareUpsetDataFromComparisons(fileNm));
  }else{
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
  paramSet <- readSet(paramSet, "paramSet");

  # If dataName is empty, use the current dataName from paramSet
  if(dataName == "" || is.null(dataName)){
    dataName <- paramSet$dataName;
  }

  dataSet <- readDataset(dataName);

  if (is.null(dataSet) || is.null(dataSet$comp.res.list) ||
      length(dataSet$comp.res.list) == 0) {
    AddErrMsg("No DE comparison results found.");
    return(0);
  }

  comp_list <- dataSet$comp.res.list
  comp_names <- names(comp_list)
  if (is.null(comp_names) || length(comp_names) == 0) {
    comp_names <- paste0("Comparison_", seq_along(comp_list))
  }
  if (!is.null(paramSet$upset.comp.nms) && length(paramSet$upset.comp.nms) > 0) {
    keep <- comp_names %in% paramSet$upset.comp.nms
    comp_list <- comp_list[keep]
    comp_names <- names(comp_list)
  }
  if (is.null(comp_names) || length(comp_names) == 0) {
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

  use_fdr <- if (is.null(paramSet$use.fdr)) TRUE else isTRUE(paramSet$use.fdr)

  sel.dats <- list()
  for (i in seq_along(comp_list)) {
    rt <- comp_list[[i]]
    if (is.null(rt) || nrow(rt) == 0) {
      sel.dats[[comp_names[[i]]]] <- character(0)
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

    sig_genes <- as.character(rownames(rt)[deg_pass & lfc_pass])
    sel.dats[[comp_names[[i]]]] <- sig_genes
  }

  sel.dats <- sel.dats[sapply(sel.dats, length) > 0]
  if (length(sel.dats) == 0) {
    AddErrMsg("No significant features for any comparison.");
    return(0);
  }

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
  paramSet <- readSet(paramSet, "paramSet");
  if (!exists("comp.vec")) {
    paramSet$upset.comp.nms <- NULL
    saveSet(paramSet, "paramSet")
    return(0)
  }
  paramSet$upset.comp.nms <- comp.vec
  rm("comp.vec", envir = .GlobalEnv)
  saveSet(paramSet, "paramSet")
  return(1)
}

#Record upset intersection mode for report
SetUpsetMode <- function(mode){
      paramSet <- readSet(paramSet, "paramSet");
      paramSet$upsetMode <- mode;
  saveSet(paramSet, "paramSet");
}


#' Plot a static UpSet PNG summarizing the pairwise sig-gene intersections.
#'
#' For multi-group / multi-contrast designs (anal.type == "default"), the
#' analysis produces N pairwise comparisons stored in
#' \code{dataSet$comp.res.list}. A single omnibus F-test answers "does
#' expression vary across the design?" but doesn't reveal WHICH dose
#' pairs each gene responds to. An UpSet plot of the per-comparison
#' significant-gene sets surfaces this cross-contrast pattern at a
#' glance — e.g., genes that respond only at high doses vs.
#' monotonically across doses vs. specific dose-pair-specific.
#'
#' @param dataName  ExpressAnalyst per-dataset key (e.g. "bmd_BRBZ_2w_entrez.txt").
#' @param imgName   Output file basename (without extension).
#' @param p.lvl     P-value cutoff for sig (matches GetSigGenes).
#' @param fc.lvl    Absolute logFC cutoff (matches GetSigGenes).
#' @param FDR       "true" → filter on adj.P.Val; "false" → raw P.Value.
#' @param max.sets  Cap on number of contrasts shown (UpSetR \code{nsets}).
#'                  With 6 dose levels, all 15 pairwise contrasts fit
#'                  easily; tighten this for designs with >8 groups.
#' @param max.inter Cap on number of intersection bars shown
#'                  (UpSetR \code{nintersects}); top-N by size.
#' @param dpi,format passed to Cairo.
#'
#' @return Path to the written PNG, or NULL on failure / no data.
PlotPairwiseUpsetPNG <- function(dataName  = "",
                                 imgName   = "upset_pairwise_0",
                                 p.lvl     = 0.05,
                                 fc.lvl    = 0,
                                 FDR       = "true",
                                 max.sets  = 5,
                                 max.inter = 20,
                                 dpi       = 150,
                                 format    = "png") {
  tryCatch({
    if (!requireNamespace("UpSetR", quietly = TRUE)) {
      message("[PlotPairwiseUpsetPNG] UpSetR package not installed — skipping")
      return(invisible(NULL))
    }
    require(Cairo)

    # Remove any stale plot from a prior run first, so every skip path below
    # leaves NO misleading old PNG for the manifest emitter to pick up.
    stale.img <- paste0(imgName, ".", format)
    if (file.exists(stale.img)) unlink(stale.img)

    dataSet <- readDataset(dataName)
    if (is.null(dataSet) || is.null(dataSet$comp.res.list) ||
        length(dataSet$comp.res.list) < 2) {
      message("[PlotPairwiseUpsetPNG] comp.res.list missing or only 1 contrast — skipping")
      return(invisible(NULL))
    }

    use.fdr <- (FDR == "true") || isTRUE(FDR)
    comp.list  <- dataSet$comp.res.list
    comp.names <- names(comp.list)
    if (is.null(comp.names) || any(comp.names == "")) {
      comp.names <- paste0("C", seq_along(comp.list))
    }

    # Build named list of significant gene IDs per contrast — UpSetR's
    # fromList() expects this shape: list(set_name = c("gene1", ...), ...)
    sig.sets <- vector("list", length(comp.list))
    names(sig.sets) <- comp.names
    for (i in seq_along(comp.list)) {
      rt <- comp.list[[i]]
      if (is.null(rt) || nrow(rt) == 0) {
        sig.sets[[i]] <- character(0); next
      }
      pcol <- if (use.fdr && "adj.P.Val" %in% names(rt)) "adj.P.Val"
              else if ("P.Value" %in% names(rt)) "P.Value"
              else if ("padj" %in% names(rt)) "padj"
              else "P.Value"
      lfc_col <- if ("logFC" %in% names(rt)) "logFC"
                 else if ("log2FoldChange" %in% names(rt)) "log2FoldChange"
                 else NULL
      pvec <- suppressWarnings(as.numeric(rt[[pcol]]))
      pass.p <- is.finite(pvec) & (pvec <= p.lvl)
      pass.fc <- if (!is.null(lfc_col)) {
        lf <- suppressWarnings(as.numeric(rt[[lfc_col]]))
        is.finite(lf) & (abs(lf) >= fc.lvl)
      } else TRUE
      sig.sets[[i]] <- rownames(rt)[pass.p & pass.fc]
    }

    # Drop contrasts with zero sig features — UpSetR fromList chokes on
    # empty sets and the resulting plot is just visual noise anyway.
    keep <- vapply(sig.sets, function(x) length(x) > 0, logical(1))
    if (sum(keep) < 2) {
      message("[PlotPairwiseUpsetPNG] fewer than 2 contrasts have sig genes — skipping")
      return(invisible(NULL))
    }
    sig.sets <- sig.sets[keep]

    # Cap to the top `max.sets` contrasts by number of significant features. A
    # many-group design (e.g. 8 groups -> 28 pairwise contrasts) otherwise yields
    # dozens of squeezed, near-empty columns; keep the most informative ones.
    if (length(sig.sets) > max.sets) {
      ord     <- order(vapply(sig.sets, length, integer(1)), decreasing = TRUE)
      dropped <- length(sig.sets) - max.sets
      sig.sets <- sig.sets[ord[seq_len(max.sets)]]
      message("[PlotPairwiseUpsetPNG] ", dropped, " smaller contrast(s) hidden — ",
              "showing top ", max.sets, " by significant-feature count")
    }

    # Build input matrix. UpSetR::fromList returns a binary
    # data.frame (genes × contrasts).
    upset.df <- UpSetR::fromList(sig.sets)
    if (nrow(upset.df) == 0) {
      message("[PlotPairwiseUpsetPNG] no genes in any sig set — skipping")
      return(invisible(NULL))
    }
    # Degenerate guard: if no feature is shared across >=2 contrasts (every
    # significant feature is unique to a single contrast) the UpSet has no real
    # intersections to show and is misleading. Skip drawing; WfEmitUpsetManifest
    # surfaces a short "no shared features" note in its place.
    if (max(rowSums(upset.df > 0), 0L) < 2) {
      message("[PlotPairwiseUpsetPNG] no feature shared across >=2 contrasts — ",
              "skipping degenerate UpSet")
      return(invisible(NULL))
    }
    # Method-standard invariant: persist the data behind the UpSet plot (the binary
    # gene x contrast membership matrix) so the figure can be re-plotted / regenerated
    # in any tool. Helper lives in wf_method.R; guard so the public
    # package still runs standalone.
    if (exists("WfSaveFigureData"))
      tryCatch(WfSaveFigureData("upset", upset.df), error = function(e) NULL)

    # Sizing — width scales with N contrasts, height with N intersections.
    n.sets <- ncol(upset.df)
    w.in <- max(8, 0.45 * n.sets + 6)
    h.in <- 6

    imgPath <- paste0(imgName, ".", format)
    Cairo::Cairo(file = imgPath, unit = "in", dpi = dpi,
                 width = w.in, height = h.in, type = format, bg = "white")
    tryCatch({
      # UpSetR::upset returns a grid object; print it to render under Cairo.
      print(UpSetR::upset(
        upset.df,
        sets         = colnames(upset.df),
        nsets        = n.sets,
        nintersects  = max.inter,
        order.by     = "freq",
        keep.order   = FALSE,
        mainbar.y.label = "Genes in intersection",
        sets.x.label    = "Sig. genes per contrast",
        point.size   = 2.4,
        line.size    = 0.8,
        text.scale   = c(1.3, 1.2, 1.2, 1.0, 1.1, 1.0)
        # text.scale order:
        #   1: intersection size title
        #   2: intersection size tick labels
        #   3: set size title
        #   4: set size tick labels
        #   5: set names
        #   6: numbers above bars
      ))
    }, error = function(e) {
      message("[PlotPairwiseUpsetPNG] UpSetR::upset render error: ", conditionMessage(e))
    })
    dev.off()

    n.total <- nrow(upset.df)
    message("[PlotPairwiseUpsetPNG] wrote ", imgPath,
            " (", length(sig.sets), " contrasts, ",
            n.total, " unique sig genes across the union)")
    return(invisible(imgPath))
  }, error = function(e) {
    message("[PlotPairwiseUpsetPNG] FAILED: ", conditionMessage(e))
    tryCatch(dev.off(), error = function(x) {})
    return(invisible(NULL))
  })
}