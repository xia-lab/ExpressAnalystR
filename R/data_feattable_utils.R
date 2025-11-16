#' Get Significant Genes from Analysis
#'
#' This function retrieves significant genes from the DE analysis based on the specified parameters.
#'
#' @param dataName A character string specifying the name of the dataset.
#' @param res.nm A character string specifying the name of the output result table
#' @param p.lvl The significance threshold for p-values.
#' @param fc.lvl The fold change threshold.
#' @param inx The index for comparison (e.g., in case of multiple comparisons).
#'
#' @return A list of information including the filename, significant gene details, and counts.
#'
#' @author Guangyan Zhou \email{guangyan.zhou@mail.mcgill.ca}
#' @details Additional details about the function, if needed.
#'
#' @examples
#' \dontrun{
#' GetSigGenes(dataName = "MyData", res.nm = "result_A", p.lvl = 0.05,
#'             fc.lvl = 1, inx = 1)
#' }
#'
#' @export
#' @license MIT License
#'
GetSigGenes <-function(dataName="", res.nm="nm", p.lvl=0.05, fc.lvl=1, inx=1, FDR = "true"){
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");
  analSet <- readSet(analSet, "analSet");
  dataSet <- readDataset(dataName);

  paramSet$use.fdr <- as.logical(FDR);
  total <- nrow(dataSet$comp.res);
  resTable <- dataSet$comp.res;
  filename <- dataSet$filename;
  if(is.null(dataSet$fc.lvl)){
      dataSet$fc.lvl <- 0;
  }
  filename <- paste(filename, "_", res.nm, "_fc_" , dataSet$fc.lvl, ".csv", sep="");
  current.msg <- "";
  
  if (is.null(resTable) || nrow(resTable) == 0){
    msgSet$current.msg <- paste(msgSet$current.msg, "No significant genes were identified using the given design and cutoff."); 
  }
  # now rank by logFC, note, the logFC for each comparisons 
  # are returned in resTable before the AveExpr columns 
  # for two-class, only one column, multiple columns can be involved
  # for > comparisons - in this case, use the largest logFC among all comparisons
  # further filter by logFC
  if (dataSet$de.method=="deseq2"){
    hit.inx <- which(colnames(resTable) == "baseMean"); 
    dataSet$comp.res <- dataSet$comp.res.list[[inx]];
    resTable <- dataSet$comp.res;
   } else if (dataSet$de.method=="limma" || dataSet$de.method=="wtt" ){
    hit.inx <- match("AveExpr", colnames(resTable));
    dataSet$comp.res <- dataSet$comp.res.list[[inx]];
    resTable <- dataSet$comp.res;

    if (is.na(hit.inx) && dataSet$de.method == "wtt") {
      ave.expr <- rowMeans(
        dataSet$data.norm[rownames(resTable), , drop = FALSE],
        na.rm = TRUE)
      resTable$AveExpr <- ave.expr
      hit.inx <- match("AveExpr", colnames(resTable))
      dataSet$comp.res <- resTable
      dataSet$comp.res.list[[inx]] <- resTable
    }
  } else {
    hit.inx <- which(colnames(resTable) == "logCPM");
    dataSet$comp.res <- dataSet$comp.res.list[[inx]];
    resTable <- dataSet$comp.res;
  }

  if(length(hit.inx) == 0){
    hit.inx <- 1;
  }

  resTable <- resTable[!is.na(resTable[,1]),]
  orig.resTable <- resTable;
  # select based on p-value
  if(FDR == "true"){
      hit.inx.p <- resTable$adj.P.Val <= p.lvl; 
  } else {
      hit.inx.p <- resTable$P.Value <= p.lvl; 
  } 
  
  resTable<-resTable[hit.inx.p,];
  
  if (is.na(hit.inx) || hit.inx < 2) {
    maxFC.inx <- 1
  } else {
    maxFC.inx <- hit.inx - 1
  }

  if (ncol(resTable) == 0) {
    logfc.mat <- matrix(0, nrow = nrow(resTable), ncol = 1,
                        dimnames = list(rownames(resTable), "logFC"))
  } else {
    cols_to_take <- seq_len(min(maxFC.inx, ncol(resTable)))
    logfc.mat <- resTable[, cols_to_take, drop = FALSE];
  }
  if(paramSet$oneDataAnalType == "dose"){
    pos.mat <- abs(logfc.mat);
    fc.vec <- apply(pos.mat, 1, max);   # for > comparisons - in this case, use the largest logFC among all comparisons
    hit.inx.fc <- fc.vec >= fc.lvl;
    resTable <- resTable[hit.inx.fc,];
  } else if (dataSet$de.method=="deseq2" || dataSet$de.method=="edger" || dataSet$de.method=="limma" || dataSet$de.method=="wtt"){
    pos.mat <- abs(logfc.mat);
    fc.vec <- pos.mat[,1];
    hit.inx.fc <- fc.vec >= fc.lvl;
    resTable <- resTable[hit.inx.fc,];
  }else {
    pos.mat <- abs(logfc.mat[,inx]);
    fc.vec <- pos.mat;
    hit.inx.fc <- fc.vec >= fc.lvl;
    resTable <- resTable[hit.inx.fc,];
  }
  
  if (nrow(resTable) == 0){
    msgSet$current.msg <- paste(msgSet$current.msg, "No significant genes were identified using the given design and cutoff.");
  }
  
  ### Note, rowname of resTable must be entrez ID
  # calculate differently if dose-response
  de.Num <- nrow(resTable);

  non.de.Num <- nrow(dataSet$data.norm) - de.Num;
  
  # may need to update data, class and meta.info
  data <- dataSet$data.norm;
  cls <- dataSet$cls; 
  meta.info <- dataSet$meta.info;
  grp.nms <- levels(cls);
  
  hit.inx <- cls %in% grp.nms;
  if (sum(hit.inx) < length(hit.inx)){
    msgSet$current.msg <- paste(msgSet$current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse=", "), "are included.");
    cls <- factor(cls[hit.inx]);
    cls.lvls <- levels(cls);
    data <- data[,hit.inx];
    meta.info <- dataSet$meta.info[hit.inx,];
  }
  qs::qsave(data, file="data.stat.qs");
o <- with(dataSet$comp.res, order(P.Value, -abs(logFC), na.last = TRUE))
dataSet$comp.res <- dataSet$comp.res[o, , drop = FALSE]
dataSet$comp.res <- dataSet$comp.res[
                      !(rownames(dataSet$comp.res) %in% rownames(resTable)), ]
dataSet$comp.res <- rbind(resTable, dataSet$comp.res)
  
  
  dataSet$sig.mat <- resTable;
  
  if (dataSet$annotated){ # annotated to entrez
    anot.id <- rownames(dataSet$comp.res);
    gene.anot <- doEntrezIDAnot(anot.id, paramSet$data.org, paramSet$data.idType)
    fast.write(cbind(EntrezID=anot.id, signif (dataSet$comp.res,5), Symbols = gene.anot$symbol,  Name=gene.anot$name), row.names=F, file=filename);
  } else if (file.exists("annotation.qs")){ # annotation information available
    anot.id <- qs::qread("annotation.qs");
    feature.vec <- rownames(dataSet$comp.res);
    entrez.vec <- anot.id[feature.vec];
    gene.anot <- doEntrezIDAnot(entrez.vec, paramSet$data.org, paramSet$data.idType);
    fast.write(cbind(signif (dataSet$comp.res,5), EntrezID=entrez.vec, Symbols = gene.anot$symbol,  Name=gene.anot$name), row.names=F, file=filename);
    rownames(gene.anot) <- feature.vec;
  } else {
    gene.anot <- NULL;
    fast.write(signif(resTable,5), file=filename);
  }
  if(is.null(gene.anot)){
    analSet$comp.genes.symbols <- rownames(dataSet$comp.res); # use the id provided
  }else{
    analSet$comp.genes.symbols <- gene.anot$symbol;
  }
  dataSet$cls.stat <- cls;
  dataSet$meta.stat <- meta.info;
  
  # now do protein mapping for network only applicable for annotated
  
  gene <- rownames(resTable);
  
  logFC <- unname(logfc.mat[,1]);
  geneList <- paste(gene, logFC, collapse="\n");
  up <- nrow(resTable[which(logfc.mat[,paramSet$selectedFactorInx]> fc.lvl),])
  down <- nrow(resTable[which(logfc.mat[,paramSet$selectedFactorInx]< -fc.lvl),])
  saveSet(msgSet, "msgSet");
  
  data.norm <- dataSet$data.norm
  colnames(data.norm) <- NULL
  lst <- list(colnames(dataSet$data.norm),data.norm, dataSet$meta.info, dataSet$comp.res, rownames(data.norm), org=paramSet$data.org)
  json.obj <- rjson::toJSON(lst);
  sink("ExpressAnalyst_matrix.json");
  cat(json.obj);
  sink();

        if (dataSet$de.method %in% c("deseq2", "edger", "limma", "wtt")) {

  significant_gene_table <- list()    # holds one data-frame per comparison

  for (inx in seq_along(dataSet$comp.res.list)) {

    resTable <- dataSet$comp.res.list[[inx]]

    resTable <- resTable[!is.na(resTable[, 1]), ]

    deg.pass <- if (FDR == "true")  resTable$adj.P.Val <= p.lvl
                else                resTable$P.Value   <= p.lvl

    lfc.pass <- abs(resTable[ , "logFC"]) >= fc.lvl

    all_pass <- deg.pass & lfc.pass
    if (!any(all_pass)) {                            # nothing passed
      significant_gene_table[[inx]] <- data.frame() # keep list length constant
      next
    }

    sig <- resTable[all_pass, ]
    sig$GeneID      <- rownames(sig)                # preserve raw ID
    sig$Comparison  <- names(dataSet$comp.res.list)[[inx]]

    res.anot <- doEntrezIDAnot(sig$GeneID,
                               paramSet$data.org,
                               paramSet$data.idType)
    sig$Symbol <- res.anot$symbol
    sig$Name   <- res.anot$name

    significant_gene_table[[inx]] <- sig
  }

  final_table <- do.call(rbind, significant_gene_table)  # may have duplicates

  output_file <- paste0(dataName, "_logFC_",format(as.numeric(fc.lvl), digits = 2, nsmall = 0, trim = TRUE, scientific = FALSE),
                        "_Significant_Genes.csv")

  if (nrow(final_table) > 0) {
    write.csv(final_table[ , setdiff(names(final_table), "GeneID")],
              file = output_file, row.names = FALSE)

    all_significant_genes <- unique(final_table$GeneID)  # de-duplicate here
    de.Num.total          <- length(all_significant_genes)

    message("Significant genes table exported to: ", output_file)
  } else {
    de.Num.total <- 0
    message("No significant genes identified to export.")
  }

  ## ---------- bookkeeping -----------------------------------------------
  if (de.Num.total == 0) {
    msgSet$current.msg <- paste(
      msgSet$current.msg,
      "No significant genes were identified using the given design and cutoff."
    )
  }
  analSet$sig.gene.count.total <- de.Num.total


    ## 2) Build & write the binary 0/1 incidence table (one row per gene)
    comp_list  <- dataSet$comp.res.list
    comp_names <- names(comp_list)
    use_fdr    <- (FDR == "true") || isTRUE(FDR)

    # collect DE gene sets per comparison
    pass_sets <- vector("list", length(comp_list))
    names(pass_sets) <- comp_names
    for (i in seq_along(comp_list)) {
      rt <- comp_list[[i]]
      if (is.null(rt) || nrow(rt) == 0) { pass_sets[[i]] <- character(0); next }

      pcol <- if (use_fdr && "adj.P.Val" %in% names(rt)) "adj.P.Val" else "P.Value"
      lfc_col <- if ("logFC" %in% names(rt)) "logFC"
                 else if ("log2FoldChange" %in% names(rt)) "log2FoldChange"
                 else NULL

      deg_pass <- is.finite(rt[[pcol]]) & (rt[[pcol]] <= p.lvl)
      lfc_pass <- if (!is.null(lfc_col)) is.finite(rt[[lfc_col]]) & (abs(rt[[lfc_col]]) >= fc.lvl) else TRUE

      pass_sets[[i]] <- rownames(rt)[deg_pass & lfc_pass]
    }

    genes <- sort(unique(unlist(pass_sets)))

    # derive a second filename from your existing output_file
    out_bin_file <- sub("_Significant_Genes\\.csv$", "_DE_binary_matrix.csv", output_file)
    if (identical(out_bin_file, output_file)) {
      out_bin_file <- paste0(tools::file_path_sans_ext(output_file), "_DE_binary_matrix.csv")
    }

    if (length(genes) > 0) {
      de_mat <- matrix(0L, nrow = length(genes), ncol = length(comp_names),
                       dimnames = list(genes, comp_names))
      for (i in seq_along(pass_sets)) {
        if (length(pass_sets[[i]]) > 0) de_mat[pass_sets[[i]], i] <- 1L
      }

      # optional annotation; if it fails, keep Symbol/Name as NA
      annot <- tryCatch(
        doEntrezIDAnot(genes, paramSet$data.org, paramSet$data.idType),
        error = function(e) data.frame(symbol = rep(NA_character_, length(genes)),
                                       name   = rep(NA_character_, length(genes)))
      )

      de_df <- data.frame(
        GeneID = genes,
        Symbol = annot$symbol,
        Name   = annot$name,
        as.data.frame(de_mat, check.names = FALSE),
        check.names = FALSE
      )
      de_df$DE_Count <- rowSums(de_mat)
      de_df$Comparisons <- apply(de_mat, 1, function(z) {
        hits <- which(z == 1L)
        if (length(hits) == 0) "" else paste(comp_names[hits], collapse = ";")
      })

      write.csv(de_df, file = out_bin_file, row.names = FALSE)
      message("Binary DE matrix written to: ", out_bin_file)
    } else {
      # write an empty shell (optional)
      write.csv(data.frame(GeneID = character(0)), file = out_bin_file, row.names = FALSE)
      message("No genes passed any comparison; wrote empty binary matrix to: ", out_bin_file)
    }
}
dataSet$comp.res.list <- lapply(dataSet$comp.res.list, function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return(tbl)

  pcol <- if ( "P.Value" %in% names(tbl)) "P.Value"
          else if ("padj" %in% names(tbl)) "padj"
          else "P.Value"
  pvec <- suppressWarnings(as.numeric(tbl[[pcol]]))

  # tie-break by |logFC|
  if ("logFC" %in% names(tbl)) {
    lfc <- abs(suppressWarnings(as.numeric(tbl$logFC)))
  } else if ("log2FoldChange" %in% names(tbl)) {
    lfc <- abs(suppressWarnings(as.numeric(tbl$log2FoldChange)))
  } else {
    lfc <- rep(0, nrow(tbl))
  }

  sig_idx <- (pvec <= p.lvl) & (lfc >= fc.lvl)

  sig  <- tbl[sig_idx,  , drop = FALSE]
  rest <- tbl[!sig_idx, , drop = FALSE]

  if (nrow(sig)  > 0) sig  <- sig [order(pvec[sig_idx],  -lfc[sig_idx],  na.last = NA),  , drop = FALSE]
  if (nrow(rest) > 0) rest <- rest[order(pvec[!sig_idx],            na.last = TRUE), , drop = FALSE]

  rbind(sig, rest)
})

## Rebuild the combined table with “sig first” and sorted by p
# p column for the combined table
pcol <- if (paramSet$use.fdr && "adj.P.Val" %in% names(resTable)) "adj.P.Val"
        else if ("padj" %in% names(resTable)) "padj"
        else "P.Value"

# ensure the sig block (resTable) itself is ordered
resTable <- resTable[order(as.numeric(resTable[[pcol]]),
                           -abs(suppressWarnings(as.numeric(resTable[[if ("logFC" %in% names(resTable)) "logFC" else 1]]))),
                           na.last = NA),
                     , drop = FALSE]

# the remainder (“other”), ordered too
other <- dataSet$comp.res[!(rownames(dataSet$comp.res) %in% rownames(resTable)), , drop = FALSE]
if (nrow(other) > 0 && pcol %in% names(other)) {
  # pick tie-break column for 'other'
  tie_col <- if ("logFC" %in% names(other)) "logFC" else 1
  other <- other[order(as.numeric(other[[pcol]]),
                       -abs(suppressWarnings(as.numeric(other[[tie_col]]))),
                       na.last = TRUE),
                 , drop = FALSE]
}

dataSet$comp.res <- rbind(resTable, other)

  output_file <- paste0(dataName, "_logFC_",format(as.numeric(fc.lvl), digits = 2, nsmall = 0, trim = TRUE, scientific = FALSE),
                        "_Significant_Genes.csv")
    write.csv(dataSet$comp.res,
              file = output_file, row.names = FALSE)

  analSet$sig.gene.count <- de.Num;
  saveSet(analSet, "analSet");

  dataSet$pval <- p.lvl;
  dataSet$fc.val <- fc.lvl;
  dataSet$comp.res.filename <- filename;
  res <- RegisterData(dataSet);
  saveSet(paramSet, "paramSet");
  return(c(output_file, de.Num, geneList, total, up, down, non.de.Num));
}
