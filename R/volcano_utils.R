  ##################################################
## R scripts for ExpressAnalyst
## Functions related to volcano plot
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#' Read volcano data from Arrow file (memory-efficient, Java-compatible)
#' @description Loads volcano data from volcano_data.arrow for zero-copy Java access
#' @return volcano list object reconstructed from Arrow data
#' @export
readVolcano <- function() {
  arrow_file <- "volcano_data.arrow"
  meta_file <- "volcano_meta.qs"

  if (file.exists(arrow_file)) {
    # Read main data from Arrow
    df <- arrow::read_feather(arrow_file)

    # Reconstruct volcano list from data frame
    vcn <- list(
      fc.symb = df$gene_id,
      fc.log = df$fc_log,
      fc.log.uniq = df$fc_log_uniq,
      p.log = df$p_log,
      p.raw = df$p_raw,
      inx.up = as.logical(df$inx_up),
      inx.down = as.logical(df$inx_down),
      inx.p = as.logical(df$inx_p)
    )
    names(vcn$fc.log) <- df$gene_id

    # Load metadata (small, ~1KB)
    if (file.exists(meta_file)) {
      meta <- ov_qs_read(meta_file)
      vcn <- c(vcn, meta)
    }
    return(vcn)
  }

  # Fallback: check analSet (legacy compatibility)
  analSet <- readSet(analSet, "analSet")
  if (!is.null(analSet$volcano)) {
    return(analSet$volcano)
  }
  return(NULL)
}

#' Save volcano data to Arrow format (Java-compatible)
#' @description Saves volcano data to Arrow for zero-copy Java access
#' @param volcano The volcano list object
#' @export
saveVolcanoArrow <- function(volcano) {
  # Create data frame with main vectors (bulk data ~11MB)
  df <- data.frame(
    gene_id = as.character(volcano$fc.symb),
    fc_log = as.numeric(volcano$fc.log),
    fc_log_uniq = as.numeric(volcano$fc.log.uniq),
    p_log = as.numeric(volcano$p.log),
    p_raw = as.numeric(volcano$p.raw),
    inx_up = as.integer(volcano$inx.up),
    inx_down = as.integer(volcano$inx.down),
    inx_p = as.integer(volcano$inx.p),
    stringsAsFactors = FALSE
  )

  # Add gene annotation if available
  if (!is.null(volcano$conv)) {
    df$symbol <- as.character(volcano$conv$symbol)
    df$gene_name <- as.character(volcano$conv$name)
  }

  # Save main data to Arrow (zero-copy for Java)
  arrow_path <- "volcano_data.arrow"
  if (file.exists(arrow_path)) {
    unlink(arrow_path)
    Sys.sleep(0.01)
  }
  arrow::write_feather(df, arrow_path, compression = "uncompressed")

  # Save small metadata to qs (~1KB)
  meta <- list(
    raw.threshx = volcano$raw.threshx,
    raw.threshy = volcano$raw.threshy,
    paired = volcano$paired,
    thresh.y = volcano$thresh.y,
    sig.mat = volcano$sig.mat,
    analType = volcano$analType,
    org = volcano$org,
    dat.opt = volcano$dat.opt,
    naviString = volcano$naviString
  )
  ov_qs_save(meta, "volcano_meta.qs")

  return(arrow_path)
}

#'Prepare data for volcano plot visualization
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
Volcano.Anal <- function(dataName="", fileNm="name", paired=FALSE, fcthresh=0, threshp=0.05, analType="NA", inx=1, dpi=default.dpi, format="png"){
  #save.image('volc.RData');

  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  anal.type <- paramSet$anal.type;

  inx <- 1
  #print("Prepare volcano anal");
  limit_fc <- T; #whether limit to -10 to 10 fc
  if(anal.type == "metadata"){
    if(paramSet$selDataNm == "meta_default"){
       if(is.null(paramSet$fc.thresh)){
         paramSet$fc.thresh <- 0; 
       }

      data <- ov_qs_read("allMeta.mat.qs")    
      p.value <- data[, 2]
      data <- cbind(unname(analSet$meta.avgFC[rownames(data)]), data);
      fcthresh <- paramSet$fc.thresh;
      threshp <- paramSet$BHth;
      inx <- 1;
    }else{
      dataSet <- readDataset(paramSet$selDataNm);
      
      data <- as.matrix(analSet$inmex.ind[paramSet$selDataNm][[1]])
      p.value <- data[, "Pval"]
    }
    limit_fc <-F;
  }else{
    dataSet <- readDataset(dataName);
    data <- as.matrix(dataSet$comp.res);
    if(is.null(paramSet$use.fdr) || paramSet$use.fdr){
        p.value <- data[, "adj.P.Val"];
    }else{
        p.value <- data[, "P.Value"];
    }
  }

  paramSet$fcthreshu <- fcthresh
  inx.p <- p.value <= threshp;
  
  zero.inx <- unname(p.value) == 0;
  p.value[zero.inx] <- min(p.value[!zero.inx])/10;
 
  p.log <- -log10(p.value);
  if(paramSet$data.org == "omk"){
    anot.id <- rownames(data);
    gene.anot <- data.frame(gene_id=anot.id, symbol=anot.id, name=anot.id, stringsAsFactors=FALSE)
  }else if (anal.type == "metadata" || isTRUE(dataSet$annotated) ){ # annotated to entrez
    anot.id <- rownames(data);
    gene.anot <- doEntrezIDAnot(anot.id, paramSet$data.org, paramSet$data.idType);
  }else{
    anot.id <- rownames(data);
    gene.anot <- data.frame(gene_id=anot.id, symbol=anot.id, name=anot.id, stringsAsFactors=FALSE)
    paramSet$init.lib <- "NA"
  }
  
  #gene symbol to be used for boxplot

  # ── Effect-size metric for the X-axis ─────────────────────────────
  # Detect an omnibus F-test result (multi-group ANOVA produced by
  # .perform_limma_edger when ncol(contrast.matrix) > 1 and
  # comp.type == "default"). In that case `data` has columns
  # `logFC.1, ..., logFC.N, AveExpr, F, P.Value, adj.P.Val` and `data[, inx]`
  # (inx=1) would read ONLY the first pairwise contrast's logFC — making
  # the volcano look like a "0 vs 25"-only plot even though the analysis
  # is a 6-way ANOVA. Replace with signed max-|logFC| across all
  # pairwise contrasts: each gene's X coordinate is the strongest
  # (largest absolute) pairwise logFC, sign preserved. This gives a
  # proper two-wing volcano where the "up" / "down" coloring reflects
  # the direction of the dominant effect across the design.
  fc.cols <- grep("^logFC", colnames(data), value = TRUE)
  is.omnibus <- length(fc.cols) > 1 && "F" %in% colnames(data)
  if (is.omnibus) {
    fc.mat  <- as.matrix(data[, fc.cols, drop = FALSE])
    storage.mode(fc.mat) <- "double"
    abs.mat <- abs(fc.mat)
    max.idx <- max.col(abs.mat, ties.method = "first")
    fc.log  <- fc.mat[cbind(seq_len(nrow(fc.mat)), max.idx)]
    names(fc.log) <- rownames(data)
  } else {
    fc.log <- data[, inx];
  }
  if(limit_fc){
    hit.maxPos <- (which(fc.log> 10) )
    hit.maxNeg <- (which(fc.log< -10) )
    fc.log[hit.maxPos] <- 10;
    fc.log[hit.maxNeg] <- -10;
  }
  #fc.all <- res$fc.all;
  
  if(fcthresh != 0){
    inx.up <- fc.log > fcthresh & p.value < threshp;
    inx.down <- fc.log < -fcthresh & p.value < threshp;
  }else{
    inx.up <- fc.log > 0 & p.value < threshp;
    inx.down <- fc.log < 0 & p.value < threshp;
  }
  
  # create named sig table for display
  inx.imp <- (inx.up | inx.down) & inx.p;
  sig.var <- cbind(fc.log[inx.imp,drop=F], p.value[inx.imp,drop=F], p.log[inx.imp,drop=F]);
  colnames(sig.var) <- c("log2(FC)", "p.value", "-log10(p)");
  # first order by log(p), then by log(FC)
  ord.inx <- order(sig.var[,3], abs(sig.var[,1]), decreasing=T);
  sig.var <- sig.var[ord.inx,, drop=F];
  
  sig.var <- signif (sig.var, 5);
  sig.var1 <- sig.var;
  sig.var1 <- cbind(rownames(sig.var), sig.var);
  colnames(sig.var1) <- c("name", "log2(FC)", "p.value", "-log10(p)");
  
  ###########################
  ## for Volcano data
  ##########################
  
  if(paramSet$init.lib != "NA"){
    saveSet(paramSet, "paramSet");
    PerformVolcanoEnrichment(dataName, "enrichment_result", paramSet$init.lib, "null", "all", inx)
    paramSet <- readSet(paramSet, "paramSet");
    msgSet <- readSet(msgSet, "msgSet");
  }
  
  fileName <- "volcano.csv";
  jsonNm <- "volcano.json";
  json.obj <- rjson::toJSON(sig.var1);
  sink(jsonNm);
  cat(json.obj);
  sink();
  fast.write(signif (sig.var,5),file=fileName);
  colnames(gene.anot)[1] <- "anot.id"
  
  volcano <- list (
    raw.threshx = fcthresh,
    raw.threshy = threshp,
    paired = paired,
    thresh.y = -log10(threshp),
    fc.symb =rownames(data),
    fc.log = fc.log,
    fc.log.uniq = jitter(fc.log),
    inx.up = inx.up,
    inx.down = inx.down,
    p.log = p.log,
    p.raw = p.value,
    inx.p = inx.p,
    sig.mat = sig.var,
    conv = gene.anot,
    analType = anal.type,
    org=paramSet$data.org,
    dat.opt = paramSet$selDataNm,
    naviString = "Volcano Plot"
  );
  
  # Save volcano to Arrow for zero-copy Java access (saves ~11 MB memory)
  saveVolcanoArrow(volcano)
  analSet$volcano <- NULL  # Clear from memory, use readVolcano() to access
  saveSet(analSet, "analSet");

  # Get IDs (these read from volcano_data.arrow)
  sigDownIds <- GetVolcanoUpLftIDs();
  sigUpIds <- GetVolcanoUpRgtIDs();
  nonSigIds <- GetVolcanoDnIDs();

  # Add IDs to local volcano object for JSON output
  volcano[["sigDownIds"]] <- sigDownIds;
  volcano[["sigUpIds"]] <- sigUpIds;
  volcano[["nonSigIds"]] <- nonSigIds;
  
  jsonNm <- paste0(fileNm, ".json");
  json.obj <- rjson::toJSON(volcano);
  sink(jsonNm);
  cat(json.obj);
  sink();
  
  if(paramSet$init.lib == "NA"){
    enr.mat <- "NA"
  }else{
    enr.mat <- ov_qs_read("enr.mat.qs");
    #fast.write(enr.mat, file="enrichment_result.csv", row.names=T);
  }
  sink("enrichment_result.json");
  cat(json.obj);
  sink();
  paramSet$jsonNms["volcano"] <- fileNm;

  # Generate volcano data frame
  # Detect x-axis semantics:
  #   is.omnibus → multi-group ANOVA omnibus: fc.log is the signed
  #     max-|logFC| across all pairwise contrasts (computed above);
  #     X-axis label communicates this.
  #   is.fstat  → kept for forward compat if someone wires the F-stat
  #     column directly as the X axis (current code doesn't); label
  #     "F statistic", suppress symmetric vertical guides.
  #   else      → standard two-group / pairwise logFC.
  xCol <- if (is.null(colnames(data))) "logFC" else colnames(data)[inx]
  is.fstat <- !is.null(xCol) && grepl("^F([. ]|stat|val|$)", xCol, ignore.case = TRUE)

  volcano_data <- data.frame(
    gene = gene.anot$symbol,
    log2FoldChange = fc.log,
    pValue = p.value,
    negLog10PValue = p.log,
    significant = ifelse(inx.up & inx.p, "Up",
                         ifelse(inx.down & inx.p, "Down", "Non-sig")),
    stringsAsFactors = FALSE
  )

  # Hover text retained for plotly interactive widget
  volcano_data$hover_text <- with(volcano_data, paste("Gene:", gene,
                                                      "<br>Log2 FC:", signif(log2FoldChange, 3),
                                                      "<br>P-value:", format(pValue, scientific = TRUE, digits = 3)))

  # Label the top-N most significant up/down hits (by -log10(p), |logFC| as tiebreaker).
  labelNum <- 8
  volcano_data$label <- NA_character_

  up_idx <- which(volcano_data$significant == "Up")
  if(length(up_idx) > 0) {
    up_ord <- up_idx[order(volcano_data$negLog10PValue[up_idx],
                           abs(volcano_data$log2FoldChange[up_idx]),
                           decreasing = TRUE)]
    top_up <- head(up_ord, labelNum)
    volcano_data$label[top_up] <- volcano_data$gene[top_up]
  }

  down_idx <- which(volcano_data$significant == "Down")
  if(length(down_idx) > 0) {
    down_ord <- down_idx[order(volcano_data$negLog10PValue[down_idx],
                               abs(volcano_data$log2FoldChange[down_idx]),
                               decreasing = TRUE)]
    top_down <- head(down_ord, labelNum)
    volcano_data$label[top_down] <- volcano_data$gene[top_down]
  }

  # Stable factor order so the legend always reads Up / Down / Non-sig
  volcano_data$significant <- factor(volcano_data$significant, levels = c("Up", "Down", "Non-sig"))

  # Hit counts annotate the title for quick scanning
  nUp   <- sum(volcano_data$significant == "Up",   na.rm = TRUE)
  nDown <- sum(volcano_data$significant == "Down", na.rm = TRUE)
  nTot  <- nrow(volcano_data)

  library(ggplot2)
  library(plotly)

  # Categorical palette: warm red for up, cool blue for down, mid-grey for the
  # non-significant cloud. Sourced from the same RdBu/diverging family used by
  # MetaboAnalyst's volcano gradient endpoints.
  vc_cols <- c("Up" = "#d6604d", "Down" = "#2166ac", "Non-sig" = "#9ca3af")

  x_lab <- if (is.fstat) {
             expression(italic(F)~statistic)
           } else if (is.omnibus) {
             expression(log[2]~"(Max-pairwise Fold Change)")
           } else {
             expression(log[2]~"(Fold Change)")
           }
  y_lab <- expression(-log[10]~"(P-value)")

  # Layer non-sig points first, then up/down on top so significant hits aren't buried.
  ns_df  <- volcano_data[volcano_data$significant == "Non-sig", , drop = FALSE]
  sig_df <- volcano_data[volcano_data$significant != "Non-sig", , drop = FALSE]

  gg_volcano <- ggplot() +
    geom_hline(yintercept = -log10(threshp),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    {if (!is.fstat && fcthresh > 0)
       geom_vline(xintercept = c(-fcthresh, fcthresh),
                  linetype = "dashed", color = "grey40", linewidth = 0.4)} +
    geom_point(data = ns_df,
               aes(x = log2FoldChange, y = negLog10PValue,
                   color = significant, text = hover_text),
               alpha = 0.45, size = 1.4, shape = 16) +
    geom_point(data = sig_df,
               aes(x = log2FoldChange, y = negLog10PValue,
                   color = significant, text = hover_text),
               alpha = 0.85, size = 1.9, shape = 16) +
    scale_color_manual(name = NULL, values = vc_cols, drop = FALSE) +
    labs(x = x_lab, y = y_lab,
         title = sprintf("Up: %d   Down: %d   Total: %d", nUp, nDown, nTot),
         subtitle = if (is.omnibus)
                      "Multi-group ANOVA: omnibus F-test p-value, signed max-pairwise logFC"
                    else NULL) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
          plot.title = element_text(size = 11, color = "grey30", hjust = 0),
          axis.title = element_text(size = 13),
          axis.text  = element_text(size = 11, color = "grey20"),
          legend.position = "right",
          legend.text = element_text(size = 11),
          legend.key = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  # Static (PNG) variant gets text labels for the top sig hits.
  gg_volcano_labeled <- gg_volcano +
    ggrepel::geom_text_repel(
      data = volcano_data[!is.na(volcano_data$label), , drop = FALSE],
      aes(x = log2FoldChange, y = negLog10PValue, label = label),
      size = 3.2, max.overlaps = Inf,
      box.padding = 0.45, point.padding = 0.25,
      min.segment.length = 0,
      segment.color = "grey55", segment.size = 0.3,
      color = "grey15", show.legend = FALSE)

  # Save outputs
  imgSet <- readSet(imgSet, "imgSet");
  imgSet$volcanoPlot <- paste0(fileNm, ".png");

  # Create interactive plotly only when node count <= 5000
  if (nrow(volcano_data) <= 5000) {
    pwidget <- ggplotly(gg_volcano, tooltip = "text") %>% layout(hovermode = 'closest')
    imgSet$volcanoPlotly <- paste0(fileNm, ".rda");
    save(pwidget, file = imgSet$volcanoPlotly);
  } else {
    imgSet$volcanoPlotly <- NULL;
  }

  # PNG sizing: previous 13.9 x 11.1 in @ 96 dpi produced a stretched, low-res
  # image (~1334 x 1066 px). 7 x 6 in @ dpi parameter (default ~150) yields a
  # crisper, properly-proportioned figure (~1050 x 900 px at 150 dpi).
  Cairo::Cairo(file = imgSet$volcanoPlot, unit = "in", dpi = dpi,
               width = 7, height = 6, type = format, bg = "white");
  print(gg_volcano_labeled)
  dev.off()

  saveSet(imgSet, "imgSet");
  saveSet(paramSet, "paramSet");
  saveSet(analSet, "analSet");

  return(1);
}


GetVolcanoDnMat <- function(){
  vcn <- readVolcano();
  imp.inx <- (vcn$inx.up | vcn$inx.down) & vcn$inx.p;
  blue.inx <- which(!imp.inx);
  
  if(sum(blue.inx)>0){
    xs <- vcn$fc.log[blue.inx]
    ys <- vcn$p.log[blue.inx];
    return(as.matrix(cbind(xs, ys)));
  }else{
    return(as.matrix(cbind(-1, -1)));
  }
}


GetVolcanoUpLftMat <- function(){
  vcn <- readVolcano();
  imp.inx <- vcn$inx.down & vcn$inx.p;
  red.inx <- which(imp.inx);
  if(sum(red.inx)>0){
    xs <- vcn$fc.log[red.inx]
    ys <- vcn$p.log[red.inx];
    return(as.matrix(cbind(xs, ys)));
  }else{
    return(as.matrix(cbind(-1, -1)));
  }
}

GetVolcanoUpRgtMat <- function(){
  vcn <- readVolcano();
  imp.inx <- vcn$inx.up & vcn$inx.p;
  red.inx <- which(imp.inx);
  if(sum(red.inx)>0){
    xs <- vcn$fc.log[red.inx]
    ys <- vcn$p.log[red.inx];
    return(as.matrix(cbind(xs, ys)));
  }else{
    return(as.matrix(cbind(-1, -1)));
  }
}

GetVolcanoUpLftIDs <- function(){
  vcn <- readVolcano();
  imp.inx <- vcn$inx.down & vcn$inx.p;
  red.inx <- which(imp.inx);
  if(sum(red.inx)>0){
    return(names(vcn$fc.log)[red.inx]);
  }else{
    return("NA");
  }
}

GetVolcanoUpRgtIDs <- function(){
  vcn <- readVolcano();
  imp.inx <- vcn$inx.up & vcn$inx.p;
  red.inx <- which(imp.inx);
  if(sum(red.inx)>0){
    return(names(vcn$fc.log)[red.inx]);
  }else{
    return("NA");
  }
}

GetVolcanoDnIDs <- function(){
  vcn <- readVolcano();
  imp.inx <- (vcn$inx.up | vcn$inx.down) & vcn$inx.p;
  blue.inx <- which(!imp.inx);
  if(sum(blue.inx)>0){
    return(names(vcn$fc.log)[blue.inx]);
  }else{
    return("NA");
  }
}


PerformVolcanoEnrichment<-function(dataName="", file.nm, fun.type, IDs, type, inx){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  fcthreshu <- paramSet$fcthreshu;
  fcthreshu <- as.numeric(fcthreshu);
  anal.type <- paramSet$anal.type;
  inx <- as.numeric(inx)
  if(anal.type == "onedata"){
    if(dataSet$type == "array"){
      sigmat <- dataSet$sig.mat
    } else {
      sigmat <- dataSet$sig.mat
    }
  }else{
    if(paramSet$selDataNm == "meta_default"){
      sigmat <- analSet$meta.mat
      sigmat <- cbind(unname(analSet$meta.avgFC[rownames(sigmat)]), sigmat);
      inx <- 1;
    }else{
      sigmat <- analSet$inmex.ind[paramSet$selDataNm][[1]][which(analSet$inmex.ind[paramSet$selDataNm][[1]][,'Pval'] < as.numeric(paramSet$pvalu)),,drop=FALSE];
    }
  }

  if(!is.matrix(sigmat)) sigmat <- as.matrix(sigmat);

  if(type == "focus"){
    gene.vec <- unlist(strsplit(IDs, "; "));
  }else if(nrow(sigmat) == 0){
    gene.vec <- character(0);
  }else if(type == "all"){
    gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),,drop=FALSE]);
    gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),,drop=FALSE]);
    gene.vec <- c(gene.vecup, gene.vecdown);
  }else if(type == "up"){
    gene.vec <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),,drop=FALSE]);
  }else{
    gene.vec <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),,drop=FALSE]);
  }
  sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
  names(gene.vec) <- sym.vec;
  res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, gene.vec, "volcano");
  return(res);
}


# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformVolcanoBatchEnrichment <- function(dataName="", file.nm, fun.type, IDs, inx){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  msgSet <- readSet(msgSet, "msgSet");
  anal.type <- paramSet$anal.type;
  # prepare lib
  inx <- as.numeric(inx);
  if(anal.type == "onedata"){
    if(dataSet$type == "array"){
      sigmat <- dataSet$sig.mat
    } else {
      sigmat <- dataSet$sig.mat
    }
  }else{
    sigmat <- analSet$inmex.ind[paramSet$selDataNm][[1]][which(analSet$inmex.ind[paramSet$selDataNm][[1]][,'Pval'] < as.numeric(paramSet$pvalu)),,drop=FALSE];
  }

  if(!is.matrix(sigmat)) sigmat <- as.matrix(sigmat);

  one.path.vec <- unlist(strsplit(IDs, "; "));

  if(nrow(sigmat) == 0){
    gene.vecup <- character(0);
    gene.vecdown <- character(0);
  }else{
    gene.vecup <- rownames(sigmat[which(sigmat[,inx] > paramSet$fcthreshu),,drop=FALSE]);
    gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -paramSet$fcthreshu),,drop=FALSE]);
  }
  ora.vec <- c(gene.vecup, gene.vecdown);
  
  
  sym.vec <- doEntrez2SymbolMapping(ora.vec, paramSet$data.org, paramSet$data.idType);
  names(ora.vec) <- sym.vec;
  
  current.geneset <- list()
  current.geneset[["Set"]] <- one.path.vec
  current.geneset[["Set2"]] <- one.path.vec
  
  # prepare query
  ora.nms <- names(ora.vec);
  
  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval", "FDR");
  
  # need to cut to the universe covered by the pathways, not all genes
  if(paramSet$universe.opt == "library"){
    current.universe <- unique(unlist(current.geneset));     
  }else{
    # cut to the universe to uploaded genes
    if(paramSet$anal.type == "onedata"){
      data.anot <- .get.annotated.data();
      current.universe <- rownames(data.anot); 
    }else if(paramSet$anal.type == "metadata"){
      inmex <- ov_qs_read("inmex_meta.qs");
      current.universe <- rownames(inmex$data); 
    }else{
      if(!is.null(paramSet$backgroundUniverse)){
        current.universe <- paramSet$backgroundUniverse;
      }else{
        current.universe <- unique(unlist(current.geneset)); 
      }
    }
  }
  
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];
  
  q.size<-length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.geneset, 
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  
  ov_qs_save(hits.query, "hits_query.qs");
  
  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(unique(x))}), use.names=FALSE);
  
  # total unique gene number
  #uniq.count <- length(current.universe);
  uniq.count <- nrow(dataSet$data.norm);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));
  
  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  # Replace NaN values with 1
  raw.pvals[is.nan(raw.pvals)] <- 1

  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- raw.pvals;
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  if(nrow(res.mat)> 0){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    #res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    imp.inx <- res.mat[,4] <= 0.05;
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }else{
    return(0);
  }
  
  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  res.mat[,"Hits"] <- res.mat[,"Hits"]

  # Check for and handle duplicate row names in enr.mat
  if(any(duplicated(rownames(res.mat)))) {
    res.mat <- res.mat[!duplicated(rownames(res.mat)), ]
    hits.query <- hits.query[match(rownames(res.mat), names(hits.query))]
  } else {
    res.mat <- res.mat
  }

  ov_qs_save(res.mat, "enr.mat.qs");
  msgSet$current.msg <- "Functional enrichment analysis was completed";
  
  # write json
  fun.anot <- hits.query; 
  total <- resTable$Total; if(length(total) ==1) { total <- matrix(total) };
  fun.pval <- resTable$Pval; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj <- resTable$FDR; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  hit.num <- paste0(resTable$Hits,"/",resTable$Total); if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(names(fun.anot)); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = "",
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total= total
  );
  json.mat <- rjson::toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  fun.hits <<- hits.query;
  fun.pval <<- resTable$Pval;
  hit.num <<- resTable$Hits;
  csv.nm <- paste(file.nm, ".csv", sep="");    
  fast.write(resTable, file=csv.nm, row.names=F);
  
  saveSet(msgSet, "msgSet");
  return(1);
}