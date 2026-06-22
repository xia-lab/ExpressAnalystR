##################################################
## R script for ExpressAnalyst
## Description: Compute Ridgeline plot
## Authors: 
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
## Jessica Ewald, jessica.ewald@mail.mcgill.ca
##################################################

compute.ridgeline <- function(dataSet, imgNm = "abc", dpi=default.dpi, format="png", fun.type = "kegg", ridgeType = "ora", ridgeColor = "teal",rankOpt="fc", sigLevel = 0.05, pwNum=20, inx = 1){
  
  #save.image("ridge.RData");
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");
  analSet <- readSet(analSet, "analSet");
  imageName <- paste0(imgNm, "dpi" , dpi, ".", format);
  anal.type <- paramSet$anal.type;
  require("dplyr");
  suppressPackageStartupMessages({
    require("reshape");
    require("ggplot2");
    require("ggridges");
  })
  # process colors
  if(ridgeColor == "teal"){
    high.col = "#3C7C60";
    low.col = "#C9E1D6";
  } else if (ridgeColor == "orange"){
    high.col = "#ffa34d";
    low.col = "#994a00";    
  } else {
    high.col = "#00da00";
    low.col = "#005000";   
  }
  
  #get pw library
  setres <- .loadEnrichLib(fun.type, paramSet);
  current.geneset <- setres$current.geneset;
  # get DEGs
  if(anal.type == "genelist"){
    if(paramSet$numOfLists > 1){
      dataSet <- readDataset(paramSet$selDataNm);    
    }
    sigmat <- as.data.frame(dataSet$prot.mat)
    sigmat$entrez <- rownames(sigmat);
    universe <- unique(unlist(current.geneset)); 
    expr.vec <- sigmat[,1];
    
    if(sum(expr.vec) == 0){
      msgSet$current.msg <- "Uploaded gene list needs to contain fold-change value to perform Ridgeline analysis!";
      saveSet(msgSet, "msgSet");
      return(-2);
    }
  }else if(anal.type == "onedata"){
    if(ridgeType == "ora"){
      sigmat <- dataSet$sig.mat;
    }else{
      sigmat <- dataSet$comp.res;
      allmat <- dataSet$comp.res;
    }
    sigmat$entrez <- rownames(sigmat);
    universe <- rownames(dataSet$data.norm);
  }else{
    meta.avgFC <- analSet$meta.avgFC;
    inx <- 1;
    if(paramSet$selDataNm == "meta_default"){
      if(ridgeType == "ora"){
        sigmat <- analSet$meta.mat;
        allmat <- ov_qs_read("meta.resTable.qs");
        sigmat <- cbind(unname(meta.avgFC[rownames(sigmat)]), sigmat);
        
      }else{
        allmat <- ov_qs_read("meta.resTable.qs");
        sigmat <- allmat;
        sigmat <- cbind(unname(meta.avgFC[rownames(sigmat)]), sigmat);
        
      }
      allmat$logFC <- unname(meta.avgFC[rownames(allmat)]);
      universe <- rownames(allmat);
    }else{
      dataSet <- readDataset(paramSet$selDataNm);
      if(ridgeType == "ora"){
        sigmat <- dataSet$sig.mat;
      }else{
        sigmat <- dataSet$comp.res;
        allmat <- dataSet$comp.res;
      }
      sigmat <- as.data.frame(sigmat);
      sigmat$entrez <- rownames(sigmat);
      universe <- rownames(dataSet$data.norm);
    }
  }

  # Bail out cleanly if there are no significant features — otherwise the
  # downstream data.frame(rownames(sigmat), sigmat[,inx]) and reshape::melt
  # calls fire mismatched-row-count errors. Same shape as the volcano /
  # enrichment guard: write a clear message and return 0.
  if(is.null(sigmat) || nrow(sigmat) == 0){
    msgSet$current.msg <- "No significant features available for the Ridgeline plot — adjust thresholds or check earlier steps.";
    saveSet(msgSet, "msgSet");
    return(0);
  }

  if(ridgeType == "ora"){
    gene.vec <- rownames(sigmat);
    if(is.null(gene.vec) || length(gene.vec) == 0){
      msgSet$current.msg <- "No significant features for ridgeline enrichment.";
      saveSet(msgSet, "msgSet");
      return(0);
    }
    sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
    if(!is.null(sym.vec) && length(sym.vec) == length(gene.vec)){
      names(gene.vec) <- sym.vec;
    }
    .performEnrichAnalysis(dataSet, imgNm, fun.type, gene.vec, "ridgeline")
    res <- ov_qs_read("enr.mat.qs");
    # Defensive: a stale enr.mat.qs written before the drop=FALSE fix in
    # .performEnrichAnalysis can round-trip as a length-5 vector (the single
    # surviving pathway row was dropped to a vector). Coerce it back to a 1x5
    # matrix so colnames()<- below does not error with
    # "length of 'dimnames' [2] not equal to array extent".
    if (is.null(dim(res)) && length(res) == 5L) {
      res <- matrix(res, nrow = 1L);
    }
    colnames(res) <- c("size", "expected", "overlap", "pval", "padj");

    res <- res[,c(4,5,3,1,2)]
    res <- as.data.frame(cbind(pathway=rownames(res), res));
    res$padj <- as.numeric(res$padj)
    res$pval <- as.numeric(res$pval)
  } else {

    rankedVec<- ComputeRankedVec(dataSet, rankOpt, paramSet$selectedFactorInx);

    # Guard against NULL / empty ranked vector — fgsea would otherwise
    # blow up with "array(x, ...) 'data' must be of a vector type, was
    # 'NULL'" when DE results / annotation are missing.
    if (is.null(rankedVec) || length(rankedVec) == 0
        || all(is.na(rankedVec)) || !any(is.finite(rankedVec))) {
        message("[ridgeline] no ranked features available — skipping fgsea (DE results / annotation missing).");
        return(0);
    }

  gene.vec <- universe;
  sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
  gene.nms <- sym.vec;

  current.geneset.symb <- lapply(current.geneset, 
                       function(x) {
                         gene.nms[gene.vec%in%unlist(x)];
  }
  );

names(rankedVec) <- doEntrez2SymbolMapping(names(rankedVec), paramSet$data.org, paramSet$data.idType);


    use_nperm <- fun.type %in% c("go_bp", "go_mf", "go_cc")
    bridge_in_rl <- ov_bridge_file("in")
    bridge_out_rl <- sub("_in.qs2", "_out.qs2", bridge_in_rl)
    ov_qs_save(list(geneset = current.geneset.symb, ranked = rankedVec, use_nperm = use_nperm),
              bridge_in_rl, preset = "fast")
    on.exit(unlink(c(bridge_in_rl, bridge_out_rl)), add = TRUE)

    run_func_via_rsclient(
      func = function(wd, bridge_in, bridge_out) {
        setwd(wd)
        require(fgsea)
        set.seed(123)
        input <- ov_qs_read(bridge_in)
        result <- if (input$use_nperm) {
          fgsea::fgsea(pathways = input$geneset, stats = input$ranked,
                       minSize = 5, maxSize = 500, scoreType = "std", nperm = 10000)
        } else {
          fgsea::fgsea(pathways = input$geneset, stats = input$ranked,
                       minSize = 5, maxSize = 500, scoreType = "std")
        }
        ov_qs_save(result, bridge_out, preset = "fast")
      },
      args = list(wd = getwd(), bridge_in = bridge_in_rl, bridge_out = bridge_out_rl),
      timeout_sec = 600
    )

    response <- if (file.exists(bridge_out_rl)) ov_qs_read(bridge_out_rl) else NULL
    if (is.null(response)) { msgSet$current.msg <- "fGSEA ridgeline analysis failed in child process"; saveSet(msgSet, "msgSet"); return(0) }
    res <- response

  }
  res <- .signif_df(res, 4);
  res <- res[order(res$pval),];
  resTable <- res;
  # process results;
  res <- res[,c(1,2,3)];
  colnames(res) <- c("name", "pval", "adj.pval");
  res.sig <- res;
  totalSigPws <- dim(res[res$pval < sigLevel, ])[1];
  
  if(pwNum != -1){
    if(dim(res.sig)[1] > pwNum){ # limit size if too many sig results
      res.sig <- res.sig[1:pwNum, ]
    }
  }
  
  # Ridge density + jittered points should show the FULL pathway membership
  # (every measured member gene's fold change) with the DE-significant members
  # HIGHLIGHTED, mirroring the GSEA ridgeline. ORA's `sigmat` is the significant
  # subset used as the enrichment INPUT; for the PLOT use all genes' fold changes
  # (comp.res) so non-significant members aren't dropped by the na.omit below.
  # Only onedata has a separate sig.mat; gene-list / meta flows plot their set.
  sig.entrez <- if (anal.type == "onedata" && !is.null(dataSet$sig.mat))
                  rownames(dataSet$sig.mat) else rownames(sigmat)
  if (anal.type == "onedata" && ridgeType == "ora" && !is.null(dataSet$comp.res)) {
    plotmat <- as.data.frame(dataSet$comp.res); plotmat$entrez <- rownames(plotmat);
  } else {
    plotmat <- sigmat;
  }

  # prepare data for plotting
  degs.plot <- data.frame(entrez = plotmat$entrez, log2FC = plotmat[,inx]);
  degs.plot <- suppressMessages(reshape::melt(degs.plot));
  colnames(degs.plot)[1] <- "entrez";

  gs.plot <- suppressMessages(reshape::melt(current.geneset));
  colnames(gs.plot) <- c("entrez", "name");
  
  df <- merge(res.sig, gs.plot, by = "name", all.x = TRUE, all.y = FALSE);
  df <- merge(df, degs.plot, by = "entrez", all.x = TRUE, all.y = FALSE);
  df <- na.omit(df)

  # ── Min-hits filter (a ridge is a DISTRIBUTION) ──────────────────────────
  # A gene set needs at least 3 measured members ("hits") to form a ridge;
  # fewer renders as a flat line / lone point. Keep only sets with >=3 members
  # — the same minimum-set-size rule used by ma/MSEA.
  RIDGE_MIN_HITS <- 3L
  .set.hits <- tapply(as.character(df$entrez), as.character(df$name),
                      function(x) length(unique(x)))
  .keep.sets <- names(.set.hits)[!is.na(.set.hits) & .set.hits >= RIDGE_MIN_HITS]
  df <- df[as.character(df$name) %in% .keep.sets, , drop = FALSE]
  if (nrow(df) == 0) {
    try(AddErrMsg(paste0("No gene set has at least ", RIDGE_MIN_HITS,
      " measured members, so no ridgeline distribution can be drawn.")), silent = TRUE)
    if (file.exists(jsonNm)) unlink(jsonNm)
    return(0)
  }
  df$sig <- df$entrez %in% sig.entrez;   # DE-significant members -> highlighted
  
  # calculate the mean fold change to order the pathways in the plot
  means <- aggregate(df$value, by = list(df$name), mean);
  means <- means[order(means$x, decreasing = FALSE), ];
  df$name <- factor(df$name, levels = means$Group.1);
  
  # make the plot
  rp <- ggplot(df, aes(x = value, y = name, fill = adj.pval)) +
    geom_density_ridges(
      aes(point_color = sig),
      jittered_points = TRUE, point_shape = "|", point_size = 5,
      alpha = 0.6,                       # semi-transparent ridges
      color = "white",
      scale = 1.5, rel_min_height = .02, size = 0.25,
      position = position_points_jitter(height = 0)) +
    scale_discrete_manual(aesthetics = "point_color",
                          values = c(`FALSE` = "#B0B0B0", `TRUE` = "#CB181D"),
                          name = "DE significant", labels = c(`FALSE` = "no", `TRUE` = "yes")) +
    geom_vline(xintercept = 0, color = "red") +
    scale_y_discrete(expand = c(0, 0), name = "Gene Set",
                     labels = function(x) ifelse(nchar(x) > 45L,
                                                 paste0(substr(x, 1L, 42L), "..."), x)) +
    scale_x_continuous(expand = c(0, 0), name = "log2FC") +
    scale_fill_gradient("adj. pval",
                        low = high.col, high = low.col) + 
    coord_cartesian(clip = "off") +
    theme_ridges(center = TRUE) +
    theme(legend.position = "right",
          text = element_text(size=12, color = "black"),
          axis.title = element_text(size=12, face = "bold"),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(size=12,color = "black"))
  
  Cairo::Cairo(file=imageName, width=10, height=8, type=format, bg="white", dpi=dpi, unit="in");
  suppressMessages(print(rp));
  dev.off();
  
  ##interative ridge json data
  ridge_bw <- rp$layers[[1]]$computed_stat_params$bandwidth;
  jsonNm <- paste0(imgNm, ".json");
  
  symb <- doEntrez2SymbolMapping(df$entrez, paramSet$data.org, paramSet$data.idType);
  df$symbol <- symb;
  
  data.list <- list();
  gene.list <- list();
  pval.list <- list();
  col.list <- list();
  
  for(i in 1:length(levels(df$name))){
    nm <- as.character(levels(df$name)[i]);
    data.list[[ nm ]] <- as.vector(unlist(df[which(df$name == nm), "value"]));
    gene.list[[ nm ]] <- as.vector(unlist(df[which(df$name == nm), "symbol"]));
    pval.list[[ nm ]] <- unname(unlist(res[which(res$name == nm), "pval"]));
  }
  
  minFc <- min(df$value);
  maxFc <- max(df$value);
  minPval <- min(df$pval);
  maxPval <- max(df$pval);
  
  #get hits per pathway
  hits.query <- lapply(current.geneset, 
                       function(x) {
                         x[x %in% universe];
                       }
  );
  hits.query <- hits.query[resTable$pathway];
  
  # enr result for table display
  fun.anot <- hits.query
  if(ridgeType == "ora"){
    total <- resTable[,5]; if(length(total) ==1) { total <- matrix(total) };
    fun.pval <- resTable[,"pval"]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
    fun.padj <- resTable[,"padj"]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
    if(ridgeType == "ora"){
      hit.num <- resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }else{
      hit.num <- resTable[,"size"]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }  
  }else{
    total <- as.list(resTable[,5])[[1]]; if(length(total) ==1) { total <- matrix(total) };
    fun.pval <- as.list(resTable[,"pval"])[[1]]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
    fun.padj <- as.list(resTable[,"padj"])[[1]]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
    if(ridgeType == "ora"){
      hit.num <- as.list(resTable[,4])[[1]]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }else{
      hit.num <- as.list(resTable[,"size"])[[1]]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }  
    
  }
  fun.ids <- as.vector(setres$current.setids[names(fun.anot)]); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  enr.res <- list(
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total= total
  );
  
  if(ridgeType == "gsea"){
    enr.res[["ES"]] <- unname(unlist(resTable[,"ES"]));
  }
  
  res.list <- list(data=data.list, 
                   genelist=gene.list, 
                   df=df, minPval=minPval, 
                   maxPval=maxPval, 
                   min=minFc, 
                   max=maxFc, 
                   minPval = min(res.sig$pval),
                   maxPval = max(res.sig$pval),
                   bandwidth=ridge_bw,
                   pathwayPvals = pval.list, 
                   pathwayCols = col.list, 
                   enrRes = enr.res,
                   dat.opt = paramSet$selDataNm,
                   naviString="ridge");

  if(ridgeType == "gsea"){
  csv.nm <- paste0(imgNm, ".csv");
  fast.write(resTable, file=csv.nm);
  }

  analSet$ridgeline <- res.list;
  saveSet(analSet, "analSet");
  
  json.obj <- rjson::toJSON(res.list);
  sink(jsonNm);
  cat(json.obj);
  sink();
  
  #for link sharing
  paramSet$jsonNms$ridge <- jsonNm
  paramSet$partialToBeSaved <- c( paramSet$partialToBeSaved, c(jsonNm));
  saveSet(paramSet, "paramSet");
  
  imgSet <- readSet(imgSet, "imgSet");
  rownames(resTable) <- NULL;
  
  imgSet$compute.ridgeline <- imageName;
    saveSet(imgSet, "imgSet");
  
  return(totalSigPws)
}