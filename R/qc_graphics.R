##################################################
## R script for ExpressAnalyst
## Description: functions for quality check boxplot
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

PlotDataBox <- function(fileName, boxplotName, dpi, format){
  dataSet <- readDataset(fileName);
  qc.boxplot(as.matrix(dataSet$data.norm), boxplotName, dpi, format, F);
  return("NA");
}

qc.boxplot <- function(dat, imgNm, dpi=72, format="png", interactive=F){
  dpi <- as.numeric(dpi)
  require('ggplot2')
  require('lattice');
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="");
  imgNm <- paste0(fileNm, format, sep="");  subgene <- 10000;

  if(class(dat)[1] == "data.frame"){
    dat <- as.matrix(dat);
  }

  if (nrow(dat)>subgene) {
    set.seed(28051968);
    sg  <- sample(nrow(dat), subgene);
    Mss <- dat[sg,,drop=FALSE];
  } else {
    Mss <- dat;
  }
  
  subsmpl=100;
  if (ncol(Mss)>subsmpl) {
    set.seed(28051968);
    ss  <- sample(ncol(Mss), subsmpl)
    Mss <- Mss[,ss,drop=FALSE]
  } else {
    Mss <- Mss
  }
  # OPTIMIZED: Direct data.frame construction to eliminate 3 intermediate copies
  df <- data.frame(
    values = as.numeric(Mss),
    sample_id = factor(rep(seq_len(ncol(Mss)), each = nrow(Mss)))
  )

  # OPTIMIZED: Single quantile call instead of duplicate calculations
  q_vals <- quantile(df$values, probs = c(0.01, 0.99), na.rm = TRUE)
  xlower <- unname(q_vals[1])
  xupper <- unname(q_vals[2])
  height <- length(unique(df$sample_id)) *20;
  if(height<450){
    height <- 450
  }
  bp <- ggplot(df, aes(sample_id, values)) +
    ylab("Values") + 
    xlab("Samples") + 
    scale_x_discrete(labels=colnames(dat)) + 
    ylim(xlower, xupper) + 
    stat_boxplot(geom = "errorbar", color="black") + 
    geom_boxplot(outlier.size=0.5, outlier.alpha=0.4) +
    theme_bw()
  bp <- bp + coord_flip();


    str <- "NA"

  if(interactive){
    library(plotly);
        m <- list(
                l = 50,
                r = 50,
                b = 20,
                t = 20,
                pad = 0.5
            )
    if(length(dataSet$meta.info) == 2){
    w=1000;
    }else{
    w=800;
    }
    ggp_build <- layout(ggplotly(bp), autosize = FALSE, width = w, height = 600, margin = m)
    return(ggp_build);
  }else{
  imgSet <- readSet(imgSet, "imgSet");
  if(grepl("norm", imgNm)){
    imgSet$qc_norm_boxplot <- imgNm;
  }else{
    imgSet$qc_boxplot <- imgNm;
  }
    saveSet(imgSet, "imgSet");
  if(dpi == 72){
  dpi <- dpi *1.34
  }
  Cairo(file=imgNm, width=600*dpi/72, height=height*dpi/72, unit="px",dpi=dpi, type=format, bg="white");
  print(bp);
  dev.off();
  return("NA")
  }
}


PlotDataDensity <- function(fileName, imgNm, dpi,format){
  dataSet <- readDataset(fileName);
  res <- qc.density(dataSet, imgNm, dpi, format, FALSE);
  return(res);
}

qc.density<- function(dataSet, imgNm="abc", dpi=72, format, interactive){
  require("ggplot2")
  dat <- dataSet$data.norm
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="");
  imgNm <- paste0(fileNm, format, sep="");
  dpi <- as.numeric(dpi)

  ######check.names=F is important for names not following conventions (i.e. with -, starts with number)
  df <- data.frame(dataSet$data.norm, stringsAsFactors = FALSE, check.names = FALSE)
  df <- stack(df)
  sampleNms <-colnames(dataSet$data.norm)
  if(length(dataSet$meta.info) == 2){

    # OPTIMIZED: Single merge instead of sequential merges to eliminate intermediate copies
    Factor1 <- dataSet$meta.info[,1]
    Factor2 <- dataSet$meta.info[,2]
    factorNm1 <- colnames(dataSet$meta.info)[1]
    factorNm2 <- colnames(dataSet$meta.info)[2]

    # Build combined metadata once
    conv <- data.frame(
      ind = sampleNms,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    conv[[factorNm1]] <- Factor1
    conv[[factorNm2]] <- Factor2

    # Single merge operation
    df1 <- merge(df, conv, by="ind")

    # Create separate plots for each factor
    d1 <- ggplot(df1, aes(x=values)) +
      geom_line(aes(color=.data[[factorNm1]], group=ind), stat="density", alpha=0.6) +
      labs(color = factorNm1) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
      )

    d2 <- ggplot(df1, aes(x=values)) +
      geom_line(aes(color=.data[[factorNm2]], group=ind), stat="density", alpha=0.6) +
      labs(color = factorNm2) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
      )

    if (requireNamespace("patchwork", quietly = TRUE)) {
      g <- patchwork::wrap_plots(d1, d2, ncol = 2)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
      g <- cowplot::plot_grid(d1, d2, ncol = 2)
    } else {
      stop("Missing dependency: install 'patchwork' or 'cowplot' for side-by-side density plots with separate legends.")
    }

    width <- 12
    height <- 5
  }else{
    Conditions <- dataSet$meta.info[,1];
    conv <- data.frame(ind=sampleNms, Conditions=Conditions, stringsAsFactors = FALSE, check.names = FALSE)
    df1 <- merge(df, conv, by="ind")
    
    g = ggplot(df1, aes(x=values)) + 
      geom_line(aes(color=Conditions, group=ind), stat="density", alpha=0.6) +
      theme_bw()
    
    width <- 8
    height <- 6
  }
  
  if(interactive){
    library(plotly);
    m <- list(
      l = 50,
      r = 50,
      b = 20,
      t = 20,
      pad = 0.5
    )
    if(length(dataSet$meta.info) == 2){
      w=1000;
    }else{
      w=800;
    }
    ggp_build <- layout(ggplotly(g), autosize = FALSE, width = w, height = 600, margin = m)
    return(ggp_build);
  }else{
    imgSet <- readSet(imgSet, "imgSet");
    imgSet$qc.density_norm <- imgNm;
    saveSet(imgSet, "imgSet");
  if(dpi == 72){
  dpi <- dpi *1.34
  }
    Cairo(file=imgNm, width=width, height=height, type=format, bg="white", dpi=dpi, unit="in");
    print(g)
    dev.off();
    return("NA")
  }
}


PlotLibSizeView<-function(fileName, imgNm,dpi=72, format="png"){
  library("ggrepel");
  require("ggplot2");

  dataSet <- readDataset(fileName);
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="");
  imgNm <- paste0(fileNm, format, sep="");
  dpi <- as.numeric(dpi);

  data.anot <- .get.annotated.data();
  data_bef<-data.matrix(data.anot);
  
  smpl.sums <- colSums(data_bef);
  
  names(smpl.sums) <- sampleNms <- colnames(data_bef);
  df <- data.frame(count=smpl.sums,ind=colnames(data_bef))
  
  if(length(dataSet$meta.info) == 2){
    Factor1 <- as.vector(dataSet$meta.info[,1])
    factor1Nm <- colnames(dataSet$meta.info)[1]
    conv <- data.frame(ind=sampleNms, Factor1=Factor1)
    colnames(conv) <- c("ind", factor1Nm)
    df1 <- merge(df, conv, by="ind")
    Factor2 <- as.vector(dataSet$meta.info[,2])
    factor2Nm <- colnames(dataSet$meta.info)[2]
    conv <- data.frame(ind=sampleNms, Factor2=Factor2)
    colnames(conv) <- c("ind", factor2Nm)
    df1 <- merge(df1, conv, by="ind")

    # Create separate plots for each factor
    if(length(df1$ind)>20){
      l1 <- ggplot(df1, aes(x = .data[[factor1Nm]], y = count, fill=.data[[factor1Nm]], label=ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',position = position_dodge(), dotsize=0.7) +
        geom_text() +
        ylab("Sum") +
        xlab(factor1Nm) +
        labs(fill = factor1Nm) +
        theme_bw() +
        theme(
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.position = "bottom"
        )

      plotData <- ggplot_build(l1)
      l1$layers[[2]] = NULL;

      l2 <- ggplot(df1, aes(x = .data[[factor2Nm]], y = count, fill=.data[[factor2Nm]], label=ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',position = position_dodge(), dotsize=0.7) +
        geom_text() +
        ylab("Sum") +
        xlab(factor2Nm) +
        labs(fill = factor2Nm) +
        theme_bw() +
        theme(
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.position = "bottom"
        )

      plotData <- ggplot_build(l2)
      l2$layers[[2]] = NULL;
    }else{
      l1 <- ggplot(df1, aes(x = .data[[factor1Nm]], y = count, fill=.data[[factor1Nm]], label=ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), dotsize=0.7) +
        geom_text_repel(force=5) +
        ylab("Sum") +
        xlab(factor1Nm) +
        labs(fill = factor1Nm) +
        theme_bw() +
        theme(
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.position = "bottom"
        )

      l2 <- ggplot(df1, aes(x = .data[[factor2Nm]], y = count, fill=.data[[factor2Nm]], label=ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), dotsize=0.7) +
        geom_text_repel(force=5) +
        ylab("Sum") +
        xlab(factor2Nm) +
        labs(fill = factor2Nm) +
        theme_bw() +
        theme(
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.position = "bottom"
        )
    }

    if (requireNamespace("patchwork", quietly = TRUE)) {
      g <- patchwork::wrap_plots(l1, l2, ncol = 2)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
      g <- cowplot::plot_grid(l1, l2, ncol = 2)
    } else {
      stop("Missing dependency: install 'patchwork' or 'cowplot' for side-by-side library size plots with separate legends.")
    }

    width <- 12
    height <- 6
    
  }else{
    Conditions= as.character(dataSet$meta.info[,1]);
    conv <- data.frame(ind=sampleNms, Conditions=Conditions)
    df1 <- merge(df, conv, by="ind")
    if(length(df1$ind)>20){
      g <- ggplot(df1, aes(x = Conditions, y = count, fill=Conditions, label= ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), dotsize=0.7) + 
        xlab("Sum") + 
        geom_text() +
        theme_bw()

      plotData <- ggplot_build(g)
      g$layers[[2]] = NULL;
    }else{

      g <- ggplot(df1, aes(x = Conditions, y = count, label=ind, fill=Conditions))+
        geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), dotsize=0.7) + 
        geom_text_repel(force=5) + 
        xlab("Sum") +
        theme_bw()

      plotData <- ggplot_build(g)
      
    }
    width <- 8
    height <- 6
    
  }
  if(dpi == 72){
  dpi <- dpi *1.34
  }
  Cairo(file=imgNm, width=width, height=height, type=format, bg="white", unit="in", dpi=dpi);
  print(g);
  dev.off();
  str <- "NA"

  imgSet <- readSet(imgSet, "imgSet");
  imgSet$libsize <- imgNm;
    saveSet(imgSet, "imgSet");

  return(str);
}

PlotDataMeanStd <- function(fileName, imgName, dpi,format){
  dataSet <- readDataset(fileName);
  if(grepl("_norm", imgName)){
    res <- qc.meanstd(dataSet$data.norm, imgName, dpi, format);
  }else{
    data.anot <- .get.annotated.data();
    res <- qc.meanstd(data.anot, imgName, dpi, format);
  }
  return(res);
}

qc.meanstd <- function(dat, imgNm,dpi=72, format="png"){
  dpi <- as.numeric(dpi)
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="");
  imgNm <- paste0(fileNm, format, sep="");
  #print(format)
  if(dpi == 72){
  dpi <- dpi *1.34
  }
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
  plot <- meanSdPlot(dat, ranks=FALSE) 
  dev.off();
  str <- "NA"
  
  imgSet <- readSet(imgSet, "imgSet");
  imgSet$qc.meanstd <- imgNm;
    saveSet(imgSet, "imgSet");

  return(str);
}

meanSdPlot <- function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
                       ylab = "sd", pch, plot = TRUE, bins = 50, ...) {
  
  stopifnot(is.logical(ranks), length(ranks) == 1, !is.na(ranks))
  
  n <- nrow(x)
  if (n == 0L) {
    warning("In 'meanSdPlot': input matrix 'x' has 0 rows. There is nothing to be done.")
    return()
  }
  if (!missing(pch)) {
    warning("In 'meanSdPlot': 'pch' is ignored.")
  }
  
  px   = rowMeans(x, na.rm = TRUE)
  py   = sqrt(rowV(x, mean = px, na.rm = TRUE))
  rpx  = rank(px, na.last = FALSE, ties.method = "random")
  
  ## run median with centers at dm, 2*dm, 3*dm,... and width 2*dm
  dm        = 0.025
  midpoints = seq(dm, 1-dm, by = dm)
  within    = function(x, x1, x2) { (x >= x1) & (x <= x2) }
  mediwind  = function(mp) median(py[within(rpx/n, mp - 2*dm, mp + 2*dm)], na.rm = TRUE)
  rq.sds    = sapply(midpoints, mediwind)
  
  res = if(ranks) {
    list(rank = midpoints*n, sd = rq.sds, px = rpx, py = py)
  } else {
    list(quantile = quantile(px, probs = midpoints, na.rm = TRUE), sd = rq.sds, px = px, py = py)
  }
  
  fmt = function() function(x) format(round(x, 0), nsmall = 0L, scientific = FALSE)
  
  res$gg = ggplot(data.frame(px = res$px, py = res$py),
                  aes_string(x = "px", y = "py")) + 
    xlab(xlab) + ylab(ylab) +
    geom_hex(bins = bins, ...) +
    scale_fill_gradient(name = "count", trans = "log", labels = fmt()) + 
    geom_line(aes_string(x = "x", y = "y"),
              data = data.frame(x = res[[1]], y = res$sd), color = "red") +
    theme_bw();
  
  if (plot) print(res$gg)
  
  return(invisible(res))
}

rowV = function(x, mean, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<1]  = NA
  if(missing(mean))
    mean=rowMeans(x, ...)
  return(rowSums(sqr(x-mean), ...)/(n-1))
}



PlotDataPCA <- function(fileName, imgName, dpi, format){
  dataSet <- readDataset(fileName);
  if(grepl("_norm", imgName)){
    # qc.pcaplot modifies dataSet$meta.info by filtering, so we capture the return value
    dataSet <- qc.pcaplot(dataSet, dataSet$data.norm, imgName, dpi, format, F);
        if (paramSet$oneDataAnalType == "dose") {
    qc.pcaplot.outliers.json(dataSet, dataSet$data.norm, imgName);
}else{
    qc.pcaplot.json(dataSet, dataSet$data.norm, imgName);

}

  }else{
    data.anot <- .get.annotated.data();
    # qc.pcaplot modifies dataSet$meta.info by filtering, so we capture the return value
    dataSet <- qc.pcaplot(dataSet, data.anot, imgName, dpi, format, F);
        if (paramSet$oneDataAnalType == "dose") {

   qc.pcaplot.outliers.json(dataSet, data.anot, imgName);
}else{
    qc.pcaplot.json(dataSet, data.anot, imgName);

}

  }
  return("NA");
}


qc.pcaplot <- function(dataSet, x, imgNm, dpi=72, format="png", interactive=FALSE) {
  dpi <- as.numeric(dpi)
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="")
  imgNm <- paste0(fileNm, format, sep="")

  require('lattice')
  require('ggplot2')
  require('reshape')
  require('see')
  require('ggrepel')
  
  get_divergent_palette <- function(n) {
    n <- max(1, n)
    grDevices::hcl.colors(n, "Dark 3")
  }

  # Force garbage collection before PCA to free memory
  gc(verbose = FALSE, full = TRUE)

  # Use memory-efficient PCA for large matrices to prevent malloc errors
  # First remove NA rows, then transpose (samples become rows)
  x_clean <- na.omit(x)
  data_matrix <- t(x_clean)

  # CRITICAL: Preserve sample names as rownames
  # In the original data, samples are columns, so colnames become rownames after transpose
  sample_names <- colnames(x_clean)
  rownames(data_matrix) <- sample_names
  matrix_size <- nrow(data_matrix) * ncol(data_matrix)

  # For matrices larger than 100k elements, use irlba (truncated SVD)
  if (matrix_size > 100000) {
    # Try to load irlba package
    if (!requireNamespace("irlba", quietly = TRUE)) {
      install.packages("irlba", repos = "https://cloud.r-project.org/")
    }

    require('irlba')

    # Use truncated SVD - compute only first 10 PCs (sufficient for visualization)
    n_components <- min(10, ncol(data_matrix), nrow(data_matrix) - 1)
    pca <- prcomp_irlba(data_matrix, n = n_components, center = TRUE, scale. = FALSE)

    # CRITICAL: Restore sample names to PCA results (irlba doesn't preserve them)
    rownames(pca$x) <- sample_names

    # Calculate variance importance for irlba (not included by default)
    variance_explained <- pca$sdev^2
    total_variance <- sum(variance_explained)
    prop_var <- variance_explained / total_variance
    cum_var <- cumsum(prop_var)

    # Create importance matrix compatible with standard prcomp
    imp.pca <- rbind(
      "Standard deviation" = pca$sdev,
      "Proportion of Variance" = prop_var,
      "Cumulative Proportion" = cum_var
    )
    colnames(imp.pca) <- paste0("PC", 1:ncol(imp.pca))

  } else {
    pca <- prcomp(data_matrix)

    # Ensure sample names are preserved (prcomp usually does this, but make sure)
    if (is.null(rownames(pca$x))) {
      rownames(pca$x) <- sample_names
    }

    imp.pca <- summary(pca)$importance
  }

  # Clean up temporary matrix
  rm(data_matrix)
  gc(verbose = FALSE)
  xlabel <- paste0("PC1"," (", 100 * round(imp.pca[2,][1], 3), "%)")
  ylabel <- paste0("PC2"," (", 100 * round(imp.pca[2,][2], 3), "%)")
  pca.res <- as.data.frame(pca$x)
  pca.res <- pca.res[, c(1, 2)]

  # Store original rownames before any filtering
  original_rownames <- rownames(pca.res)
  if ("newcolumn" %in% colnames(dataSet$meta.info)) {
    dataSet$meta.info <- data.frame(dataSet$meta.info[, -which(colnames(dataSet$meta.info) == "newcolumn")])
  }

  # First align with metadata BEFORE filtering non-finite values
  # Find common rows between PCA and metadata
  common_rows <- intersect(rownames(pca.res), rownames(dataSet$meta.info))
  if (length(common_rows) == 0) {
    stop("[qc.pcaplot] ERROR: No common samples between PCA results and metadata!")
  }

  # Align both datasets to common samples
  pca.res <- pca.res[common_rows, , drop = FALSE]
  dataSet$meta.info <- dataSet$meta.info[common_rows, , drop = FALSE]

  # NOW check for non-finite values (after alignment)
  if (any(!is.finite(pca.res$PC1)) || any(!is.finite(pca.res$PC2))) {
    # Filter both PCA and metadata together to keep them aligned
    valid_rows <- is.finite(pca.res$PC1) & is.finite(pca.res$PC2)
    pca.res <- pca.res[valid_rows, , drop = FALSE]
    dataSet$meta.info <- dataSet$meta.info[valid_rows, , drop = FALSE]
  }

  # Increase xlim and ylim for text label
  xlim <- GetExtendRange(pca.res$PC1)
  ylim <- GetExtendRange(pca.res$PC2)
  
  if (length(dataSet$meta.info) == 2) {
    # OPTIMIZED: Get column names once
    meta_colnames <- colnames(dataSet$meta.info)
    Factor1 <- as.vector(dataSet$meta.info[, 1])
    factorNm1 <- meta_colnames[1]
    pca.res[, factorNm1] <- Factor1
    Factor2 <- as.vector(dataSet$meta.info[, 2])
    factorNm2 <- meta_colnames[2]
    pca.res[, factorNm2] <- Factor2
    pca.rest <- reshape::melt(pca.res, measure.vars = c(factorNm1, factorNm2))
    colnames(pca.rest)[4] <- "Conditions"
    pca.rest$names <- rep(rownames(pca.res), times = 2)

    # Remove rows with NA or non-finite values before aggregation
    pca.rest.clean <- pca.rest[is.finite(pca.rest$PC1) & is.finite(pca.rest$PC2) & !is.na(pca.rest$Conditions), ]

    if (nrow(pca.rest.clean) == 0) {
      # Use original data without centroid/outlier detection
      pca.rest.clean <- pca.rest
      pca.rest.clean$outlier <- FALSE
    } else {
      # Calculate group centroids
      centroids <- aggregate(. ~ Conditions, data = pca.rest.clean[, c("PC1", "PC2", "Conditions")], mean)

      # Save names column before merge (merge will reorder/modify data)
      names_col <- pca.rest.clean$names

      # Merge centroids back to the pca.rest dataframe
      pca.rest.clean <- merge(pca.rest.clean, centroids, by = "Conditions", suffixes = c("", "_centroid"))

      # CRITICAL: Restore rownames using the names column after merge
      # The merge operation changes row order and drops rownames
      rownames(pca.rest.clean) <- paste0(pca.rest.clean$names, "_", seq_len(nrow(pca.rest.clean)))

      # Calculate the distance to the centroid
      pca.rest.clean$distance <- sqrt((pca.rest.clean$PC1 - pca.rest.clean$PC1_centroid)^2 + (pca.rest.clean$PC2 - pca.rest.clean$PC2_centroid)^2)
      # Identify outliers based on variance threshold (20% here)
      threshold <- 0.2 * mean(pca.rest.clean$distance, na.rm = TRUE)
      pca.rest.clean$outlier <- pca.rest.clean$distance > threshold
    }

    # Use cleaned data for plotting
    pca.rest <- pca.rest.clean
    
    p1_df <- pca.rest[pca.rest$variable == factorNm1, , drop = FALSE]
    p2_df <- pca.rest[pca.rest$variable == factorNm2, , drop = FALSE]
    pal1 <- get_divergent_palette(length(unique(p1_df$Conditions)))
    pal2 <- get_divergent_palette(length(unique(p2_df$Conditions)))
    
    p1 <- ggplot(p1_df, aes(x = PC1, y = PC2, color = Conditions, label = names)) +
      geom_point(size = 3, alpha = 0.6) +
      xlim(xlim) +
      ylim(ylim) +
      xlab(xlabel) +
      ylab(ylabel) +
      theme_bw() +
      scale_color_manual(values = pal1) +
      guides(color = guide_legend(title = factorNm1)) +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
      )

    p2 <- ggplot(p2_df, aes(x = PC1, y = PC2, color = Conditions, label = names)) +
      geom_point(size = 3, alpha = 0.6) +
      xlim(xlim) +
      ylim(ylim) +
      xlab(xlabel) +
      ylab(ylabel) +
      theme_bw() +
      scale_color_manual(values = pal2) +
      guides(color = guide_legend(title = factorNm2)) +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
      )
    
    if (requireNamespace("patchwork", quietly = TRUE)) {
      pcafig <- patchwork::wrap_plots(p1, p2, ncol = 2)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
      pcafig <- cowplot::plot_grid(p1, p2, ncol = 2)
    } else {
      stop("Missing dependency: install 'patchwork' or 'cowplot' for side-by-side PCA plots with separate legends.")
    }
    width <- 12
    height <- 6
  } else {
    Factor <- dataSet$meta.info[, 1]
    pca.res$Conditions <- Factor
    pca.res$names <- rownames(pca.res)

    # Save original rownames before merge (merge drops rownames!)
    original_sample_names <- rownames(pca.res)

    # Calculate group centroids
    centroids <- aggregate(. ~ Conditions, data = pca.res[, c("PC1", "PC2", "Conditions")], mean)
    # Merge centroids back to the pca.res dataframe
    pca.res <- merge(pca.res, centroids, by = "Conditions", suffixes = c("", "_centroid"))

    # CRITICAL: Restore original rownames after merge (they were lost!)
    rownames(pca.res) <- pca.res$names

    # Calculate the distance to the centroid
    pca.res$distance <- sqrt((pca.res$PC1 - pca.res$PC1_centroid)^2 + (pca.res$PC2 - pca.res$PC2_centroid)^2)
    # Identify outliers based on variance threshold (20% here)
    threshold <- 0.2 * mean(pca.res$distance, na.rm = TRUE)
    pca.res$outlier <- pca.res$distance > threshold
    
      pcafig <- ggplot(pca.res, aes(x = PC1, y = PC2, color = Conditions)) +
        geom_point(size = 3, alpha = 0.6) +
        xlim(xlim) +
        ylim(ylim) +
        xlab(xlabel) +
        ylab(ylabel) +
        theme_bw() +
        theme(
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)
        )
    
    width <- 8
    height <- 6
    pal <- get_divergent_palette(length(unique(pca.res$Conditions)))
    pcafig <- pcafig + scale_fill_manual(values = pal) + scale_color_manual(values = pal)
  }
  

  # Print the outliers
  outliers <- pca.res[pca.res$outlier, "names"]
  if(!is.null(outliers)){
    paramSet <- readSet(paramSet, "paramSet")
    paramSet$pca.outliers <- outliers;
  }else{
    paramSet$pca.outliers <- c("NA");
  }

 # permanova_results ="";

  # Run PERMANOVA with error handling to prevent crashes
  #trash <- capture.output({
  #  permanova_results <- tryCatch({
  #    ComputePERMANOVA(pca.res$PC1,
  #                     pca.res$PC2,
  #                     dataSet$meta.info[,1],
  #                     999)
  #  }, error = function(e) {
  #    print(paste("[qc.pcaplot] WARNING: PERMANOVA failed:", e$message));
  #    print("[qc.pcaplot] Continuing without PERMANOVA statistics");
  #    # Return empty result
  #    list(
  #      stat.info = "PERMANOVA test skipped due to error",
  #      stat.info.vec = c(`F-value` = NA, `R-squared` = NA, `p-value` = NA),
  #      pair.res = NULL
  #    )
  #  })
  #}, file = NULL)  # discard output

  analSet <- readSet(analSet, "analSet");

  # pca is a large object. only save those required for json
  # never use more than top 3. Update this if you require more PCs

  # Use the importance matrix we already computed (works for both prcomp and irlba)
  imp     <- imp.pca[2, 1:2]
  xlabel  <- sprintf("PC1 (%.1f%%)", 100 * imp[1])
  ylabel  <- sprintf("PC2 (%.1f%%)", 100 * imp[2])

  # Select top 3 PCs (or fewer if less than 3 PCs were computed)
  # CRITICAL: Save the FILTERED pca results (matching the filtered metadata)
  # pca.res already contains filtered PC1 and PC2, we need to filter pca$x the same way
  n_pcs_to_save <- min(3, ncol(pca$x))

  # Get the samples that remain after filtering (these are in pca.res)
  filtered_samples <- rownames(pca.res)

  # Subset pca$x to only include the filtered samples
  # Safely subset - check both dimensions
  actual_n_pcs <- min(n_pcs_to_save, ncol(pca$x), length(filtered_samples) - 1)  # PCA can have at most n-1 components

  filtered_pca_x <- pca$x[filtered_samples, 1:actual_n_pcs, drop = FALSE]

  my.pca <- list(
        x = filtered_pca_x,  # Use filtered PCA results, not original
        xlabel = xlabel,
        ylabel = ylabel
    );

  analSet$pca <- my.pca;
  #analSet$permanova.res <-permanova_results;
  saveSet(analSet, "analSet");
  saveSet(paramSet, "paramSet");
  
  if (interactive) {
    library(plotly)
    m <- list(l=50, r=50, b=20, t=20, pad=0.5)
    w <- if (length(dataSet$meta.info)==2) 1000 else 800
    ggp_build <- layout( ggplotly(pcafig), autosize=FALSE, width=w, height=600, margin=m )

    return(dataSet)  # Return modified dataSet even in interactive mode
  } else {
  if(dpi == 72){
  dpi <- dpi *1.34
  }
  imgSet <- readSet(imgSet, "imgSet")
  if (grepl("norm", imgNm)) {
    imgSet$qc_norm_pca <- imgNm
  } else {
    imgSet$qc_pca <- imgNm
  }
  saveSet(imgSet, "imgSet")
    Cairo(file = imgNm, width=width, height=height, type=format, bg="white", unit="in", dpi=dpi)
    if (length(dataSet$meta.info) == 2) {
      print(pcafig)
    } else {
      print(pcafig)
    }
    dev.off()
    return(dataSet)  # Return modified dataSet with filtered metadata
  }
}

GetPcaOutliers <- function(){
    paramSet <- readSet(paramSet, "paramSet")
    return(paramSet$pca.outliers);
}

PlotDataNcov5 <- function(fileName, imgName, dpi, format){
  dataSet <- readDataset(fileName)
  if (is.null(dataSet$summary_df)) {
    stop("summary_df not found in dataSet. Please run SummarizeQC first.")
  }

  ncov5_df <- dataSet$summary_df[, c("Sample", "HighCoverageGeneCount")]
  
  ## Compute outlier limits
  Q1  <- quantile(ncov5_df$HighCoverageGeneCount, 0.25)
  Q3  <- quantile(ncov5_df$HighCoverageGeneCount, 0.75)
  IQRv <- IQR(ncov5_df$HighCoverageGeneCount)
  lower <- Q1 - 3 * IQRv
  upper <- Q3 + 3 * IQRv

  ncov5_df$Status <- ifelse(ncov5_df$HighCoverageGeneCount < lower |
                            ncov5_df$HighCoverageGeneCount > upper,
                            "Outlier", "Normal")

  qc.ncov5.plot(ncov5_df, imgName, lower, upper, dpi, format);
  qc.ncov5plot.json(ncov5_df, imgName, lower, upper);
  return("NA")
}

qc.ncov5.plot <- function(ncov5_df,
                          imgNm = "NCov5_plot",
                          lower,
                          upper,
                          dpi = 72,
                          format = "png",
                          interactive = FALSE) {
  require(ggplot2)
  require(ggrepel)
  require(Cairo)
  
  dpi <- as.numeric(dpi)
  if (dpi <= 0) stop("DPI must be a positive number.")

  g <- ggplot(ncov5_df, aes(x = "", y = HighCoverageGeneCount)) +
    geom_boxplot(outlier.shape = NA, fill = "grey80") +
    geom_jitter(aes(color = Status), width = 0.25, height = 0) +
    geom_hline(yintercept = c(lower, upper),
               linetype = "dashed", color = "blue") +
    geom_text_repel(data = subset(ncov5_df, Status == "Outlier"),
                    aes(label = Sample), nudge_x = 0.35, size = 3) +
    scale_color_manual(values = c(Normal = "grey40", Outlier = "red"),
                       name = "Sample status") +
    theme_minimal(base_size = 11) +
    labs(x = NULL,
         y = "Genes with > 5 uniquely mapped reads") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
  
  width  <- 8
  height <- 6
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  imgNm  <- paste0(fileNm, format)

  if (interactive) {
    require(plotly)
    m <- list(l = 50, r = 50, b = 20, t = 20, pad = 0.5)
    return(layout(plotly::ggplotly(g),
                  autosize = FALSE, width = 1000, height = 600, margin = m))
  } else {
    if (dpi == 72) dpi <- dpi * 1.34
    Cairo(file = imgNm, width = width, height = height,
          type = format, bg = "white", dpi = dpi, unit = "in")
    print(g)
    dev.off()
    return("NA")
  }
}


PlotDataNsig <- function(fileName, imgName, dpi, format){
  dataSet <- readDataset(fileName)
  if (is.null(dataSet$summary_df)) {
    stop("summary_df not found in dataSet. Please run SummarizeQC first.")
  }

  nsig_df <- dataSet$summary_df[, c("Sample", "NSig80")]


  ## identify outliers (± 3×IQR)
  Q1  <- quantile(nsig_df$NSig80, 0.25)
  Q3  <- quantile(nsig_df$NSig80, 0.75)
  IQRv <- IQR(nsig_df$NSig80)
  lower <- Q1 - 3 * IQRv
  upper <- Q3 + 3 * IQRv

  nsig_df$outlier <- ifelse(nsig_df$NSig80 < lower | nsig_df$NSig80 > upper,
                            "Outlier", "Normal")

  qc.nsig.plot(nsig_df, imgName, lower, upper, dpi, format)
  qc.nsigplot.json(nsig_df, imgName, lower, upper); 
  return("NA")
}

qc.nsig.plot <- function(nsig_df,
                         imgNm = "NSig80_plot",
                         lower,
                         upper,
                         dpi = 72,
                         format = "png",
                         interactive = FALSE) {
  require("ggplot2")
  require("Cairo")
  require("ggrepel")

  dpi <- as.numeric(dpi)
  if (dpi <= 0) stop("DPI must be a positive number.")

  g <- ggplot(nsig_df, aes(x = "", y = NSig80)) +
    geom_boxplot(outlier.shape = NA, fill = "grey80") +
    geom_jitter(aes(color = outlier), width = 0.25, height = 0) +
    scale_color_manual(values = c("Normal" = "grey40", "Outlier" = "red")) +
    geom_text_repel(data = subset(nsig_df, outlier == "Outlier"),
                    aes(label = Sample), nudge_x = 0.35, size = 3) +
    geom_hline(yintercept = c(lower, upper), linetype = "dashed",
               color = "blue") +
    theme_minimal(base_size = 11) +
    labs(x = NULL,
         y = "NSig80 (genes reaching 80 % of signal)",
         color = "Sample Status") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  width  <- 8
  height <- 6
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  imgNm  <- paste0(fileNm, format)

  if (interactive) {
    require("plotly")
    m <- list(l = 50, r = 50, b = 20, t = 20, pad = 0.5)
    return(layout(plotly::ggplotly(g),
                  autosize = FALSE, width = 1000, height = 600, margin = m))
  } else {
    if (dpi == 72) dpi <- dpi * 1.34
    Cairo(file = imgNm, width = width, height = height,
          type = format, bg = "white", dpi = dpi, unit = "in")
    print(g)
    dev.off()
    return("NA")
  }
}

PlotDataDendrogram <- function(fileName, imgName, threshold, dpi, format){
  dataSet <- readDataset(fileName)

  if (is.null(dataSet$summary_df)) {
    stop("summary_df not found in dataSet. Please run SummarizeQC first.")
  }
  dendro_df <- dataSet$summary_df[, c("Sample", "Dendrogram_Distance")]
  dendro_df$Status <- ifelse(dendro_df$Dendrogram_Distance > threshold, "Outlier", "Normal")

  ## Decide label set
  out_idx <- which(dendro_df$Status == "Outlier")
  label_idx <- if (length(out_idx) <= 20) {
    out_idx
  } else {
    out_idx[order(dendro_df$Dendrogram_Distance[out_idx], decreasing = TRUE)[1:20]]
  }
  dendro_df$LabelMe <- FALSE
  dendro_df$LabelMe[label_idx] <- TRUE


  qc.dendrogram.plot(dendro_df, threshold, imgName, dpi, format)
  qc.dendrogram.json(dendro_df, imgName);
  return("NA")
}

qc.dendrogram.plot <- function(dendro_df,
                               threshold = 0.1,
                               imgNm = "Dendrogram_plot",
                               dpi = 72,
                               format = "png",
                               interactive = FALSE) {
  require(ggplot2)
  require(ggrepel)
  require(Cairo)

  dpi <- as.numeric(dpi)
  if (dpi <= 0) stop("DPI must be positive.")

  set.seed(1)  # For reproducible jitter
  dendro_df$xj <- jitter(rep(1, nrow(dendro_df)), amount = 0.25)

  g <- ggplot(dendro_df, aes(x = xj, y = Dendrogram_Distance)) +
    geom_boxplot(aes(x = 1), outlier.shape = NA,
                 width = 0.4, fill = "grey80") +
    geom_point(aes(color = Status), size = 2.2) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "blue") +
    geom_text_repel(data = dendro_df[dendro_df$LabelMe, ],
                    aes(label = Sample),
                    max.overlaps = Inf,
                    box.padding = 0.35,
                    segment.size = 0.2,
                    size = 4.2) +
    scale_color_manual(values = c(Normal = "grey40", Outlier = "red"),
                       name = "Sample status") +
    theme_minimal(base_size = 12) +
    labs(x = NULL, y = "Max pair-wise distance (1 − Pearson ρ)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  outFile <- paste0(imgNm, "dpi", dpi, ".", format)

  if (interactive) {
    require(plotly)
    m <- list(l = 50, r = 50, b = 20, t = 20, pad = 0.5)
    return(layout(plotly::ggplotly(g),
                  autosize = FALSE, width = 1000, height = 600, margin = m))
  } else {
    if (dpi == 72) dpi <- dpi * 1.34
    Cairo(file = outFile, width = 8, height = 6,
          type = format, bg = "white", dpi = dpi, unit = "in")
    print(g)
    dev.off()
    return(outFile)
  }
}



GetSummaryTable <- function(dataName){
  dataSet <- readDataset(dataName)
  df <- dataSet$summary_df;
  df_rounded <- df;
  df_rounded[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)], signif, digits = 3)
  return(df_rounded);
}

calculate_gini <- function(x) {
  n <- length(x)
  sorted_x <- sort(x)
  index <- 1:n
  gini <- (2 * sum(index * sorted_x) / sum(sorted_x)) - (n + 1)
  return(gini / n)
}


ComputePERMANOVA <- function(pc1, pc2, cls, numPermutations = 999) {
  # Combine PC1 and PC2 scores into a matrix
  pc.mat <- cbind(pc1, pc2)
  
  # Calculate PERMANOVA significance
  res <- .calculateDistSig(pc.mat, cls)
  
  # Extract the main results
  resTab <- res[[1]][1, ]
  
  # Format and create the PERMANOVA summary statistics
  stat.info <- paste("[PERMANOVA] F-value: ", signif(resTab$F, 5),
                     "; R-squared: ", signif(resTab$R2, 5),
                     "; p-value (based on ", numPermutations, " permutations): ",
                     signif(resTab$Pr, 5), sep = "")
  
  # Create a named vector for the statistics
  stat.info.vec <- c(F_value = signif(resTab$F, 5), 
                     R_squared = signif(resTab$R2, 5), 
                     p_value = signif(resTab$Pr, 5))
  names(stat.info.vec) <- c("F-value", "R-squared", "p-value");

  # Extract pairwise PERMANOVA results if available
  pair.res <- res[[2]]
  
  # Return the results as a list
  list(
    stat.info = stat.info,
    stat.info.vec = stat.info.vec,
    pair.res = pair.res
  )
}

# use a PERMANOVA to partition the euclidean distance by groups based on current score plot:
.calculateDistSig <- function(pc.mat, grp){

    # Force garbage collection before computationally intensive PERMANOVA
    print("[calculateDistSig] Starting PERMANOVA test...");
    gc(verbose = FALSE, full = TRUE)

    # Use fewer permutations for large datasets to prevent memory issues
    n_samples <- nrow(pc.mat)
    n_perms <- if (n_samples > 100) 99 else 999  # Reduce permutations for large datasets
    print(paste("[calculateDistSig] Using", n_perms, "permutations for", n_samples, "samples"));

    # Calculate distance matrix
    data.dist <- dist(as.matrix(pc.mat), method = 'euclidean');

    # Wrap adonis2 in error handling to catch malloc errors
    res <- tryCatch({
      vegan::adonis2(formula = data.dist ~ grp, permutations = n_perms)
    }, error = function(e) {
      print(paste("[calculateDistSig] ERROR in adonis2:", e$message));
      print("[calculateDistSig] Returning dummy result due to error");
      # Return a dummy result structure
      dummy_df <- data.frame(
        Df = c(1, n_samples - 1),
        SumOfSqs = c(0, 0),
        R2 = c(0, 1),
        F = c(0),
        Pr = c(1)  # p-value = 1 (not significant)
      )
      rownames(dummy_df) <- c("grp", "Residual")
      return(dummy_df)
    })

    # pairwise for multi-grp
    if(length(levels(grp)) > 2){
      pair.res <- tryCatch({
        .permanova_pairwise(x = data.dist, grp, permutations = n_perms);
      }, error = function(e) {
        print(paste("[calculateDistSig] ERROR in pairwise PERMANOVA:", e$message));
        return(NULL)
      })

      if (!is.null(pair.res)) {
        rownames(pair.res) <- pair.res$pairs;
        pair.res$pairs <- NULL;
        pair.res <- signif(pair.res,5);
        fast.write.csv(pair.res, file="pca_pairwise_permanova.csv");
      }
    }else{
      pair.res <- NULL;
    }

    print("[calculateDistSig] PERMANOVA test complete");
    gc(verbose = FALSE)  # Clean up after PERMANOVA

    return(list(res, pair.res));
}

###adopted from ecole package https://rdrr.io/github/phytomosaic/ecole/
.permanova_pairwise <- function(x,
                                 grp,
                                 permutations = 999,
                                 method = 'bray',
                                 padj = 'fdr', ...) {
  f     <- grp
  if (!all(table(f) > 1)) warning('factor has singletons! perhaps lump them?')
  co    <- combn(unique(as.character(f)),2)
  nco   <- NCOL(co)
  out   <- data.frame(matrix(NA, nrow=nco, ncol=5))
  dimnames(out)[[2]] <- c('pairs', 'SumOfSqs', 'F.Model', 'R2', 'pval')
  if (!inherits(x, 'dist')) {
    D <- vegan::vegdist(x, method=method)
  } else {
    D <- x
  }
  #cat('Now performing', nco, 'pairwise comparisons. Percent progress:\n')
  for(j in 1:nco) {
    cat(round(j/nco*100,0),'...  ')
    ij  <- f %in% c(co[1,j],co[2,j])
    Dij <- as.dist(as.matrix(D)[ij,ij])
    fij <- data.frame(fij = f[ij])
    a   <- vegan::adonis2(Dij ~ fij, data=fij, permutations = permutations, ...);
    out[j,1] <- paste(co[1,j], 'vs', co[2,j])
    out[j,2] <- a$SumOfSqs[1]
    out[j,3] <- a$F[1]
    out[j,4] <- a$R2[1]
    out[j,5] <- a$`Pr(>F)`[1]
  }
  #cat('\n')
  out$p.adj <- p.adjust(out$pval, method=padj)
  out$SumOfSqs <-NULL
  #attr(out, 'p.adjust.method') <- padj
  #cat('\np-adjust method:', padj, '\n\n');
  return(out)
}

qc.pcaplot.json <- function(dataSet, x, imgNm) {
  jsonFile <- paste0(imgNm, ".json")

  # libs (only jsonlite is really needed now)
  suppressMessages({
    require(jsonlite)
  })

  # ---------- Load PCA & metadata ----------
  analSet <- readSet(analSet, "analSet")
  pca     <- analSet$pca
  xlabel  <- pca$xlabel
  ylabel  <- pca$ylabel

  pca.res <- as.data.frame(pca$x)[, 1:2, drop = FALSE]
  colnames(pca.res) <- c("PC1", "PC2")

  # Find common samples between PCA and metadata (matching PNG generation logic)
  common_samples <- intersect(rownames(pca.res), rownames(dataSet$meta.info))
  if (length(common_samples) == 0) {
    stop("[qc.pcaplot.json] ERROR: No common samples between PCA results and metadata!")
  }

  # Align both PCA and metadata to common samples (instead of using match which can create NAs)
  pca.res <- pca.res[common_samples, , drop = FALSE]
  meta.info.aligned <- dataSet$meta.info[common_samples, , drop = FALSE]

  # metadata1 → color (group)
  pca.res$group  <- as.character(meta.info.aligned[[1]])
  pca.res$sample <- rownames(pca.res)

  # ---------- Detect 2nd metadata for shapes ----------
  doShape <- FALSE
  shape.levels <- character(0)
  shape.map <- NULL

  if (ncol(meta.info.aligned) >= 2) {
    second <- meta.info.aligned[[2]]
    # treat non-numeric as discrete for shapes
    isDisc  <- !is.numeric(second)
    levs    <- unique(as.character(second))
    if (isDisc && length(levs) <= 8) {
      doShape <- TRUE
      pca.res$shape <- as.character(second)
      symbols <- c("circle","square","diamond",
                   "cross","x","triangle-up",
                   "triangle-down","star")
      shape.map    <- stats::setNames(symbols[seq_along(levs)], levs)
      shape.levels <- levs
    }
  }

  # ---------- Color mapping (dose-aware like qc.pcaplot) ----------
  paramSet    <- readSet(paramSet, "paramSet")
  unique_grps <- unique(pca.res$group)

  if (grepl("norm", imgNm) &&
      !is.null(paramSet$oneDataAnalType) &&
      paramSet$oneDataAnalType == "dose") {
    pal <- grDevices::colorRampPalette(c("#2196F3", "#DE690D"))(length(unique_grps))
  } else {
    okabe <- c("#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7")
    pal <- rep(okabe, length.out = length(unique_grps))
  }
  col.map <- stats::setNames(pal, unique_grps)

  # ---------- Build traces ----------
  traces <- list()

  if (doShape) {
    # one trace per (group × shape)
    combos <- unique(pca.res[, c("group","shape")])
    for (i in seq_len(nrow(combos))) {
      g  <- combos$group[i]
      sh <- combos$shape[i]

      df <- pca.res[pca.res$group == g & pca.res$shape == sh, , drop = FALSE]
      if (nrow(df) == 0) next

      mkr <- list(
        color = unname(col.map[g]),
        size  = 8,
        line  = list(color = "white", width = 0.5),
        symbol = unname(shape.map[[sh]]) # scalar symbol per trace
      )

      traces[[length(traces) + 1]] <- list(
        x            = df$PC1,
        y            = df$PC2,
        type         = "scatter",
        mode         = if (nrow(df) > 20) "markers" else "markers+text",
        name         = paste0(g, " • ", sh),
        legendgroup  = g,          # groups align in legend
        marker       = mkr,
        text         = if (nrow(df) <= 20) df$sample else NULL,
        hovertext    = df$sample,  # Always include sample names in hover
        hoverinfo    = "text",
        textposition = "top center"
      )
    }
  } else {
    # one trace per color group (no shapes)
    for (g in unique_grps) {
      df <- pca.res[pca.res$group == g, , drop = FALSE]
      if (nrow(df) == 0) next

      mkr <- list(
        color = unname(col.map[g]),
        size  = 8,
        line  = list(color = "white", width = 0.5)
      )

      traces[[length(traces) + 1]] <- list(
        x            = df$PC1,
        y            = df$PC2,
        type         = "scatter",
        mode         = if (nrow(df) > 20) "markers" else "markers+text",
        name         = g,
        legendgroup  = g,
        marker       = mkr,
        text         = if (nrow(df) <= 20) df$sample else NULL,
        hovertext    = df$sample,  # Always include sample names in hover
        hoverinfo    = "text",
        textposition = "top center"
      )
    }
  }

  # ---------- Layout ----------
  layout <- list(
    title = "",
    xaxis = list(title = xlabel, zeroline = FALSE),
    yaxis = list(title = ylabel, zeroline = FALSE),
    legend = list(
      orientation = "v",
      x           = 1.02,
      y           = 1,
      xanchor     = "left",
      yanchor     = "top"
    )
  )

  # ---------- Dump JSON ----------
  plot_data <- list(data = traces, layout = layout)
  json.obj  <- jsonlite::toJSON(plot_data, auto_unbox = TRUE, null = "null", digits = NA)
  writeLines(json.obj, jsonFile)

  return("NA")
}


PlotDataGini <- function(fileName, imgName, threshold, dpi, format){
  dataSet <- readDataset(fileName)
  if (is.null(dataSet$summary_df)) {
    stop("summary_df not found in dataSet. Please run SummarizeQC first.")
  }
  
  # Select Gini data
  gini_df <- dataSet$summary_df[, c("Sample", "Gini")]
  gini_df$Status <- ifelse(gini_df$Gini > threshold, "Outlier", "Normal")
  
  ## Plot
  qc.gini.plot(gini_df, imgName, threshold, dpi, format)
  qc.giniplot.json(gini_df, imgName);
  return("NA")
}

qc.gini.plot <- function(gini_df,
                         imgNm   = "Gini_plot",
                         threshold = 0.95,
                         dpi     = 72,
                         format  = "png",
                         interactive = FALSE) {
  require(ggplot2)
  require(ggrepel)
  require(Cairo)
  
  dpi <- as.numeric(dpi)
  if (dpi <= 0) stop("DPI must be a positive number.")
  
  g <- ggplot(gini_df, aes(x = "", y = Gini)) +
    geom_boxplot(outlier.shape = NA, fill = "grey80") +
    geom_jitter(aes(color = Status), width = 0.25, height = 0) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "blue") +
    geom_text_repel(data = subset(gini_df, Status == "Outlier"),
                    aes(label = Sample), nudge_x = 0.35, size = 3) +
    scale_color_manual(values = c(Normal = "grey40", Outlier = "red"),
                       name = "Sample status") +
    theme_minimal(base_size = 11) +
    labs(x = NULL, y = "Gini coefficient") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
  
  width  <- 8
  height <- 6
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  imgNm  <- paste0(fileNm, format)
  
  if (interactive) {
    require(plotly)
    m <- list(l = 50, r = 50, b = 20, t = 20, pad = 0.5)
    return(layout(plotly::ggplotly(g),
                  autosize = FALSE, width = 1000, height = 600, margin = m))
  } else {
    if (dpi == 72) dpi <- dpi * 1.34
    Cairo(file = imgNm, width = width, height = height,
          type = format, bg = "white", dpi = dpi, unit = "in")
    print(g)
    dev.off()
    return("NA")
  }
}

SummarizeQC <- function(fileName, imgNameBase, threshold = 0.1) {
  # save.image("summarize.RData");
  dataSet <- readDataset(fileName)

  summary_df <- data.frame(
    Sample = character(),
    HighCoverageGeneCount = numeric(),
    NSig80 = numeric(),
    Gini = numeric(),
    Dendrogram_Distance = numeric(),
    Outlier_HighCoverageGeneCount = numeric(),
    Outlier_NSig80 = numeric(),
    Outlier_Gini = numeric(),
    Outlier_Dendrogram = numeric(),
    stringsAsFactors = FALSE
  )

  if (grepl("norm", imgNameBase)) {
    data <- dataSet$data.norm
  } else {
    data <- .get.annotated.data()
  }

  ## --- basic per-sample metrics ---
  HighCoverageGeneCount <- colSums(data > 5, na.rm = TRUE)
  ncov5_df <- data.frame(
    Sample = names(HighCoverageGeneCount),
    HighCoverageGeneCount = as.numeric(HighCoverageGeneCount),
    stringsAsFactors = FALSE
  )

  NSig80 <- apply(data, 2, function(col) {
    col[is.na(col)] <- 0
    s <- sum(col)
    if (s <= 0) return(0)
    sum(cumsum(sort(col, decreasing = TRUE)) <= 0.8 * s)
  })
  nsig_df <- data.frame(
    Sample = names(NSig80),
    NSig80 = as.numeric(NSig80),
    stringsAsFactors = FALSE
  )

  gini_scores <- apply(data, 2, calculate_gini)
  gini_df <- data.frame(
    Sample = colnames(data),
    Gini = as.numeric(gini_scores),
    stringsAsFactors = FALSE
  )

  ## --- correlation distance (0..1) & within-group mean distance ---
  pearson_corr <- suppressWarnings(cor(data, method = "pearson", use = "pairwise.complete.obs"))
  # keep diagonal sane and replace NA correlations
  if (!is.null(pearson_corr)) {
    diag(pearson_corr) <- 1
    pearson_corr[is.na(pearson_corr)] <- 0
  }
  distance_matrix <- as.dist((1 - pearson_corr) / 2)  # 0..1
  dist_mat <- as.matrix(distance_matrix)

  group_info <- dataSet$meta.info[, 1]
  names(group_info) <- rownames(dataSet$meta.info)

  mean_distances <- sapply(colnames(data), function(sample) {
    sample_group <- group_info[sample]
    same_group_samples <- names(group_info)[group_info == sample_group]
    same_group_samples <- setdiff(same_group_samples, sample)
    # ensure they exist in the distance matrix
    same_group_samples <- intersect(same_group_samples, colnames(data))
    if (length(same_group_samples) == 0 || is.null(dist_mat)) return(NA_real_)
    mean(dist_mat[sample, same_group_samples], na.rm = TRUE)
  })

  dendrogram_df <- data.frame(
    Sample = names(mean_distances),
    Dendrogram_Distance = as.numeric(mean_distances),
    stringsAsFactors = FALSE
  )

  ## --- Merge first (align by Sample) ---
  summary_df <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE),
                       list(ncov5_df, nsig_df, gini_df, dendrogram_df))

  ## --- Outlier calls on the merged columns ---
  # NSig80 IQR rule (3*IQR)
  Q1_nsig <- quantile(summary_df$NSig80, 0.25, na.rm = TRUE)
  Q3_nsig <- quantile(summary_df$NSig80, 0.75, na.rm = TRUE)
  IQR_nsig <- IQR(summary_df$NSig80, na.rm = TRUE)
  summary_df$Outlier_NSig80 <- as.integer(
    summary_df$NSig80 < (Q1_nsig - 3 * IQR_nsig) |
      summary_df$NSig80 > (Q3_nsig + 3 * IQR_nsig)
  )

  # HighCoverageGeneCount IQR rule (3*IQR)
  Q1_cov <- quantile(summary_df$HighCoverageGeneCount, 0.25, na.rm = TRUE)
  Q3_cov <- quantile(summary_df$HighCoverageGeneCount, 0.75, na.rm = TRUE)
  IQR_cov <- IQR(summary_df$HighCoverageGeneCount, na.rm = TRUE)
  summary_df$Outlier_HighCoverageGeneCount <- as.integer(
    summary_df$HighCoverageGeneCount < (Q1_cov - 3 * IQR_cov) |
      summary_df$HighCoverageGeneCount > (Q3_cov + 3 * IQR_cov)
  )

  # Gini hard cutoff
  summary_df$Outlier_Gini <- as.integer(summary_df$Gini > 0.95)

  # Dendrogram distance threshold (0..1 scale)
  summary_df$Outlier_Dendrogram <- as.integer(summary_df$Dendrogram_Distance > threshold)

  ## --- finalize & register ---
  # Optional: stable ordering by Sample
  summary_df <- summary_df[order(summary_df$Sample), , drop = FALSE]

  dataSet$summary_df <- summary_df
  RegisterData(dataSet)

  return(1)
}


# -------------------------------------------------------------------------
#  qc.giniplot.json()
#  ------------------------------------------------------------------------
#  gini_df      data.frame with columns: Sample, Gini, Status
#  imgNm        stem for the JSON file (".json" is appended automatically)
#  threshold    horizontal dashed reference line
#  jitter.w     half-width of horizontal jitter (0–0.5 recommended)
# -------------------------------------------------------------------------
qc.giniplot.json <- function(gini_df,
                             imgNm     = "Gini_plot",
                             threshold = 0.95,
                             jitter.w  = 0.45) {

  stopifnot(all(c("Sample", "Gini", "Status") %in% names(gini_df)))

  ## 1 · Tukey fences & statistical-outlier flag -------------------------
  stats      <- boxplot.stats(gini_df$Gini, coef = 1.5)$stats
  q1         <- stats[2]; q3 <- stats[4]; iqr <- q3 - q1
  lowFence   <- q1 - 1.5 * iqr
  highFence  <- q3 + 1.5 * iqr
  gini_df$stat_out <- with(gini_df, Gini < lowFence | Gini > highFence)

  ## 2 · Semantic palette (Normal / Outlier) -----------------------------
  status_cols <- c(Normal = "#666666", Outlier = "#E41A1C")

  ## 3 · Traces ----------------------------------------------------------
  # 3a · Box built from in-fence values only
  tr_box <- list(
    x              = rep(0, sum(!gini_df$stat_out)),
    y              = I(gini_df$Gini[!gini_df$stat_out]),
    quartilemethod = "linear",
    type           = "box",
    width          = 0.8,
    name           = "",
    boxpoints      = FALSE,
    fillcolor      = "rgba(200,200,200,0.6)",
    line           = list(color = "#000000"),
    hoverinfo      = "skip",
    showlegend     = FALSE
  )

  # 3b · Invisible all-points trace (for autoscale)
  set.seed(1)
  tr_all <- list(
    x          = I(runif(nrow(gini_df), -jitter.w, jitter.w)),
    y          = I(gini_df$Gini),
    type       = "scatter",
    mode       = "markers",
    marker     = list(color = "rgba(0,0,0,0)", size = 0),
    hoverinfo  = "skip",
    showlegend = FALSE
  )

  # 3c · Visible points (semantic colour, stat-outlines)
  show_labels <- nrow(gini_df) <= 20
  set.seed(2)
  points_trace <- list(
    x    = I(runif(nrow(gini_df), -jitter.w, jitter.w)),
    y    = I(gini_df$Gini),
    type = "scatter",
    mode = if (show_labels) "markers+text" else "markers",
    text = if (show_labels) gini_df$Sample else "",
    textposition = "right",
    name = "Samples",
    hoverinfo = "text",
    hovertext = paste0(
      "Sample: ", gini_df$Sample,
      "<br>Gini: ", signif(gini_df$Gini, 3),
      "<br>Status: ", gini_df$Status
    ),
    marker = list(
      color = status_cols[gini_df$Status],
      size  = 8,
      line  = list(
        color = ifelse(gini_df$stat_out, "black", "rgba(0,0,0,0)"),
        width = ifelse(gini_df$stat_out, 1, 0)
      )
    ),
    showlegend = FALSE
  )

  traces <- list(tr_box, tr_all, points_trace)

  ## 4 · Layout ----------------------------------------------------------
  layout <- list(
    plot_bgcolor  = "#FFFFFF",
    paper_bgcolor = "#FFFFFF",
    xaxis = list(
      title          = "",
      range          = c(-jitter.w - 0.1, jitter.w + 0.1),
      zeroline       = FALSE,
      showticklabels = FALSE,
      showline       = TRUE,
      linecolor      = "#000000"
    ),
    yaxis = list(
      title     = list(text = "Gini coefficient"),
      zeroline  = FALSE,
      ticks     = "outside",
      showline  = TRUE,
      linecolor = "#000000",
      showgrid  = TRUE,
      gridcolor = "rgba(200,200,200,0.4)"
    ),
    shapes = list(list(
      type  = "line",
      xref  = "paper", x0 = 0, x1 = 1,
      yref  = "y",     y0 = threshold, y1 = threshold,
      line  = list(color = "#0026FF", dash = "dot")
    )),
    legend = list(
      title       = list(text = "Sample Status"),
      orientation = "v",
      x = 1.02, y = 1,
      xanchor = "left", yanchor = "top"
    ),
    margin = list(l = 60, r = 110, t = 20, b = 40)
  )

  ## 5 · Write JSON ------------------------------------------------------
  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE, digits = 16
  )

  invisible("NA")
}


# -------------------------------------------------------------------------
#  dendro_df   data.frame with columns:
#              Sample, Dendrogram_Distance, Status (Normal / Outlier),
#              LabelMe (TRUE/FALSE → label on plot)
#  imgNm       stem for JSON file ("<imgNm>.json")
#  threshold   horizontal dashed cut-off
#  jitter.w    half-width for horizontal jitter of points
# -------------------------------------------------------------------------
qc.dendrogram.json <- function(dendro_df,
                               imgNm     = "Dendrogram_plot",
                               threshold = 0.10,
                               jitter.w  = 0.45) {

  stopifnot(all(c("Sample", "Dendrogram_Distance", "Status", "LabelMe") %in% names(dendro_df)))

  ## ── 1 · Tukey fences and statistical-outlier flag -------------------
  stats     <- boxplot.stats(dendro_df$Dendrogram_Distance, coef = 1.5)$stats
  q1        <- stats[2];  q3 <- stats[4];  iqr <- q3 - q1
  lowFence  <- q1 - 1.5 * iqr
  highFence <- q3 + 1.5 * iqr
  dendro_df$stat_out <- with(dendro_df,
                             Dendrogram_Distance < lowFence |
                             Dendrogram_Distance > highFence)

  ## ── 2 · Semantic palette -------------------------------------------
  status_cols <- c(Normal = "#666666", Outlier = "#E41A1C")

  ## ── 3 · Traces ------------------------------------------------------
  # 3a · Box from in-fence points
  tr_box <- list(
    x              = rep(0, sum(!dendro_df$stat_out)),
    y              = I(dendro_df$Dendrogram_Distance[!dendro_df$stat_out]),
    quartilemethod = "linear",
    type           = "box",
    width          = 0.8,
    name           = "",
    boxpoints      = FALSE,
    fillcolor      = "rgba(200,200,200,0.6)",
    line           = list(color = "#000000"),
    hoverinfo      = "skip",
    showlegend     = FALSE
  )

  # 3b · Invisible all-points scatter (autoscale helper)
  set.seed(1)
  tr_all <- list(
    x          = I(runif(nrow(dendro_df), -jitter.w, jitter.w)),
    y          = I(dendro_df$Dendrogram_Distance),
    type       = "scatter",
    mode       = "markers",
    marker     = list(color = "rgba(0,0,0,0)", size = 0),
    hoverinfo  = "skip",
    showlegend = FALSE
  )

  # 3c · Visible points (semantic colouring, stat outlines)
  set.seed(2)
  points_trace <- list(
    x    = I(runif(nrow(dendro_df), -jitter.w, jitter.w)),
    y    = I(dendro_df$Dendrogram_Distance),
    type = "scatter",
    mode = "markers+text",
    text = ifelse(dendro_df$LabelMe, dendro_df$Sample, ""),
    textposition = "right",
    name = "Samples",
    hoverinfo = "text",
    hovertext = paste0(
      "Sample: ", dendro_df$Sample,
      "<br>Distance: ", signif(dendro_df$Dendrogram_Distance, 3),
      "<br>Status: ", dendro_df$Status
    ),
    marker = list(
      color = status_cols[dendro_df$Status],
      size  = 8,
      line  = list(
        color = ifelse(dendro_df$stat_out, "black", "rgba(0,0,0,0)"),
        width = ifelse(dendro_df$stat_out, 1, 0)
      )
    ),
    showlegend = FALSE
  )

  traces <- list(tr_box, tr_all, points_trace)

  ## ── 4 · Layout ------------------------------------------------------
  layout <- list(
    plot_bgcolor  = "#FFFFFF",
    paper_bgcolor = "#FFFFFF",
    xaxis = list(
      title          = "",
      range          = c(-jitter.w - 0.1, jitter.w + 0.1),
      zeroline       = FALSE,
      showticklabels = FALSE,
      showline       = TRUE,
      linecolor      = "#000000"
    ),
    yaxis = list(
      title     = list(text = "Max pair-wise distance (1 \u2212 Pearson \u03c1)"),
      zeroline  = FALSE,
      ticks     = "outside",
      showline  = TRUE,
      linecolor = "#000000",
      showgrid  = TRUE,
      gridcolor = "rgba(200,200,200,0.4)"
    ),
    shapes = list(list(
      type  = "line",
      xref  = "paper", x0 = 0, x1 = 1,
      yref  = "y",     y0 = threshold, y1 = threshold,
      line  = list(color = "#0026FF", dash = "dot")
    )),
    legend = list(
      title       = list(text = "Sample Status"),
      orientation = "v",
      x = 1.02, y = 1,
      xanchor = "left", yanchor = "top"
    ),
    margin = list(l = 70, r = 110, t = 20, b = 40)
  )

  ## ── 5 · Write JSON ---------------------------------------------------
  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE, digits = 16
  )
  invisible("NA")
}

qc.ncov5plot.json <- function(ncov5_df,
                              imgNm    = "NCov5_plot",
                              lower,
                              upper,
                              jitter.w = 0.45) {

  stopifnot(all(c("Sample", "HighCoverageGeneCount", "Status") %in% names(ncov5_df)),
            is.numeric(lower), length(lower) == 1,
            is.numeric(upper), length(upper) == 1)

  ## ── 1 · Tukey fences & statistical-outlier flag --------------------
  stats      <- boxplot.stats(ncov5_df$HighCoverageGeneCount, coef = 1.5)$stats
  q1         <- stats[2]; q3 <- stats[4]; iqr <- q3 - q1
  lowFence   <- q1 - 1.5 * iqr
  highFence  <- q3 + 1.5 * iqr
  ncov5_df$stat_out <- with(ncov5_df,
                             HighCoverageGeneCount < lowFence |
                             HighCoverageGeneCount > highFence)

  ## ── 2 · palettes ----------------------------------------------------
  stat_cols <- c(Normal = "#666666", Outlier = "#E41A1C")

  ## ── 3 · Traces ------------------------------------------------------
  # 3a · box (only in-fence points)
  tr_box <- list(
    x              = rep(0, sum(!ncov5_df$stat_out)),
    y              = I(ncov5_df$HighCoverageGeneCount[!ncov5_df$stat_out]),
    quartilemethod = "linear",
    type           = "box",
    width          = 0.8,
    name           = "",
    boxpoints      = FALSE,
    fillcolor      = "rgba(200,200,200,0.6)",
    line           = list(color = "#000000"),
    hoverinfo      = "skip",
    showlegend     = FALSE
  )

  # 3b · invisible all-points scatter (for autoscale)
  set.seed(1)
  tr_all <- list(
    x          = I(runif(nrow(ncov5_df), -jitter.w, jitter.w)),
    y          = I(ncov5_df$HighCoverageGeneCount),
    type       = "scatter",
    mode       = "markers",
    marker     = list(color = "rgba(0,0,0,0)", size = 0),
    hoverinfo  = "skip",
    showlegend = FALSE
  )

  # 3c · labelled points (semantic status colouring, stat-outlines)
  # 3c · visible points (one trace for all samples)
set.seed(2)
points_trace <- list(
  x    = I(runif(nrow(ncov5_df), -jitter.w, jitter.w)),
  y    = I(ncov5_df$HighCoverageGeneCount),
  type = "scatter",
  mode = if (any(ncov5_df$stat_out)) "markers+text" else "markers",
  text = ifelse(ncov5_df$stat_out, ncov5_df$Sample, ""),
  textposition = "right",
  name = "Samples",
  hoverinfo = "text",
  hovertext = paste0(
    "Sample: ", ncov5_df$Sample,
    "<br>Count: ", ncov5_df$HighCoverageGeneCount,
    "<br>Status: ", ncov5_df$Status
  ),
  marker = list(
    color = stat_cols[ncov5_df$Status],           # vector OK
    size  = 8,
    line  = list(
      color = ifelse(ncov5_df$stat_out, "black", "rgba(0,0,0,0)"),
      width = ifelse(ncov5_df$stat_out, 1, 0)
    )
  ),
  showlegend = FALSE
)

  traces <- list(tr_box, tr_all, points_trace)

  ## ── 4 · Layout ------------------------------------------------------
  layout <- list(
    plot_bgcolor  = "#FFFFFF",
    paper_bgcolor = "#FFFFFF",
    xaxis = list(
      title          = "",
      range          = c(-jitter.w - 0.1, jitter.w + 0.1),
      zeroline       = FALSE,
      showticklabels = FALSE,
      showline       = TRUE,
      linecolor      = "#000000"
    ),
    yaxis = list(
      title     = list(text = "Genes with > 5 uniquely mapped reads"),
      zeroline  = FALSE,
      ticks     = "outside",
      showline  = TRUE,
      linecolor = "#000000",
      showgrid  = TRUE,
      gridcolor = "rgba(200,200,200,0.4)"
    ),
    shapes = list(
      list(type="line", xref="paper", x0=0, x1=1,
           yref="y", y0=lower, y1=lower,
           line=list(color="#0026FF", dash="dot")),
      list(type="line", xref="paper", x0=0, x1=1,
           yref="y", y0=upper, y1=upper,
           line=list(color="#0026FF", dash="dot"))
    ),
    legend = list(
      title       = list(text = "Sample Status"),
      orientation = "v",
      x = 1.02, y = 1,
      xanchor = "left", yanchor = "top"
    ),
    margin = list(l = 70, r = 110, t = 20, b = 40)
  )

  ## ── 5 · Write JSON ---------------------------------------------------
  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE, digits = 16
  )

  invisible("NA")
}

qc.nsigplot.json <- function(nsig_df,
                             imgNm    = "NSig80_plot",
                             lower,
                             upper,
                             jitter.w = 0.45) {

  stopifnot(all(c("Sample", "NSig80", "outlier") %in% names(nsig_df)))

  ## ------------------------------------------------------------------
  ## 1 · Compute Tukey fences and flag statistical outliers
  ## ------------------------------------------------------------------
  stats       <- boxplot.stats(nsig_df$NSig80, coef = 1.5)$stats
  q1          <- stats[2]; q3 <- stats[4]; iqr <- q3 - q1
  lowFence    <- q1 - 1.5 * iqr
  highFence   <- q3 + 1.5 * iqr
  nsig_df$stat_out <- with(nsig_df, NSig80 < lowFence | NSig80 > highFence)

  ## ------------------------------------------------------------------
  ## 2 · Palette for semantic status (Normal / Outlier)
  ## ------------------------------------------------------------------
  status_cols <- c(Normal = "#666666", Outlier = "#E41A1C")

  ## ------------------------------------------------------------------
  ## 3 · Plotly traces
  ## ------------------------------------------------------------------
  # 3a ─ Box: only in-fence points
  tr_box <- list(
    x              = rep(0, sum(!nsig_df$stat_out)),
    y              = I(nsig_df$NSig80[!nsig_df$stat_out]),
    quartilemethod = "linear",
    type           = "box",
    width          = 0.8,
    name           = "",
    boxpoints      = FALSE,
    fillcolor      = "rgba(200,200,200,0.6)",
    line           = list(color = "#000000"),
    hoverinfo      = "skip",
    showlegend     = FALSE
  )

  # 3b ─ Invisible “all” trace for autoscaling
  set.seed(1)
  tr_all <- list(
    x          = I(runif(nrow(nsig_df), -jitter.w, jitter.w)),
    y          = I(nsig_df$NSig80),
    type       = "scatter",
    mode       = "markers",
    marker     = list(color = "rgba(0,0,0,0)", size = 0),
    hoverinfo  = "skip",
    showlegend = FALSE
  )

  # 3c ─ Points (semantic status colouring, stat-outliers outlined)
  set.seed(2)
  points_trace <- list(
    x    = I(runif(nrow(nsig_df), -jitter.w, jitter.w)),
    y    = I(nsig_df$NSig80),
    type = "scatter",
    mode = "markers+text",
    text = ifelse(nsig_df$stat_out, nsig_df$Sample, ""),
    textposition = "right",
    name = "Samples",
    hoverinfo = "text",
    hovertext = paste0(
      "Sample: ", nsig_df$Sample,
      "<br>NSig80: ", nsig_df$NSig80,
      "<br>Status: ", nsig_df$outlier
    ),
    marker = list(
      color = status_cols[nsig_df$outlier],
      size  = 8,
      line  = list(
        color = ifelse(nsig_df$stat_out, "black", "rgba(0,0,0,0)"),
        width = ifelse(nsig_df$stat_out, 1, 0)
      )
    ),
    showlegend = FALSE
  )

  traces <- list(tr_box, tr_all, points_trace)

  ## ------------------------------------------------------------------
  ## 4 · Layout (unchanged)
  ## ------------------------------------------------------------------
  layout <- list(
    plot_bgcolor  = "#FFFFFF", paper_bgcolor = "#FFFFFF",
    xaxis = list(title="", range=c(-jitter.w-0.1, jitter.w+0.1),
                 zeroline=FALSE, showticklabels=FALSE,
                 showline=TRUE, linecolor="#000"),
    yaxis = list(title=list(text="NSig80 (genes reaching 80% of signal)"),
                 zeroline=FALSE, ticks="outside", showline=TRUE,
                 linecolor="#000", showgrid=TRUE,
                 gridcolor="rgba(200,200,200,0.4)"),
    shapes = list(
      list(type="line", xref="paper", x0=0, x1=1,
           yref="y", y0=lower, y1=lower,
           line=list(color="#0026FF", dash="dot")),
      list(type="line", xref="paper", x0=0, x1=1,
           yref="y", y0=upper, y1=upper,
           line=list(color="#0026FF", dash="dot"))
    ),
    legend = list(title=list(text="Sample Status"),
                  orientation="v", x=1.02, y=1,
                  xanchor="left", yanchor="top"),
    margin = list(l=70, r=110, t=20, b=40)
  )

  ## ------------------------------------------------------------------
  ## 5 · Write JSON
  ## ------------------------------------------------------------------
  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE, digits = 16
  )
}
qc.pcaplot.outliers.json <- function(dataSet, x, imgNm,
                                     uniq_map_col = "uniq_map",
                                     min_per_dose = 2,
                                     min_vehicles = 3) {
  jsonFile <- paste0(imgNm, ".json")
  csvFile  <- paste0(imgNm, "_outliers.csv")

  require(plotly)
  require(jsonlite)  # Use jsonlite instead of rjson for proper array handling

  # ----- load PCA & align to meta -----
  analSet <- readSet(analSet, "analSet")
  pca     <- analSet$pca
  xlabel  <- pca$xlabel
  ylabel  <- pca$ylabel

  pca.res <- as.data.frame(pca$x)[, 1:2, drop = FALSE]
  colnames(pca.res) <- c("PC1","PC2")

  # Find common samples between PCA and metadata (matching PNG generation logic)
  common_samples <- intersect(rownames(pca.res), rownames(dataSet$meta.info))
  if (length(common_samples) == 0) {
    stop("[qc.pcaplot.outliers.json] ERROR: No common samples between PCA results and metadata!")
  }

  # Align both PCA and metadata to common samples (instead of using match which can create NAs)
  pca.res <- pca.res[common_samples, , drop = FALSE]
  meta <- dataSet$meta.info[common_samples, , drop = FALSE]

  pca.res$sample_id <- rownames(pca.res)
  stopifnot(nrow(meta) == nrow(pca.res))
  pca.res$group <- as.character(meta[[1]])

  doShape <- FALSE; shape.map <- NULL; shape.levels <- NULL
  if (ncol(meta) >= 2) {
    second <- meta[[2]]
    isDisc <- !is.numeric(second)
    levs   <- unique(as.character(second))
    if (isDisc && length(levs) <= 6) {
      doShape <- TRUE
      pca.res$shape <- as.character(second)
      symbols <- c("circle","square","diamond","cross","x","triangle-up","triangle-down","star")
      shape.map <- setNames(symbols[seq_along(levs)], levs)
      shape.levels <- levs
    }
  }

  nR <- nrow(meta)
  pca.res$dose <- if ("dose" %in% colnames(meta)) as.character(meta[["dose"]]) else rep(NA_character_, nR)
  pca.res$is_vehicle <- if ("is_vehicle" %in% colnames(meta)) {
    as.logical(as.character(meta[["is_vehicle"]]))
  } else rep(FALSE, nR)
  pca.res$uniq_map <- if (uniq_map_col %in% colnames(meta)) as.numeric(meta[[uniq_map_col]]) else rep(NA_real_, nR)

  # ----- color mapping (your palettes) -----
  paramSet <- readSet(paramSet, "paramSet")
  unique_grps <- unique(pca.res$group)
  if (grepl("norm", imgNm) && !is.null(paramSet$oneDataAnalType) && paramSet$oneDataAnalType == "dose") {
    pal <- colorRampPalette(c("#2196F3", "#DE690D"))(length(unique_grps))
  } else {
    okabe <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
    pal <- rep(okabe, length.out = length(unique_grps))
  }
  col.map <- setNames(pal, unique_grps)

  axis_class <- function(vals_named) {
    labs <- names(vals_named); vals <- as.numeric(vals_named)
    n <- length(vals)
    cls <- rep("none", n)

    med <- median(vals, na.rm = TRUE)
    md  <- mad(vals, constant = 1, na.rm = TRUE)
    if (is.na(md) || md == 0) {
      q <- quantile(vals, probs = c(0.10, 0.90), na.rm = TRUE, names = FALSE)
      core_idx <- vals >= q[1] & vals <= q[2]
    } else {
      core_idx <- abs(vals - med) <= 3 * md
    }
    if (sum(core_idx, na.rm = TRUE) < 3) {
      ord <- order(vals); k1 <- max(1, floor(0.10 * n)); k2 <- min(n, ceiling(0.90 * n))
      core_idx <- FALSE; core_idx[ord[k1:k2]] <- TRUE
    }

    core_vals <- vals[core_idx]
    span_core <- max(core_vals, na.rm = TRUE) - min(core_vals, na.rm = TRUE)
    if (!is.finite(span_core) || span_core == 0) span_core <- diff(range(vals, na.rm = TRUE))

    for (k in seq_len(n)) {
      xi <- vals[k]
      if (core_idx[k]) { cls[k] <- "none"; next }
      sep_core <- min(abs(xi - core_vals))
      if (sep_core > 2 * span_core)      cls[k] <- "strong"
      else if (sep_core > 1 * span_core) cls[k] <- "moderate"
      else                                cls[k] <- "none"
    }
    setNames(cls, labs)
  }

  euclid_med <- function(df) {
    cx <- median(df$PC1, na.rm = TRUE); cy <- median(df$PC2, na.rm = TRUE)
    sqrt((df$PC1 - cx)^2 + (df$PC2 - cy)^2)
  }
  within_dose_far <- function(df) {
    if (all(is.na(df$dose))) return(rep(FALSE, nrow(df)))
    unlist(
      by(df, df$dose, function(dd) {
        if (nrow(dd) < 2) return(rep(NA, nrow(dd)))
        cx <- median(dd$PC1); cy <- median(dd$PC2)
        d  <- sqrt((dd$PC1-cx)^2 + (dd$PC2-cy)^2)
        thr <- median(d, na.rm = TRUE) * 2
        d > thr
      }),
      use.names = FALSE
    )
  }

  .safe_chr <- function(x) { if (length(x) == 0 || is.na(x) || x == "NA") "" else as.character(x) }
  .safe_reason <- function(x) { if (length(x) == 0 || is.na(x) || x == "NA" || x == "") "" else as.character(x) }
  .make_customdata <- function(subdf) {
    n <- nrow(subdf)
    out <- vector("list", n)
    for (i in seq_len(n)) {
      out[[i]] <- list(
        dose       = .safe_chr(subdf$dose[i]),
        is_vehicle = .safe_chr(subdf$is_vehicle[i]),
        reason     = .safe_reason(subdf$reason[i]),
        status     = .safe_chr(subdf$.__status__[i])
      )
    }
    out
  }

  # trace appending guards
  .is_single_trace <- function(x) is.list(x) && !is.null(x$type) && !is.null(x$x) && !is.null(x$y)
  .append_trace <- function(dst, tr) {
    if (is.null(tr) || !is.list(tr)) return(dst)
    if (.is_single_trace(tr)) { dst[[length(dst) + 1]] <- tr; return(dst) }
    for (k in seq_along(tr)) {
      tk <- tr[[k]]
      if (.is_single_trace(tk)) dst[[length(dst) + 1]] <- tk
    }
    dst
  }

  df <- pca.res
  ax1 <- axis_class(setNames(df$PC1, df$sample_id))
  ax2 <- axis_class(setNames(df$PC2, df$sample_id))
  df$ax_PC1 <- ax1[df$sample_id]
  df$ax_PC2 <- ax2[df$sample_id]
  df$axis_class <- ifelse(df$ax_PC1 == "strong" | df$ax_PC2 == "strong", "strong",
                          ifelse(df$ax_PC1 == "moderate" | df$ax_PC2 == "moderate", "moderate", "none"))
  df$moderate_both_axes <- (df$ax_PC1 == "moderate" & df$ax_PC2 == "moderate")

  D  <- euclid_med(df)
  df$far_euclid <- D > (2 * median(D, na.rm = TRUE))
  df$far_repl  <- within_dose_far(df)

  df$reason  <- NA_character_
  df$exclude <- df$axis_class == "strong"
  df$reason[df$exclude] <- "Strong axis separation vs. core (>2× core span)"

  if (!all(is.na(df$dose))) {
    strong_rows <- which(df$exclude)
    if (length(strong_rows) > 1) {
      dups <- duplicated(df$dose[strong_rows]) | duplicated(df$dose[strong_rows], fromLast = TRUE)
      if (any(dups, na.rm = TRUE)) {
        df$exclude[strong_rows] <- FALSE
        df$reason[strong_rows]  <- NA_character_
      }
    }
  }

  m_idx <- which(!df$exclude & df$axis_class == "moderate")
  for (i in m_idx) {
    reasons <- character(0)
    if (isTRUE(df$moderate_both_axes[i]) && isTRUE(df$far_euclid[i]))
      reasons <- c(reasons, "Moderate on both axes with large Euclidean distance")
    if (isTRUE(df$far_repl[i]))
      reasons <- c(reasons, "Far from replicate/similar dose cluster")
    if (isTRUE(df$worse_qc[i]))
      reasons <- c(reasons, "Lower sequencing quality")
    if (length(reasons)) {
      df$exclude[i] <- TRUE
      df$reason[i]  <- paste(reasons, collapse = "; ")
    }
  }

  if (!all(is.na(df$dose))) {
    kept_by_dose <- tapply(!df$exclude & !df$is_vehicle, df$dose, sum)
    drop_doses <- names(kept_by_dose[!is.na(kept_by_dose) & kept_by_dose < min_per_dose])
    if (length(drop_doses)) {
      hit <- which(df$dose %in% drop_doses & !df$is_vehicle)
      df$exclude[hit] <- TRUE
      df$reason[hit]  <- ifelse(is.na(df$reason[hit]),
                                "Dose dropped (<2 samples after QC/outlier)",
                                paste(df$reason[hit], "Dose dropped (<2 samples)", sep="; "))
    }
  }

  veh_kept <- sum(!df$exclude & df$is_vehicle, na.rm = TRUE)
  vehicle_note <- if (any("is_vehicle" == colnames(meta)) && veh_kept < min_vehicles)
    sprintf("Warning: only %d vehicle samples kept (< %d).", veh_kept, min_vehicles) else NULL

  status_lab <- ifelse(df$exclude, "Excluded",
                       ifelse(df$axis_class == "moderate", "Moderate", "Kept"))
  df$.__status__ <- status_lab

  # Add a default reason for moderate outliers when no specific reason is set.
  mod_idx <- which(is.na(df$reason) | df$reason == "")
  mod_idx <- mod_idx[df$.__status__[mod_idx] == "Moderate"]
  if (length(mod_idx)) {
    df$reason[mod_idx] <- ifelse(df$moderate_both_axes[mod_idx],
                                 "Moderate separation on both axes vs. core",
                                 "Moderate separation on one axis vs. core")
  }

  status_styles <- list(
    Kept     = list(line = list(color = "white", width = 0.5), size = 8,  opacity = 0.9),
    Moderate = list(line = list(color = "orange", width = 2),  size = 9,  opacity = 1.0),
    Excluded = list(line = list(color = "red",    width = 3),  size = 10, opacity = 1.0)
  )

  # Create one trace per group, with per-sample marker styling for outliers
  traces <- list()
  for (g in unique_grps) {
    gdf <- df[df$group == g, , drop = FALSE]
    if (nrow(gdf) == 0) next

    # Build per-sample marker properties based on status
    marker_colors <- rep(col.map[[g]], nrow(gdf))
    marker_sizes <- rep(8, nrow(gdf))
    marker_line_colors <- rep("white", nrow(gdf))
    marker_line_widths <- rep(0.5, nrow(gdf))
    marker_opacities <- rep(0.9, nrow(gdf))
    text_labels <- rep("", nrow(gdf))

    for (i in seq_len(nrow(gdf))) {
      status <- gdf$.__status__[i]
      if (status == "Moderate") {
        marker_line_colors[i] <- "orange"
        marker_line_widths[i] <- 2
        marker_sizes[i] <- 9
        marker_opacities[i] <- 1.0
        text_labels[i] <- gdf$sample_id[i]
      } else if (status == "Excluded") {
        marker_line_colors[i] <- "red"
        marker_line_widths[i] <- 3
        marker_sizes[i] <- 10
        marker_opacities[i] <- 1.0
        text_labels[i] <- gdf$sample_id[i]
      }
    }

    # Build hover text with sample name, dose, status, and reason
    hover_texts <- sapply(seq_len(nrow(gdf)), function(j) {
      parts <- c(paste0("Sample: ", gdf$sample_id[j]))
      if (!is.na(gdf$dose[j]) && gdf$dose[j] != "") parts <- c(parts, paste0("Dose: ", gdf$dose[j]))
      if (!is.na(gdf$is_vehicle[j]) && gdf$is_vehicle[j]) parts <- c(parts, "Vehicle")
      if (gdf$.__status__[j] != "Kept") parts <- c(parts, paste0("Status: ", gdf$.__status__[j]))
      if (!is.na(gdf$reason[j]) && gdf$reason[j] != "") parts <- c(parts, paste0("Reason: ", gdf$reason[j]))
      paste(parts, collapse = "<br>")
    })

    # Determine if we should show text labels (only for outliers)
    has_outliers <- any(gdf$.__status__ != "Kept")
    mode_val <- if (has_outliers) "markers+text" else "markers"

    # Wrap single values to ensure they become arrays in JSON
    x_val <- if (length(gdf$PC1) == 1) list(gdf$PC1) else gdf$PC1
    y_val <- if (length(gdf$PC2) == 1) list(gdf$PC2) else gdf$PC2

    # Build marker list - arrays should not be auto-unboxed
    mkr <- list(
      color = if (length(unique(marker_colors)) == 1) marker_colors[1] else I(marker_colors),
      size = if (length(unique(marker_sizes)) == 1) marker_sizes[1] else I(marker_sizes),
      opacity = if (length(unique(marker_opacities)) == 1) marker_opacities[1] else I(marker_opacities),
      line = list(
        color = if (length(unique(marker_line_colors)) == 1) marker_line_colors[1] else I(marker_line_colors),
        width = if (length(unique(marker_line_widths)) == 1) marker_line_widths[1] else I(marker_line_widths)
      )
    )

    traces[[length(traces) + 1]] <- list(
      x = x_val,
      y = y_val,
      type = "scatter",
      mode = mode_val,
      name = g,
      showlegend = TRUE,
      marker = mkr,
      text = if (all(text_labels == "")) NULL else I(text_labels),
      hovertext = I(hover_texts),
      customdata = .make_customdata(gdf),
      hoverinfo = "text",
      textposition = "top center"
    )
  }

  if (doShape) {
    for (sh in shape.levels) {
      traces[[length(traces) + 1]] <- list(
        x = c(NA), y = c(NA), type = "scatter", mode = "markers",
        name = paste0("Shape: ", sh), showlegend = TRUE,
        marker = list(symbol = shape.map[[sh]], color = "black", size = 8)
      )
    }
  }

  subtitle <- if (!is.null(vehicle_note)) vehicle_note else ""
  layout <- list(
    title = "",
    xaxis = list(title = xlabel),
    yaxis = list(title = ylabel),
    legend = list(orientation = "v", x = 1.02, y = 1, xanchor = "left", yanchor = "top"),
    `shape.map` = shape.map,          # keep for JS legend helpers if you use them
    meta2Name = if (doShape) names(meta)[2] else NULL,
    annotations = if (nzchar(subtitle)) list(list(
      x = 0, y = 1.08, xref = "paper", yref = "paper",
      xanchor = "left", yanchor = "bottom",
      text = subtitle, showarrow = FALSE, font = list(size = 12)
    )) else NULL
  )

  plot_data <- list(data = traces, layout = layout)
  # Use jsonlite with auto_unbox=FALSE to ensure single-element vectors become arrays, not scalars
  json.obj  <- jsonlite::toJSON(plot_data, auto_unbox = TRUE, digits = NA, na = "null")
  writeLines(json.obj, jsonFile)

  out_cols <- c("sample_id","group","dose","is_vehicle","uniq_map",
                "PC1","PC2","ax_PC1","ax_PC2","axis_class",
                "moderate_both_axes","far_euclid","far_repl",
                "worse_qc","exclude","reason","__.__status__")
  keep_cols <- intersect(out_cols, colnames(df))
  out_tab   <- df[, keep_cols, drop = FALSE]
  colnames(out_tab)[colnames(out_tab) == ".__status__"] <- "status"
  utils::write.csv(out_tab, file = csvFile, row.names = FALSE)

  return("NA")
}
