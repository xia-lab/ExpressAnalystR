##################################################
## R script for ExpressAnalyst
## Description: functions for quality check boxplot
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################


#'Perform data normalization
#'@description Filtering and Normalizing gene expression data
#'@param norm.opt Normalization method to be used
#'@param var.thresh Variance threshold
#'@param abundance Relative abundance threshold
#'@param count.thresh Count threshold for RNA-seq data and abundance threshold for microarray data
#'@param filterUnmapped, boolean, whether to filter unmapped genes or not
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
PerformNormalization <- function(dataName, norm.opt, var.thresh, count.thresh, filterUnmapped,
                                 islog = "false", countOpt = "sum") {
  #save.image("norm.RData");
  paramSet <- readSet(paramSet, "paramSet");
  msgSet   <- readSet(msgSet, "msgSet");
  dataSet  <- readDataset(dataName);
  msg <- ""

  ds <- qs::qread("data.raw.qs");
  msgSet$current.msg <- c(msgSet$current.msg, paste0("Diagnostic before filtering: samples=", ncol(ds), 
                                                    " features=", nrow(ds),
                                                    " zeros=", sum(ds==0),
                                                    " min=", min(ds, na.rm=T),
                                                    " max=", max(ds, na.rm=T)))
  saveSet(msgSet,"msgSet")

  # reload the original annotated counts before filtering so each normalization run starts from raw data
  if (file.exists("orig.data.anot.qs")) {
    raw.annot <- qs::qread("orig.data.anot.qs")
  } else if (file.exists("data.raw.qs")) {
    raw.annot <- qs::qread("data.raw.qs")
  } else {
    raw.annot <- dataSet$data.norm
  }
  qs::qsave(raw.annot, file = "data.anot.qs")
  data <- PerformFiltering(dataSet, var.thresh, count.thresh, filterUnmapped, countOpt)
  .save.annotated.data(data)
  msg <- paste(filt.msg, msg)

  data <- sanitizeSmallNumbers(data)
  diag.filtered <- paste0("Diagnostic after filtering: rows=", nrow(data),
                          " cols=", ncol(data),
                          " zeros=", sum(data == 0, na.rm=T))
  msgSet$current.msg <- c(msgSet$current.msg, diag.filtered)
  saveSet(msgSet, "msgSet")

  if (dataSet$type == "prot") {
    diag.prot <- paste0("Diagnostic before normalization (prot): norm.opt=", norm.opt,
                        " rows=", nrow(data), " cols=", ncol(data))
    msgSet$current.msg <- c(msgSet$current.msg, diag.prot)
    saveSet(msgSet, "msgSet")
    if (islog == "true" || norm.opt == "Rlr" || norm.opt == "Loess") {
      data <- NormalizeData(data, "log", "NA", "NA")
      msg  <- paste(norm.msg, msg)
    }
  }

  paramSet$norm.opt   <- norm.opt
  paramSet$var.perc   <- var.thresh
  paramSet$abun.perc  <- count.thresh

  if (identical(norm.opt, "MORlog")) {

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      AddErrMsg("MORlog normalization requires the 'DESeq2' package. Please install it.")
      return(0)
    }

    m <- as.matrix(data)

    # basic checks for counts
    if (any(m < 0, na.rm = TRUE)) {
      AddErrMsg("MORlog expects non-negative count data.")
      return(0)
    }
    # ensure integer-like counts for DESeq2
    if (!is.integer(m)) {
      m <- round(m)
    }

    # minimal colData (no outcome needed for size factors)
    # rownames must match sample names (columns of count matrix)
    cd <- S4Vectors::DataFrame(row.names = colnames(m))

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = cd, design = ~ 1)
    dds <- DESeq2::estimateSizeFactors(dds)
    norm_counts <- DESeq2::counts(dds, normalized = TRUE)
    data <- log2(norm_counts + 1)

    msg <- paste("[MORlog] Applied DESeq2 median-of-ratios size-factor normalization and log2(x+1).", msg)
  } else {
    # ---- existing generic normalization ----
    data <- NormalizeData(data, norm.opt, "NA", "NA")
    msg  <- paste(norm.msg, msg)
  }

  if (paramSet$oneDataAnalType == "dose" && min(data) < 0) {
    add.val <- abs(min(data)) + 0.05 * abs(min(data))
    data <- data + add.val
  }

  dataSet$data.norm <- data
  fast.write(sanitizeSmallNumbers(data), file = "data_normalized.csv")
  qs::qsave(data, file = "data.stat.qs")

  msgSet$current.msg <- msg
  saveSet(msgSet,   "msgSet")
  saveSet(paramSet, "paramSet")
  return(RegisterData(dataSet))
}


PerformFiltering <- function(dataSet, var.thresh, count.thresh, filterUnmapped, countOpt){
  msg <- "";
  if(filterUnmapped == "false"){
    # need to update those with annotations
    data1 <- qs::qread("data.raw.qs");
    colnames(data1) <- colnames(dataSet$data.norm)
    anot.id <- qs::qread("annotation.qs");
    hit.inx <- !is.na(anot.id);
    rownames(data1)[hit.inx] <- anot.id[hit.inx];
    res <- RemoveDuplicates(data1, "mean", quiet=T, paramSet, msgSet);
    data1 <- res[[1]];
    msgSet <- res[[2]];
    raw.data.anot <- data1;
    msg <- "Only features with annotations are kept for further analysis.";
  }else{
    if(dataSet$type=="prot"){
     raw.data.anot <- qs::qread("data.missed.qs");
    }else{
     raw.data.anot <- qs::qread("orig.data.anot.qs");
    }
   colnames(raw.data.anot) <- colnames(dataSet$data.norm)
  }
  
  data <- raw.data.anot;
  data<- data[,which(colnames(data)%in% rownames(dataSet$meta.info))]
  # PERFORMANCE FIX (Issue #1): Use vectorized row operations instead of apply()
  # rowSums/rowMeans are 60-100x faster than apply(data, 1, sum/mean)
  # Critical for normalization - affects 100% of users
  if (dataSet$type == "count"){
    if (countOpt == "sum") {
        # Sum approach: sum counts across samples for each gene
        sum.counts <- rowSums(data, na.rm = TRUE)
        rm.inx <- sum.counts < count.thresh
        msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low counts using sum method.", collapse = " ")
    } else if (countOpt == "average") {
        # Average approach: calculate average counts across samples for each gene
        avg.counts <- rowMeans(data, na.rm = TRUE)
        rm.inx <- avg.counts < count.thresh
        msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low counts using average method.", collapse = " ")
    }
  }else{
    avg.signal <- rowMeans(data, na.rm=TRUE)
    abundance.pct <- count.thresh/100;
    p05 <- quantile(avg.signal, abundance.pct)
    all.rows <- nrow(data)
    rm.inx <- avg.signal < p05;
    msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low relative abundance (average expression signal).", collapse=" ");
  }

  if(var.thresh > 0){
  data <- data[!rm.inx,];
  filter.val <- apply(data, 1, IQR, na.rm=T);
  nm <- "Interquantile Range";
  filter.val <- -filter.val
  rk <- rank(filter.val, ties.method='random');
  # remove constant values
  good.inx <- -filter.val > 0;
  kp.pct <- (100 - var.thresh)/100;
  remain <- rk < nrow(data)*kp.pct;
    initial_gene_count <- nrow(data)
  data <- data[remain&good.inx,];
 # Calculate number of genes filtered by IQR
    filtered_by_iqr <- initial_gene_count - nrow(data)

    # Update message with correct number of filtered genes
    filt.msg <<- paste(msg, "Filtered", filtered_by_iqr, "low variance genes based on IQR.")
  }else{
  filt.msg <<- paste(msg, paste("Filtered 0 low variance genes based on IQR"), collapse=" ");
  }
  
  return(data);
}

NormalizeDataMetaMode <-function (nm, opt, colNorm="NA", scaleNorm="NA"){
  if(nm == "NA"){
    paramSet <- readSet(paramSet, "paramSet");
    mdata.all <- paramSet$mdata.all;
    sel.nms <- names(mdata.all);
    for(i in 1:length(sel.nms)){
      dataName <- sel.nms[i];
      dataSet = readDataset(dataName);
      data.filtered <- readDataQs("data.filtered.qs", paramSet$anal.type, dataName);
      data <- NormalizeData(data.filtered,opt, colNorm, scaleNorm);
      if(length(data) == 1){
        return(0);
      }
    dataSet$data.norm <- data;
    dataSet$norm.opt <- opt;
    fast.write(sanitizeSmallNumbers(dataSet$data.norm), file = "data_normalized.csv")
    RegisterData(dataSet);
    }
    return(1)
  }else{
    dataSet <- readDataset(nm);
    data.filtered <- readDataQs("data.filtered.qs", paramSet$anal.type, nm);
    data <- NormalizeData(data.filtered,opt, colNorm, scaleNorm);
    if(length(data) == 1){
      return(0);
    }
    dataSet$data.norm <- data;
    dataSet$norm.opt <- opt;
    qs::qsave(data, file="data.stat.qs");
    return(RegisterData(dataSet));
    
  }
}

NormalizeData <-function (data, norm.opt, colNorm="NA", scaleNorm="NA"){
  msg <- ""
  row.nms <- rownames(data);
  col.nms <- colnames(data);
  msgSet <- readSet(msgSet, "msgSet");
  
  data <- sanitizeSmallNumbers(data)

  # column(sample)-wise normalization
  if(colNorm=="SumNorm"){
    data<-t(apply(data, 2, SumNorm));
    rownm<-"Normalization to constant sum";
  }else if(colNorm=="MedianNorm"){
    data<-t(apply(data, 2, MedianNorm));
    rownm<-"Normalization to sample median";
  }else{
    # nothing to do
    rownm<-"N/A";
  }
  # norm.opt
  if(norm.opt=="log"){
    positiveVals <- data[data > 0]
    if(length(positiveVals) == 0){
      AddErrMsg("All values are non-positive; log normalization cannot proceed.")
      return(0)
    }
    min.pos <- max(min(positiveVals, na.rm=T)/10, 1e-6)
    numberOfNeg = sum(data<0, na.rm = TRUE) + 1; 
    totalNumber = length(data)
    if((numberOfNeg/totalNumber)>0.2){
      msg <- paste(msg, "Can't perform log2 normalization, over 20% of data are negative. Try a different method or maybe the data already normalized?", collapse=" ");
      msgSet$norm.msg <- msgSet$current.msg <- msg;
      saveSet(msgSet, "msgSet");
      return(0);
    }
    data[data<=0] <- min.pos;
    data <- log2(data);
    msg <- paste(msg, "Log2 transformation.", collapse=" ");
  }else if(norm.opt=="vsn"){
    require(limma);
    data <- normalizeVSN(data);
    msg <- paste(msg, "VSN normalization.", collapse=" ");
  }else if(norm.opt=="quantile"){
    require('preprocessCore');
    data <- normalize.quantiles(as.matrix(data), copy=TRUE);
    msg <- paste(msg, "Quantile normalization.", collapse=" ");
  }else if(norm.opt=="combined"){
    require(limma);
    data <- normalizeVSN(data);
    require('preprocessCore');
    data <- normalize.quantiles(as.matrix(data), copy=TRUE);
    msg <- paste(msg, "VSN followed by quantile normalization.", collapse=" ");
  }else if(norm.opt=="logcount"){ # for count data, do it in DE analysis, as it is dependent on design matrix
    require(edgeR);
    nf <- calcNormFactors(data, method = "none");
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- paste(msg, "Limma based on log2-counts per million transformation.", collapse=" ");
  } else if(norm.opt=="RLE"){
    suppressMessages(require(edgeR))
    nf <- calcNormFactors(data,method="RLE");
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- c(msg, paste("Performed RLE Normalization"));
  }else if(norm.opt=="TMM"){
    suppressMessages(require(edgeR))
    nf <- calcNormFactors(data,method="TMM");
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- c(msg, paste("Performed TMM Normalization"));
  }else if(norm.opt=="clr"){
    data <- apply(data, 2, clr_transform);
    msg <- "Performed centered-log-ratio normalization.";
  }else if(norm.opt=='LogNorm'){
    min.val <- min(abs(data[data!=0]))/10;
    data<-apply(data, 2, LogNorm, min.val);
  }else if(norm.opt=='CrNorm'){
    norm.data <- abs(data)^(1/3);
    norm.data[data<0] <- - norm.data[data<0];
    data <- norm.data;
  }else if(norm.opt=='Rlr'){
    norm.data <- RLRNorm(data)
    msg <- paste(msg, "Performed Linear Regression Normalization.", collapse=" ");
  }else if(norm.opt=='Loess'){
    norm.data <- LoessNorm(data)
    msg <- paste(msg, "Performed Local Regression Normalization.", collapse=" ");
  }else if(norm.opt=='EigenMS'){
     msg <- paste(msg, "Performed EigenMS Normalization.", collapse=" ");
  }else if(norm.opt=='median'){
    data<- apply(data, 2, MedianNorm);
    msg <- paste(msg, "Normalization to sample median.", collapse=" ");
  }
  
  
  # scaling
  if(scaleNorm=='MeanCenter'){
    data<-apply(data, 1, MeanCenter);
    scalenm<-"Mean Centering";
  }else if(scaleNorm=='AutoNorm'){
    data<-apply(data, 1, AutoNorm);
    scalenm<-"Autoscaling";
  }else if(scaleNorm=='ParetoNorm'){
    data<-apply(data, 1, ParetoNorm);
    scalenm<-"Pareto Scaling";
  }else if(scaleNorm=='RangeNorm'){
    data<-apply(data, 1, RangeNorm);
    scalenm<-"Range Scaling";
  }else if(scaleNorm=="colsum"){
    data <- sweep(data, 2, colSums(data), FUN="/")
    data <- data*10000000;
    msg <- c(msg, paste("Performed total sum normalization."));
  }else if(scaleNorm=="upperquartile" || norm.opt == "upperquartile"){
    suppressMessages(require(edgeR))
    nf <- calcNormFactors(data,method="upperquartile");
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- c(msg, paste("Performed upper quartile normalization"));
  }else if(scaleNorm=="CSS"){
    suppressMessages(require(metagenomeSeq))
    #biom and mothur data also has to be in class(matrix only not in phyloseq:otu_table)
    data1 <- as(data,"matrix");
    dataMR <- newMRexperiment(data1);
    data <- cumNorm(dataMR,p=cumNormStat(dataMR));
    data <- MRcounts(data,norm = T);
    msg <- c(msg, paste("Performed cumulative sum scaling normalization"));
  }else{
    scalenm<-"N/A";
  }
  
  if(scaleNorm %in% c('MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm')){
    data <- t(data)
  }
  
  norm.msg <<- msg;

  rownames(data) <- row.nms;
  colnames(data) <- col.nms;

  msgSet$current.msg <- msg;
  saveSet(msgSet, "msgSet");
  return(data)
}


########
#
#Normalization internal methods
#
########

# based on phyloseq post: https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
clr_transform <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}


# generalize log, tolerant to 0 and negative values
LogNorm<-function(x, min.val){
  log10((x + sqrt(x^2 + min.val^2))/2)
}


SumNorm<-function(x){
  1000*x/sum(x, na.rm=T);
}

# normalize by median
MedianNorm<-function(x){
  x/median(x, na.rm=T);
}


# normalize to zero mean and unit variance
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

# normalize to zero mean but variance/SE
ParetoNorm<-function(x){
  (x - mean(x))/sqrt(sd(x, na.rm=T));
}

# normalize to zero mean but variance/SE
MeanCenter<-function(x){
  x - mean(x);
}

# normalize to zero mean but variance/SE
RangeNorm<-function(x){
  if(max(x) == min(x)){
    x;
  }else{
    (x - mean(x))/(max(x)-min(x));
  }
}
########### adapted from NormalyzerDE (https://github.com/ComputationalProteomics/NormalyzerDE)

RLRNorm <- function(data) {
  
  sampleLog2Median <- apply(data, 1, median,na.rm=T)
  
  calculateRLMForCol <- function(colIndex, sampleLog2Median, data) {
    
    lrFit <- MASS::rlm(as.matrix(data[, colIndex])~sampleLog2Median, na.action=stats::na.exclude)
    coeffs <- lrFit$coefficients
    coefIntercept <- coeffs[1]
    coefSlope <- coeffs[2]
    globalFittedRLRCol <- (data[, colIndex] - coefIntercept) / coefSlope
    globalFittedRLRCol
  }
  
  globalFittedRLR <- vapply(
    seq_len(ncol(data)),
    calculateRLMForCol,
    rep(0, nrow(data)),
    sampleLog2Median=sampleLog2Median,
    data=data
  )
  
  colnames(globalFittedRLR) <- colnames(data)
  
  return(globalFittedRLR)
}

LoessNorm <- function(x, weights = NULL, span=0.7, iterations = 3){
  x <- as.matrix(x)
  n <- ncol(x)
    for (k in 1:iterations) {
      a <- rowMeans(x,na.rm=TRUE)
      for (i in 1:n){
        m <- x[,i] - a
        f <- limma::loessFit(m, a, weights=weights, span=span)$fitted
        x[,i] <- x[,i] - f
      }
    }
  
  return(x)
}


# Prepare MORlog: save expression matrix into dat.in.qs
.prepare.morlog <- function(expr) {
  di <- list(expr = expr)
  qs::qsave(di, "dat.in.qs")
  return(1L)
}


# Apply MORlog back to the active dataSet
.apply.morlog <- function(dataName) {
  dataSet <- readDataset(dataName)
  di <- qs::qread("dat.in.qs")

  dataSet$expr <- di$expr
  dataSet$norm <- di$norm

  dataSet$data.norm <- dataSet$norm
  fast.write(dataSet$data.norm, file = "data_normalized.csv")
  qs::qsave(dataSet$data.norm, file = "data.stat.qs")

  msgSet <- readSet(msgSet, "msgSet")
  msgSet$current.msg <- "[MORlog] Applied DESeq2 size-factor normalization and log2(x+1)."
  saveSet(msgSet, "msgSet")

  return(RegisterData(dataSet))
}


# Microservice entrypoint: expects dat.in.qs in working dir,
# runs DESeq2 size factor normalization + log2(x+1)
morlog_micro_run <- function(expr_field = "expr", norm_field = "norm") {
  requireNamespace("DESeq2", quietly = TRUE)

  di <- qs::qread("dat.in.qs")
  m  <- as.matrix(di[[expr_field]])

  # basic checks
  if (any(m < 0, na.rm = TRUE)) stop("MORlog expects non-negative counts")
  if (!is.integer(m)) m <- round(m)

  # minimal DESeq2 object
  cd  <- S4Vectors::DataFrame(row.names = colnames(m))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = cd, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)

  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  di[[norm_field]] <- log2(norm_counts + 1)

  qs::qsave(di, "dat.in.qs")
  return(1L)
}

#' Perform batch correction using ComBat
#' @description Apply batch correction to normalized data for single dataset
#' @param dataName Name of the dataset
#' @param batchVar Name of the metadata column representing batch variable
#' @author Guangyan Zhou \email{guangyan.zhou@mail.mcgill.ca}
#' McGill University, Canada
#' License: MIT
#' @export
PerformExpressBatchCorrection <- function(dataName, batchVar) {
  .prepare.express.batch(dataName, batchVar);
  .perform.computing();
  .finalize.express.batch(dataName);
  return(1);
}

.prepare.express.batch <- function(dataName, batchVar) {
  # Read dataset before microservice
  qsfile <- gsub("\\.csv$|\\.txt$", ".qs", dataName);
  dataSet <- qs::qread(qsfile);

  my.fun <- function() {
    require('sva');

    # Read dat.in which contains the dataset
    dat.in <- qs::qread("dat.in.qs");
    dataSet <- dat.in$dataSet;
    batchVar <- dat.in$batchVar;

    # Get the normalized data
    data <- dataSet$data.norm;
    data <- as.matrix(data);
    storage.mode(data) <- "double";

    # Get metadata from dataSet
    meta.info <- dataSet$meta.info;

    # Check if batch variable exists in metadata
    if (!batchVar %in% colnames(meta.info)) {
      stop(paste("Batch variable", batchVar, "not found in metadata"));
    }

    # Align metadata rows to data columns when possible
    if (!is.null(rownames(meta.info)) && all(colnames(data) %in% rownames(meta.info))) {
      meta.info <- meta.info[colnames(data), , drop = FALSE];
    }

    # Get batch vector
    batch <- meta.info[[batchVar]];

    if (length(batch) != ncol(data)) {
      stop("Batch vector length does not match number of samples");
    }
    if (any(is.na(batch))) {
      stop("Batch variable has NA values for some samples");
    }

    # Check if batch has at least 2 levels
    if (length(unique(batch)) < 2) {
      stop("Batch variable must have at least 2 different levels");
    }

    # Create model matrix (null model, no covariates to preserve)
    mod <- model.matrix(~1, data = data.frame(sample = colnames(data)));

    # Apply ComBat
    data.batch.corrected <- ComBat(dat = data, batch = batch, mod = mod,
                                    par.prior = TRUE, prior.plots = FALSE);

    # Track max absolute delta to confirm correction changed values
    max.abs.delta <- max(abs(data.batch.corrected - data), na.rm = TRUE);

    # Save batch-corrected data back to dataSet
    dataSet$data.norm <- data.batch.corrected;

    # Save back to dat.in with message
    dat.in$dataSet <- dataSet;
    dat.in$numBatches <- length(unique(batch));
    dat.in$maxAbsDelta <- max.abs.delta;
    qs::qsave(dat.in, "dat.in.qs");
  }

  dat.in <- list(my.fun = my.fun, dataSet = dataSet, batchVar = batchVar);
  qs::qsave(dat.in, file = "dat.in.qs");
  return(1);
}

.finalize.express.batch <- function(dataName) {
  print("BATCH CORRECTION FINALIZE: Starting");

  # Read the result from microservice
  dat.in <- qs::qread("dat.in.qs");
  dataSet.corrected <- dat.in$dataSet;
  batchVar <- dat.in$batchVar;
  numBatches <- dat.in$numBatches;

  if (!is.null(dat.in$maxAbsDelta)) {
    print(paste("BATCH CORRECTION FINALIZE: max abs delta =", dat.in$maxAbsDelta));
  }

  # Save the updated dataset back to file
  qsfile <- gsub("\\.csv$|\\.txt$", ".qs", dataName);
  qs::qsave(dataSet.corrected, qsfile);
  # Update the data.anot.qs with batch-corrected data
  qs::qsave(dataSet.corrected$data.norm, file = "data.anot.qs");

  # Update the proc data as well
  qs::qsave(dataSet.corrected$data.norm, file = "data.proc.qs");

  # IMPORTANT: Update the active dataSet in the R session
  dataSet <- readDataset(dataName);

  dataSet$data.norm <- dataSet.corrected$data.norm;

  RegisterData(dataSet);

  # Update message (NOW we can use readSet/saveSet, outside microservice)
  msgSet <- readSet(msgSet, "msgSet");
  if (is.null(numBatches) || is.na(numBatches)) {
    # Fallback if microservice did not persist numBatches.
    if (!is.null(dataSet.corrected$meta.info) && batchVar %in% colnames(dataSet.corrected$meta.info)) {
      numBatches <- length(unique(dataSet.corrected$meta.info[[batchVar]]))
    } else {
      numBatches <- "unknown number of"
    }
  }
  msgSet$current.msg <- paste0("Batch correction applied using variable: ", batchVar,
                               ". Adjusted for ", numBatches, " batches.");
  saveSet(msgSet, "msgSet");
  print(paste("BATCH CORRECTION FINALIZE:", msgSet$current.msg));
  print(paste("BATCH CORRECTION FINALIZE: Done (", numBatches, " batches )"));
  return(1);
}
