##################################################
## R script for ExpressAnalyst
## Description: Functions related to web interface
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
##################################################

GetSigGeneCount <- function(){
  analSet <- readSet(analSet, "analSet");
  return(analSet$sig.gene.count);
}

GetSigGeneCountTotal <- function(){
  analSet <- readSet(analSet, "analSet");
  return(analSet$sig.gene.count.total);
}


CheckRawDataAlreadyNormalized <- function(dataName=""){
  dataSet <- readDataset(dataName);
  #data <- dataSet$data.anot;
  data <- .get.annotated.data();
  if(sum(data > 100) > 100){ # now we think it is raw counts
    return(0);
  }else{
    return(1);
  }
}

GetMetaCol<- function(dataName=""){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  if(anal.type == "onedata"){
    colNms <- colnames(dataSet$comp.res);
    if (dataSet$de.method=="limma"){
      inx <- match("AveExpr", colNms)
      return(names(dataSet$comp.res.list));
    } else if(dataSet$de.method=="wtt"){
      inx <- match("t", colNms)
  } else if (dataSet$de.method=="deseq2"){
      inx <- match("baseMean", colNms)
      return(names(dataSet$comp.res.list));
    } else {
      inx <- match("logCPM", colNms)
      return(names(dataSet$comp.res.list));
    }
    resT <- dataSet$comp.res;
    if(inx > 2){
      resT <- resT[,1:inx-1];
      nms <- gsub("logFC.", "logFC_", colnames(resT));
      
      # if there are decimals, we don't want to replace them with "vs"
      # find number of decimals in each column name
      num.dec <- lengths(regmatches(nms, gregexpr("\\.", nms)))
      
      # when only one ".", it's easy
      nms[num.dec == 1] <- gsub("\\.", " vs ", nms[num.dec == 1])
      
      # for three ".", replace only the middle
      nms[num.dec == 3] <- sapply(strsplit(nms[num.dec == 3], "\\."), function(x) {
        g <- seq_along(x)
        g[g < 2] <- 2
        g[g > 2 + 1] <- 2+1
        paste(tapply(x, g, paste, collapse = "."), collapse = " vs ")
      })
      nms <- unlist(nms)
      
      # for two ".", difficult to know which one - just leave as is
      
      return(as.vector(nms));
    }else{
      return(dataSet$par1);
    }
  }else{
    nms <- paste(unique(dataSet$cls), collapse=" vs ");
    return(nms);
  }
}

GetSummaryData <- function(){
  msgSet <- readSet(msgSet, "msgSet");
#print(msgSet$summaryVec);
  return(msgSet$summaryVec);
}

GetMetaColLength <- function(dataName = "") {

  dataSet <- readDataset(dataName)

  # bail out early if no comparison results
  if (is.null(dataSet$comp.res) && is.null(dataSet$comp.res.list)) {
    return(0L)
  }

  method <- tolower(dataSet$de.method)

  if (method == "limma") {
    return(length(dataSet$comp.res.list))
  } else if (method == "deseq2") {
    return(length(dataSet$comp.res.list))
  } else {           
    return(length(dataSet$comp.res.list))
  }

  # if the marker column isn't present → length 0
  #if (is.na(inx) || inx <= 1) {
  #  return(0L)
  #}

  #resT <- dataSet$comp.res[, seq_len(inx - 1), drop = FALSE]
  #length(colnames(resT))
}

GetInitLib <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  init.lib <- paramSet$init.lib;
  return(init.lib)
}

GetMetaDatasets<- function(){
  paramSet <- readSet(paramSet, "paramSet");
  mdata.all <- paramSet$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  return(sel.nms);
}

SetSelMetaData<- function(selNm){
    paramSet <- readSet(paramSet, "paramSet");
    paramSet$selDataNm <- selNm;
    paramSet$jsonNms$dataName <- selNm;
    saveSet(paramSet, "paramSet");
}

# only for switching single expression data results
SetCurrentData <- function(nm){
  dataSet <- readDataset(nm);
  return(1);
}

GetOmicsDataDims <- function(dataName){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  if(paramSet$anal.type == "genelist"){
  dm <- c(nrow(dataSet$prot.mat), 0);
  naNum <- 0;
  }else{
  dm <- dim(dataSet$data.norm);
  naNum <- sum(is.na(dataSet$data.norm));
  }

  return(c(dm, naNum));
} 

# given dataSet Name, sample name, and class name, do update
# note, for multiple #class, this set which one to use in the subsequent steps
# last one wins

# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# return the total matched gene number

# obtain sample names and their class labels
GetSampleInfo <- function(dataName, clsLbl){
    dataSet <- readDataset(dataName);
    grpInfo <- dataSet$meta.info[[clsLbl]];
    grpLbls <- paste(levels(grpInfo), collapse="\n");
    smplInfo <- paste(Sample = colnames(dataSet$data.orig), "\t", Class=grpInfo, collapse="\n");
    return(c(grpLbls, smplInfo));
}

#for metadata
GetMetaSummaryData<- function(){
    paramSet <- readSet(paramSet, "paramSet");
    inmex.meta <- qs::qread("inmex_meta.qs");
    sel.nms <- unique(inmex.meta$data.lbl)
    sel.nms <- paste(sel.nms, collapse="; ")
    cls.lbls <- unique(inmex.meta$cls.lbl)
    cls.lbls <- paste(cls.lbls, collapse="; ")
    paramSet$summaryVecMeta <- c(length(colnames(inmex.meta$data)),nrow(inmex.meta$data), sel.nms, cls.lbls);
    saveSet(paramSet);
    return(paramSet$summaryVecMeta)
}

GetDatasetNamesString <- function(){
    inmex.meta <- qs::qread("inmex_meta.qs");
    paste(unique(inmex.meta$data.lbl), collapse="||");
}

##Single matrix
GetSampleNumber <-function(){
  data.orig <- qs::qread("data.raw.qs");
  return(ncol(data.orig));
}


GetFilesToBeSaved <-function(naviString){
  paramSet <- readSet(paramSet, "paramSet");
  return(unique(paramSet$partialToBeSaved));
}

GetMetaInfo <- function(dataName=""){
  paramSet <- readSet(paramSet, "paramSet");
  print(paste0("metainfo==dataname=", dataName));
  if(paramSet$anal.type == "metadata"){
  metaNms<-setdiff(colnames(paramSet$dataSet$meta.info),dataSet$rmMetaCol)
  }else{
  dataSet <- readDataset(dataName);
  metaNms<-setdiff(colnames(dataSet$meta.info),dataSet$rmMetaCol)
  }
  return(metaNms);
}

GetExpressResultGeneSymbols<-function(){
  analSet <- readSet(analSet, "analSet");
  return(analSet$comp.genes.symbols);
}

GetExpressResultGeneIDLinks <- function(dataName=""){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  ids <- rownames(dataSet$comp.res);
  
  if(paramSet$data.org == "generic"){
    if(paramSet$data.idType == "ko"){
      annots <- paste0("<a href='https://www.genome.jp/dbget-bin/www_bget?", ids, "' target='_blank'>KEGG</a>");
    } else if(paramSet$data.idType == "s2f"){
      annots <- paste0("<a href='https://www.ecoomicsdb.ca/#/query?ortho=", ids, "' target='_blank'>EODB</a>");
    } else {
      annots <- ids; # Keep as-is
    }
  } else if (paramSet$data.org == "custom"){
    annots <- ids; # Keep as-is
  } else {
    annots <- paste0("<a href='http://www.ncbi.nlm.nih.gov/gene?term=", ids, "' target='_blank'>NCBI</a>");
  }
  return(annots);  # Ensure this is a character vector, NOT a list
}


GetExpressResultColNames<-function(){
  resT <- qs::qread("express.de.res.qs");
  colnames(resT);
}

GetExpressResultGeneIDs<-function(dataName=""){
    dataSet <- readDataset(dataName);
    return(rownames(dataSet$comp.res));
}

GetExpressGeneIDType<-function(dataName=""){
  dataSet <- readDataset(dataName);
  return(dataSet$id.current);
}

GetExpressResultMatrix <- function(dataName = "", inxt) {
    dataSet  <- readDataset(dataName);
    paramSet <- readSet(paramSet, "paramSet");
    inxt     <- as.numeric(inxt)

    if (dataSet$de.method == "deseq2") {

        inx <- match("baseMean", colnames(dataSet$comp.res))
        res <- dataSet$comp.res.list[[inxt]];
    }else if (dataSet$de.method=="edger"){
        inx <- match("logCPM", colnames(dataSet$comp.res))
        res <- dataSet$comp.res.list[[inxt]];
    }else if (dataSet$de.method=="limma"){
        inx <- match("AveExpr", colnames(dataSet$comp.res))
        res <- dataSet$comp.res.list[[inxt]];
    } else{
         if (dataSet$de.method == "wtt") {
            inx <- match("t", colnames(dataSet$comp.res))
        } 

        res <- dataSet$comp.res
        res <- res[, -(1:(inx - 1)), drop = FALSE]                    # ← fixed slice
        res <- res[rownames(dataSet$comp.res), , drop = FALSE]        # ← align rows
        res <- cbind(dataSet$comp.res[, inxt], res)
        colnames(res)[1] <- colnames(dataSet$comp.res)[inxt]
    }

    dataSet$comp.res <- dataSet$comp.res[order(dataSet$comp.res$adj.P.Val), ]
    dataSet$comp.res <- dataSet$comp.res[
        !(rownames(dataSet$comp.res) %in% rownames(dataSet$sig.mat)), ]
    dataSet$comp.res <- rbind(dataSet$sig.mat, dataSet$comp.res)
    dataSet$comp.res <- dataSet$comp.res[complete.cases(dataSet$comp.res), ]

    ## --- now extract the column(s) for the return value -------
    if (dataSet$de.method %in% c("limma", "edger", "deseq2")) {
      res <- dataSet$comp.res.list[[inxt]]
    } else {
      res <- dataSet$comp.res[ , c(inxt, (inx+1):ncol(dataSet$comp.res)), drop = FALSE]
      res <- res[order(res$adj.P.Val), ]
      colnames(res)[1] <- colnames(dataSet$comp.res)[inxt]
    }

    RegisterData(dataSet)
    qs::qsave(res, "express.de.res.qs")

    return(head(signif(as.matrix(res), 5),1000))
}


###Gene list
GetNumOfLists <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$numOfLists)
}

GetNumOfGenes <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$genenum)
}


GetNumOfAnnoGenes <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$annonum)
}


GetMetaSigGeneCount <- function(){
  analSet <- readSet(analSet, "analSet");
  return(nrow(analSet$meta.mat));
}

GetCurrentJson <- function(type) {
  paramSet <- readSet(paramSet, "paramSet")
  
  # Check if the list paramSet$jsonNms contains the key 'type'
  if (!is.null(paramSet$jsonNms[[type]])) {
    return(paramSet$jsonNms[[type]])
  } else {
    return("NA")  # or return a default value if appropriate
  }
}


SelectDataSet <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  if(!exists('nm.vec')){
    AddErrMsg("No dataset is selected for analysis!");
    return(0);
  }
  mdata.all <- paramSet$mdata.all
  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <- 1;
    }else{
      mdata.all[[nm]] <- 0;
    }
  }
  
  if("meta_dat" %in% nm.vec){
    meta.selected <<- TRUE;
  }else{
    meta.selected <<- FALSE;
  }
  
  rm('nm.vec', envir = .GlobalEnv);

  paramSet$mdata.all <- mdata.all
  saveSet(paramSet, "paramSet");
  return(1);
  
}


GetFeatureNum <- function(dataName) {

  dataSet <- readDataset(dataName, quiet = TRUE)

  if (is.null(dataSet) || is.null(dataSet$data.norm)) {
    return(0L)                               # nothing loaded → report zero
  }

  nrow(dataSet$data.norm)
}

# get qualified inx with at least number of replicates
GetDiscreteInx <- function(my.dat, min.rep=2){
  good.inx <- apply(my.dat, 2, function(x){
    x <- x[x!="NA"]
    good1.inx <- length(x) > length(unique(x));
    good2.inx <- min(table(x)) >= min.rep;
    return (good1.inx & good2.inx);
  });
  return(good.inx);
}

GetNumbericalInx <- function(my.dat){
  suppressWarnings({
  good.inx <- apply(my.dat, 2, function(x){
    isNum = as.numeric(as.character(x[x!="NA"]))
    return(all(!is.na(as.numeric(as.character(isNum)))));
  });
  })
  return(good.inx);
}

.set.dataSet <- function(dataSetObj=NA){
  RegisterData(dataSetObj);
  return (1);
}

# remove data object, the current dataSet will be the last one by default 
RemoveData <- function(dataName){
  paramSet <- readSet(paramSet, "paramSet");
  mdata.all <- paramSet$mdata.all;
  if(!is.null(paramSet$mdata.all[[dataName]])){
    paramSet$mdata.all[[dataName]] <- NULL;
  }
  saveSet(paramSet, "paramSet");
}


GetCovSigFileName <-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.nm;
}

GetCovSigMat<-function(dataName){
  dataSet <- readDataset(dataName);
  drops <- c("ids","label")
  return(CleanNumber(as.matrix(dataSet$analSet$cov$sig.mat[, !(names(dataSet$analSet$cov$sig.mat) %in% drops)])));
}

GetCovSigIds<-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.mat$ids;
}

GetCovSigSymbols<-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.mat$label
}

GetCovSigColNames<-function(dataName){
  dataSet <- readDataset(dataName);
  drops <- c("ids","label");
  colnames(dataSet$analSet$cov$sig.mat[,!(names(dataSet$analSet$cov$sig.mat) %in% drops)]);
}

GetCovDENums <- function(dataName){
    deNum <- nrow(dataSet$analSet$cov$sig.mat);
    nonDeNum <- nrow(dataSet$comp.res) - deNum;
    return(c(deNum, nonDeNum));
}


#'Replace infinite numbers
#'@description Replace -Inf, Inf to 99999 and -99999
#'@param bdata Input matrix to clean numbers
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'
CleanNumber <-function(bdata){
  if(sum(bdata==Inf)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- 999999;
  }
  if(sum(bdata==-Inf)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- -999999;
  }
  bdata;
}

GetMetaMethodPVal <-function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$BHth);
}

#for enrichment analysis
SetUniverseOpt <- function(universe.opt){
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$universe.opt <- universe.opt;
  if(paramSet$universe.opt == "uploaded"){
    paramSet$universe.opt.readable <- "Uploaded Data";
  }else{
    paramSet$universe.opt.readable <- "Gene Set Library";
  }
  saveSet(paramSet);
}

SetInitLib <-function(library){
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$init.lib <- library;
  saveSet(paramSet);
}
