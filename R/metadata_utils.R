##################################################
## R script for ExpressAnalyst
## Description: functions for complex metadata
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#####'Sanity check metadata after metadata edited 
SanityCheckMeta <- function(fileName, init){
  #save.image("sanity.RData")          # debug snapshot
  msgSet   <- readSet(msgSet,   "msgSet")
  paramSet <- readSet(paramSet, "paramSet")

  sel.nms <- if (fileName == "NA") names(paramSet$mdata.all) else fileName

  for (nm in sel.nms) {
    dataSet <- readDataset(nm)
    meta    <- dataSet$meta.info            # may be NULL or 0-col

    ## ------------------------------------------------------------------
    ## 1. Guarantee at least one metadata column (“CLASS”) exists
    ## ------------------------------------------------------------------
    if (is.null(meta) || ncol(meta) == 0) {
      meta <- data.frame(
        CLASS   = factor(dataSet$cls),
        row.names = colnames(dataSet$data.norm),
        stringsAsFactors = FALSE
      )
      dataSet$meta.info <- meta
      # also update index vectors so downstream code stays happy
      dataSet$disc.inx  <- setNames(TRUE,  "CLASS")
      dataSet$cont.inx  <- setNames(FALSE, "CLASS")
    }

    ## ------------------------------------------------------------------
    ## 2. Optional QC on missing values (your original logic)
    ## ------------------------------------------------------------------
    if (init != 1 && (any(is.na(meta)) | any(meta == "") | any(meta == "NA"))) {
      return(2)
    }

    ## ------------------------------------------------------------------
    ## 3. Use the first metadata column as the class label
    ## ------------------------------------------------------------------
    cls      <- meta[meta[, 1] != "NA", 1]
    cls.lbl  <- factor(as.character(cls))

    # replicate check
    if (min(table(cls.lbl)) < 2) {
      msgSet$current.msg <- paste0(
        "No replicates detected for group ",
        names(which(table(cls.lbl) < 2)), " in ", colnames(meta)[1]
      )
      saveSet(msgSet, "msgSet")
      return(0)
    }

    ## 4. ensure all metadata columns are factors
    meta[] <- lapply(meta, as.factor)

    ## 5. stash back and persist
    dataSet$cls    <- cls.lbl
    dataSet$rmidx  <- which(meta[, 1] == "NA")
    dataSet$meta.info <- meta
    RegisterData(dataSet)
  }

  return(1)
}


# here should first try to load the original data
# the data in the memory could be changed
GetGroupNames <- function(dataName, meta="NA"){
  dataSet <- readDataset(dataName);
  if(meta == "NA"){
    grpnms = levels(factor(dataSet$meta.info[,1]));
  }else{
    grpnms =levels(factor(dataSet$meta.info[,meta]));
  }
  return(grpnms[grpnms!="NA"])
}

GetResRowNames <- function(dataName=""){
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
    return(rownames(dataSet$meta.info));
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    dataSet <- paramSet$dataSet;
    return(rownames(paramSet$dataSet$meta.info));
  }
}

GetResColNames <- function(dataName=""){
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
    colnms<- colnames(dataSet$meta.info)[colnames(dataSet$meta.info)!="newcolumn"]
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    colnms<- colnames(paramSet$dataSet$meta.info)[colnames(paramSet$dataSet$meta.info)!="newcolumn"]
  }
  return(colnms);
}

filterDatasetColumn <- function(colNames) {
  colNames[!toupper(colNames) %in% "DATASET"]
}

GetDiscMetas <- function(dataName=""){
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    dataSet <- paramSet$dataSet;
  }

  # Handle missing metadata
  if(is.null(dataSet$meta.info) || is.null(dataSet$disc.inx)){
    return(character(0));
  }

  if(length(keepVec)>0){
    keepidx <- which(keepVec %in% colnames(dataSet$meta.info))
    keepidx <- intersect(keepidx,which(dataSet$disc.inx))
  }else{
    keepidx <-  which(dataSet$disc.inx)
  }
  colnms<- colnames(dataSet$meta.info)[keepidx]
  #print(colnms)

  # Return empty character vector if no discrete metadata columns found
  if(is.null(colnms) || length(colnms) == 0){
    return(character(0));
  }

  return(filterDatasetColumn(colnms));
}

GetMetaDataCol <- function(dataName="",colnm){
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    dataSet <- paramSet$dataSet;
  }

  # Handle missing metadata or column
  if(is.null(dataSet$meta.info)){
    return(character(0));
  }

  if(!(colnm %in% colnames(dataSet$meta.info))){
    return(character(0));
  }

  cls = levels(dataSet$meta.info[,colnm]);
  if(is.null(cls) || length(cls) == 0){
    return(character(0));
  }
  return(cls[cls!="NA"]);
}

GetMetaCell <- function(dataName="",ridx=1,cidx=1){
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    dataSet <- paramSet$dataSet;
  }
  return(dataSet$meta.info[ridx,cidx]);
}

# Note R is column as a vector, operate on row 
# will lead to different factors, need to transpose
GetMetaRow <- function(dataName, ridx=1){
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    dataSet <- paramSet$dataSet;
  }
  my.meta.info <- t(dataSet$meta.info);
  return(as.character(my.meta.info[, ridx])); # now column operation
}

ResetMetaTab <- function(dataName=""){
  paramSet <- readSet(paramSet, "paramSet")
  paramSet$excluded.samples <- character(0)
  saveSet(paramSet, "paramSet")
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
     if(dataSet$type=="prot"){
       data.anot <- ov_qs_read("data.missed.qs");
    }else{
       data.anot <- ov_qs_read("orig.data.anot.qs");
    }
    dataSet$data.norm <- data.anot;
    .save.annotated.data(data.anot);
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    sel.nms <- names(paramSet$mdata.all);
    for(i in 1:length(sel.nms)){
      dataSet <- readDataset(sel.nms[i]);
      dataSet$data.norm <- dataSet$data.anot <- readDataQs("data.annotated.qs", paramSet$anal.type, sel.nms[i]);
      RegisterData(dataSet);
    }
    dataSet <- paramSet$dataSet;
  }
  dataSet$meta.info <- dataSet$metaOrig;
  dataSet$disc.inx <- dataSet$disc.inx.orig;
  dataSet$cont.inx <- dataSet$cont.inx.orig;
  RegisterData(dataSet);
}

GetResColType <- function(dataName="",colNm="NA"){
  if(dataName != "NA"){
    dat <- readDataset(dataName);
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    dat <- paramSet$dataSet;
  }

  # readDataset returns NULL when the .qs file is missing (stale dataName
  # held by a Java caller after the R session was reinitialised). Without
  # this guard, the rep(T, ncol(NULL)) below errors with
  # "invalid 'times' argument".
  if(is.null(dat) || is.null(dat$meta.info) || length(dim(dat$meta.info)) < 2){
    return(character(0));
  }

  if(colNm=="NA"){
    meta.status <- ifelse(dat$disc.inx,"disc","cont")
  }else{
    meta.status <- ifelse(dat$disc.inx[colNm],"disc","cont")
  }

  if(length(meta.status) == 0){
    meta.status <- rep(T, ncol(dat$meta.info))
    names(meta.status) <- colnames(dat$meta.info);
    dat$disc.inx <- meta.status;
    if(dataName != "NA"){
      RegisterData(dat);
    }else{
      paramSet$dataSet <- dat;
      saveSet(paramSet, "paramSet");
    }
  }
  return(meta.status);
}

UpdateMetaType <-function(dataName="", metadata="NA", type="disc"){
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
    dataSet$meta.types[metadata] = type;
  }else{
    paramSet <- readSet(paramSet, "paramSet"); 
    paramSet$dataSet$meta.types[metadata] = type;
  }
  return(.set.rdt.set(rdtSet));
}

GetUniqueMetaNames <-function(dataName="",metadata){
  paramSet <- readSet(paramSet, "paramSet"); 
  data.type <- paramSet$dataSet[["meta.types"]][metadata];
  if(data.type == "cont"){
    return("--- NA ---");
  } else {
    return(levels(as.factor(paramSet$dataSet$meta.info[,metadata])));
  }
}

UpdateMetaStatus <- function(dataName="",colNm){
  dataSet <- readDataset(dataName);
  msgSet <- readSet(msgSet, "msgSet");
  cidx <- which(colnames(dataSet$meta.info)==colNm)
  old = ifelse(dataSet$disc.inx[cidx],"Discrete","Continuous")
  if(dataSet$disc.inx[cidx]){
    if(all(is.na( as.numeric(as.character(dataSet$meta.info[,cidx]))))){
      msgSet$current.msg <- "Category metadata cannot be continuous data!"
      saveSet(msgSet, "msgSet"); 
      return(0)
    }
    dataSet$disc.inx[cidx]=FALSE;
    dataSet$cont.inx[cidx]=TRUE;
  }else{
    if(all(!duplicated(as.character(dataSet$meta.info[,cidx])))){
      msgSet$current.msg <- "No duplicates were detected! The metadata cannot be discrete!"
      saveSet(msgSet, "msgSet"); 
      return(0)
    }
    dataSet$disc.inx[cidx]=TRUE;
    dataSet$cont.inx[cidx]=FALSE;
  }
  new = ifelse(dataSet$disc.inx[cidx],"Discrete","Continuous")
  msgSet$current.msg <- paste0("Metadata type of ",colnames(dataSet$meta.info)[cidx]," has been changed to ", new, " !")
  saveSet(msgSet, "msgSet"); 
  RegisterData(dataSet);
  paramSet <- readSet(paramSet, "paramSet")
  if (is.null(paramSet$excluded.samples)) {
    paramSet$excluded.samples <- character(0)
  }
  paramSet$excluded.samples <- unique(c(paramSet$excluded.samples, samplNm))
  saveSet(paramSet, "paramSet")
  paramSet <- readSet(paramSet, "paramSet")
  if (is.null(paramSet$excluded.samples)) {
    paramSet$excluded.samples <- character(0)
  }
  paramSet$excluded.samples <- unique(c(paramSet$excluded.samples, samplNm))
  saveSet(paramSet, "paramSet")
  return(1);
}


DeleteSample <- function(dataName="",samplNm){
  #print(dataName)
  #print(samplNm)
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
    dataSet$meta.info <- dataSet$meta.info[rownames(dataSet$meta.info)!=samplNm,,drop=F]
    dataSet$data.norm <- dataSet$data.norm[,colnames(dataSet$data.norm)!=samplNm]
    RegisterData(dataSet);
  }else{
    paramSet <- readSet(paramSet, "paramSet")  
    dataName <- paramSet$dataSet$meta.info$Dataset[rownames(paramSet$dataSet$meta.info)==samplNm];
    paramSet$dataSet$meta.info <- paramSet$dataSet$meta.info[rownames(paramSet$dataSet$meta.info)!=samplNm,,drop=F];
    
    dataSet <- readDataset(dataName);
    dataSet$data.norm <- dataSet$data.norm[,colnames(dataSet$data.norm)!=samplNm];
    dataSet$meta.info <- dataSet$meta.info[rownames(dataSet$meta.info)!=samplNm,,drop=F];
    
    inmex.meta<-ov_qs_read("inmex_meta.qs");
    inmex.meta$data <- inmex.meta$data[,colnames(inmex.meta$data) !=samplNm]
    ov_qs_save(inmex.meta, "inmex_meta.qs");
    saveSet(paramSet, "paramSet");
    RegisterData(dataSet);
  }
  
  return(1);
}

GetUploadedSampleNames <- function(dataName="", type="kept"){
  if(dataName == "NA"){
    paramSet <- readSet(paramSet, "paramSet");
    if(type == "all" && !is.null(paramSet$dataSet$meta.info.full)){
      return(rownames(paramSet$dataSet$meta.info.full));
    }
    return(rownames(paramSet$dataSet$meta.info));
  }
  dataSet <- readDataset(dataName);
  if(type == "all" && !is.null(dataSet$meta.info.full)){
    return(rownames(dataSet$meta.info.full));
  }
  return(rownames(dataSet$meta.info));
}

#'Update the sample set for one dataset
#'@description Keeps only the samples named in smpl.nm.vec. Surviving rows keep
#'their (possibly edited) metadata; rows added back are restored from the copy
#'taken on first use, together with their data columns, so an exclusion is always
#'reversible. Applies to the expression matrix and the sample metadata together,
#'which is what keeps the two in step.
#'@export
UpdateUploadedSampleItems <- function(dataName=""){
  if(!exists("smpl.nm.vec")){
    current.msg <<- "Cannot find the sample names to keep!";
    return(0);
  }
  if(dataName == "NA"){
    return(.UpdateMetaSampleSet(smpl.nm.vec));
  }
  dataSet <- readDataset(dataName);

  if(is.null(dataSet$meta.info.full)){
    dataSet$meta.info.full <- dataSet$meta.info;
    dataSet$data.norm.full <- dataSet$data.norm;
  }

  full.meta <- dataSet$meta.info.full;
  keep <- intersect(rownames(full.meta), smpl.nm.vec);
  if(length(keep) < 3){
    current.msg <<- "At least three samples must remain for analysis.";
    return(0);
  }

  cur <- dataSet$meta.info;
  kept.cur <- cur[intersect(rownames(cur), keep), , drop=FALSE];
  add.back <- setdiff(keep, rownames(cur));
  if(length(add.back) > 0){
    common.cols <- intersect(colnames(cur), colnames(full.meta));
    kept.cur <- rbind(kept.cur[, common.cols, drop=FALSE],
                      full.meta[add.back, common.cols, drop=FALSE]);
  }
  ord <- rownames(full.meta)[rownames(full.meta) %in% rownames(kept.cur)];
  dataSet$meta.info <- kept.cur[ord, , drop=FALSE];

  full.data <- dataSet$data.norm.full;
  dataSet$data.norm <- full.data[, colnames(full.data) %in% ord, drop=FALSE];

  excluded <- setdiff(rownames(full.meta), keep);
  current.msg <<- paste0(length(ord), " samples kept; ", length(excluded), " excluded.");
  RegisterData(dataSet);
  return(1);
}

# Meta-analysis has no single dataset to edit: the sample set spans every uploaded
# study plus the merged matrix. Keep all four representations in step -- the combined
# metadata on paramSet, each study's own meta.info + data.norm, and the merged
# inmex_meta matrix, whose cls.lbl and data.lbl are per-COLUMN vectors and so have to
# be subset alongside the data (DeleteSample subset only the data, which left the two
# labels one element longer than the matrix after every removal).
.UpdateMetaSampleSet <- function(keep.vec){
  paramSet <- readSet(paramSet, "paramSet");
  if(is.null(paramSet$dataSet$meta.info)){
    current.msg <<- "No combined metadata is available for editing.";
    return(0);
  }
  if(is.null(paramSet$dataSet$meta.info.full)){
    paramSet$dataSet$meta.info.full <- paramSet$dataSet$meta.info;
  }
  full.meta <- paramSet$dataSet$meta.info.full;
  keep <- intersect(rownames(full.meta), keep.vec);
  if(length(keep) < 3){
    current.msg <<- "At least three samples must remain for analysis.";
    return(0);
  }
  ord <- rownames(full.meta)[rownames(full.meta) %in% keep];

  cur <- paramSet$dataSet$meta.info;
  kept.cur <- cur[intersect(rownames(cur), keep), , drop=FALSE];
  add.back <- setdiff(keep, rownames(cur));
  if(length(add.back) > 0){
    common.cols <- intersect(colnames(cur), colnames(full.meta));
    kept.cur <- rbind(kept.cur[, common.cols, drop=FALSE],
                      full.meta[add.back, common.cols, drop=FALSE]);
  }
  paramSet$dataSet$meta.info <- kept.cur[ord, , drop=FALSE];

  for(dataName in unique(as.character(full.meta$Dataset))){
    dataSet <- readDataset(dataName);
    if(is.null(dataSet)){
      next;
    }
    if(is.null(dataSet$meta.info.full)){
      dataSet$meta.info.full <- dataSet$meta.info;
      dataSet$data.norm.full <- dataSet$data.norm;
    }
    d.full.meta <- dataSet$meta.info.full;
    d.keep <- rownames(d.full.meta)[rownames(d.full.meta) %in% ord];
    if(length(d.keep) == 0){
      next;
    }
    d.cur <- dataSet$meta.info;
    d.kept <- d.cur[intersect(rownames(d.cur), d.keep), , drop=FALSE];
    d.add <- setdiff(d.keep, rownames(d.cur));
    if(length(d.add) > 0){
      d.cols <- intersect(colnames(d.cur), colnames(d.full.meta));
      d.kept <- rbind(d.kept[, d.cols, drop=FALSE],
                      d.full.meta[d.add, d.cols, drop=FALSE]);
    }
    dataSet$meta.info <- d.kept[d.keep, , drop=FALSE];
    d.full.data <- dataSet$data.norm.full;
    dataSet$data.norm <- d.full.data[, colnames(d.full.data) %in% d.keep, drop=FALSE];
    RegisterData(dataSet);
  }

  .SubsetMergedMetaMatrix("inmex_meta.qs", "inmex_meta_full.qs", ord);
  .SubsetMergedMetaMatrix("inmex.meta.orig.qs", "inmex.meta.orig_full.qs", ord);

  saveSet(paramSet, "paramSet");
  current.msg <<- paste0(length(ord), " samples kept; ",
                         length(setdiff(rownames(full.meta), keep)), " excluded.");
  return(1);
}

# Survivors keep whatever the live matrix holds (it may carry a batch correction the
# original copy does not); only added-back columns come from the pristine snapshot.
.SubsetMergedMetaMatrix <- function(live.file, full.file, keep.nms){
  if(!file.exists(live.file)){
    return(invisible(FALSE));
  }
  if(!file.exists(full.file)){
    ov_qs_save(ov_qs_read(live.file), full.file);
  }
  live <- ov_qs_read(live.file);
  full <- ov_qs_read(full.file);
  if(is.null(live$data) || is.null(full$data)){
    return(invisible(FALSE));
  }
  ord <- colnames(full$data)[colnames(full$data) %in% keep.nms];
  add.back <- setdiff(ord, colnames(live$data));
  mat <- live$data[, ord[ord %in% colnames(live$data)], drop=FALSE];
  if(length(add.back) > 0){
    feats <- intersect(rownames(mat), rownames(full$data));
    mat <- cbind(mat[feats, , drop=FALSE], full$data[feats, add.back, drop=FALSE]);
  }
  mat <- mat[, ord[ord %in% colnames(mat)], drop=FALSE];
  # Index the ORIGINAL vectors positionally rather than rebuilding them from
  # characters: cls.lbl is a factor, and re-factor()ing re-derives its levels
  # alphabetically, which silently changes the reference group and flips the sign
  # of every logFC downstream. Subsetting keeps the level order; droplevels only
  # removes a group the user has excluded entirely.
  idx <- match(colnames(mat), colnames(full$data));
  live$data <- mat;
  if(!is.null(full$cls.lbl)){
    live$cls.lbl <- droplevels(full$cls.lbl[idx]);
  }
  if(!is.null(full$data.lbl)){
    live$data.lbl <- full$data.lbl[idx];
  }
  ov_qs_save(live, live.file);
  invisible(TRUE);
}

DeleteMetaCol <- function(dataName="",metaCol){
  if(dataName != "NA"){
    sel.nms <- c(dataName);
  }else{
    paramSet <- readSet(paramSet, "paramSet")  
    sel.nms <- names(paramSet$mdata.all);
  }
  
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i]);
    idx = which(colnames(dataSet$meta.info)==metaCol)
    dataSet$meta.info <- dataSet$meta.info[,-idx,drop=F]
    dataSet$disc.inx <- dataSet$disc.inx[-idx]
    dataSet$cont.inx <- dataSet$cont.inx[-idx]
    if(!exists("rmMetaCol",dataSet)){
      dataSet$rmMetaCol <- vector()
    }
    dataSet$rmMetaCol <- unique(c(dataSet$rmMetaCol,metaCol))
    RegisterData(dataSet);
  }
  
  #for meta-anal also remove from meta.info
  if(dataName == "NA"){
    paramSet <- readSet(paramSet, "paramSet");
    idx = which(colnames(paramSet$dataSet$meta.info)==metaCol);
    paramSet$dataSet$meta.info <- paramSet$dataSet$meta.info[,-idx,drop=F];
    paramSet$dataSet$disc.inx <- paramSet$dataSet$disc.inx[-idx]
    paramSet$dataSet$cont.inx <- paramSet$dataSet$cont.inx[-idx]
  } 
  
  return(1);
}

CleanRmCol <- function(dataName=""){
  if(dataName != "NA"){
    paramSet <- readSet(paramSet, "paramSet")
    mdata.all <- paramSet$mdata.all
    sel.nms <- names(mdata.all);
  }else{
    sel.nms <- c(dataName);
  }
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i]);
    if(exists("rmMetaCol",dataSet)){
      dataSet$rmMetaCol <- vector()
    }
    RegisterData(dataSet);
  }
  return(1);
}

GetSampleNm <- function(dataName="",ridx=1){
  if(dataName != "NA"){
    dataSet <- readDataset(dataName);
  }else{
    paramSet <- readSet(paramSet, "paramSet")
    dataSet <- paramSet$dataSet;
  }
  return( rownames(dataSet$meta.info)[ridx]);
}


UpdateSampInfo <-  function(dataName="",rowNm,colNm,cell){
  if(dataName == "NA"){
    paramSet <- readSet(paramSet, "paramSet"); 
    dataName <- paramSet$dataSet$meta.info$Dataset[rownames(paramSet$dataSet$meta.info)==rowNm];
    if(colNm==""){
      rownames(paramSet$dataSet$meta.info)[rownames(paramSet$dataSet$meta.info)==rowNm] <- cell;
      saveSet(paramSet, "paramSet");
    }
  }
  
  dataSet <- readDataset(dataName);
  meta <- dataSet$meta.info
  ridx <- which(rownames(meta)==rowNm)
  if(colNm==""){
    if(rowNm != cell){
      rownames(meta)[ridx]=cell
      colnames(dataSet$data.norm)[which(colnames(dataSet$data.norm)==rowNm)] <- cell 
      #if(exists("dataSet$data.anot")){
       if(file.exists("data.anot.qs")){
        data.anot <- .get.annotated.data();
        colnames(data.anot)[which(colnames(data.anot)==rowNm)] <- cell;
        .save.annotated.data(data.anot);
      }
    }
  }else{  
    cidx<- which(colnames(meta)==colNm)
    if(cell!= as.character(meta[ridx,cidx])){
      if(cell %in% levels(meta[,cidx])){
        meta[ridx,cidx] = cell
      }else{
        levels(meta[,cidx]) <- c(levels(meta[,cidx]), cell)
        meta[ridx,cidx] = cell
      }
      meta[,cidx] <- droplevels(meta[,cidx])
    }
  }
  dataSet$meta.info = meta
  RegisterData(dataSet);
  return(1);
}


GetSelectedMetaInfo <- function(dataName="",colNm){
  dataSet <- readDataset(dataName);
  lvls <- levels(dataSet$meta.info[,colNm])
  lvls <-  lvls[lvls!="NA"]
  return(lvls);
}

UpdateMetaOrder <- function(dataName="",metacol){
  dataSet <- readDataset(dataName);
  if(length(metaVec)>0 & metacol %in% colnames(dataSet$meta.info)){
    dataSet$meta.info[,metacol] <- factor(as.character(dataSet$meta.info[,metacol]),levels = metaVec)
    
  }else{
    msgSet <- readSet(msgSet, "msgSet");
    msgSet$current.msg <- "The metadata column is empty! Please check your selection!"
    saveSet(msgSet, "msgSet"); 
    return(0)
  }
  RegisterData(dataSet);
  return(1)
}

UpdateMetaName <-  function(dataName="",oldvec,newvec){
  if(dataName == "NA"){
    sel.nms <- names(paramSet$mdata.all);
    paramSet <- readSet(paramSet, "paramSet");
    idx <- which(colnames(paramSet$dataSet$meta.info)==oldvec)
    if(length(idx)==1){
      colnames(paramSet$dataSet$meta.info)[idx] <- names(paramSet$dataSet$disc.inx)[idx] <- 
        names(paramSet$dataSet$cont.inx)[idx] <- newvec
    }else{
      return(0)
    }
    saveSet(paramSet, "paramSet");
  }else{
    sel.nms <- c(dataName);
  }
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i]);
    idx <- which(colnames(dataSet$meta.info)==oldvec)
    if(length(idx)==1){
      colnames(dataSet$meta.info)[idx] <- names(dataSet$disc.inx)[idx] <- 
        names(dataSet$cont.inx)[idx] <- newvec
    }else{
      return(0)
    }
    RegisterData(dataSet);
  }
  return(1);
}

GetMetaSummary <- function(dataName = "") {
    paramSet <- readSet(paramSet, "paramSet");

  dataSet <- readDataset(dataName)
  meta    <- dataSet$meta.info

  ## -- discrete -----------------------------------------------------------
  disc.len <- if (!is.null(dataSet$disc.inx))
                length(which(dataSet$disc.inx))
              else 0L
  
  disc.vec <- if (disc.len > 0)
                paste(names(dataSet$disc.inx)[which(dataSet$disc.inx)],
                      collapse = ", ")
              else ""
  
  ## -- continuous  (NEW: gracefully handle NULL) --------------------------
  cont.len <- if (!is.null(dataSet$cont.inx))
                length(which(dataSet$cont.inx))
              else 0L
  
  cont.vec <- if (cont.len > 0)
                paste(names(dataSet$cont.inx)[which(dataSet$cont.inx)],
                      collapse = ", ")
              else ""
  
  ## -- NA check -----------------------------------------------------------
  na.vec <- na.check(meta)
  
  ## -- build result vector -------------------------------------------------
  # The caller indexes all nine positions, so the vector must always have nine. c() DROPS a
  # zero-length element, and several of these can be NULL when meta.info is not a proper data
  # frame: ncol() on a non-matrix is NULL, and names(x)[1] on an unnamed object is NULL too.
  # Two silent drops shifted every later field left and the caller read past the end.
  .one <- function(x, default = "NA") {
    if (is.null(x) || length(x) == 0L) return(default)
    x <- as.character(x)[1]
    if (is.na(x)) default else x
  }
  .col1 <- tryCatch(meta[, 1], error = function(e) NULL)
  res <- c(
    .one(if (is.null(dim(meta))) length(meta) else ncol(meta), "0"),  # total metadata columns
    .one(disc.len, "0"),   # # discrete cols
    .one(disc.vec, ""),    # names of discrete cols
    .one(cont.len, "0"),   # # continuous cols   (== 0 if cont.inx NULL)
    .one(cont.vec, ""),    # names of continuous cols
    .one(names(meta)[1], "NA"),                                       # first column name
    .one(length(unique(.col1)), "0"),
    .one(paste(unique(.col1), collapse = ", "), ""),
    .one(na.vec, "None")
  )
  
  ## -- save & return -------------------------------------------------------
  paramSet$metadata.summary <- res
  saveSet(paramSet, "paramSet")

  return(res)
}

na.check <- function(mydata){
  # Guard against NULL / scalar / zero-column input. apply() blows up on
  # those with "dim(X) must have a positive length", which is what we see
  # when GetMetaSummary() is invoked with a stale dataName whose qs file
  # is missing (readDataset returns NULL, so dataSet$meta.info is NULL).
  if(is.null(mydata) || length(dim(mydata)) < 2 || ncol(mydata) == 0){
    return("None")
  }
  na.idx <- apply(mydata,2,function(x) "NA" %in% x)
  if(all(!na.idx)){
    return("None")
  }
  na.num <- apply(mydata,2,function(x) length(which(x=="NA")))
  naInfo <- data.frame(names(mydata)[na.idx],num = na.num[na.num>0])
  naInfo <- apply(naInfo, 1, function(x) paste0(x[1]," (",x[2],")"))
  naInfo <- paste(naInfo,collapse = ", ")
  return(naInfo)
}


UpdatePrimaryMeta <- function(fileName,primaryMeta){
  dataSet <- readDataset(fileName);
  msgSet <- readSet(msgSet,"msgSet");
  meta <- dataSet$meta.info
  if(primaryMeta %in% colnames(meta)){
    cidx <- which(colnames(meta)==primaryMeta)
    dataSet$meta.info<-cbind(meta[,cidx,drop=F],meta[,-cidx,drop=F])
    dataSet$disc.inx=c(dataSet$disc.inx[cidx],dataSet$disc.inx[-cidx])
    dataSet$cont.inx=c(dataSet$cont.inx[cidx],dataSet$cont.inx[-cidx])
  }else{
    msgSet$current.msg <- "The metadata column is empty! Please check your selection!"
    saveSet(msgSet, "msgSet"); 
    return(0)
  }
  RegisterData(dataSet);
  return(1)
}


GetMetaDataGroups <- function(dataName){
  paramSet <- readSet(paramSet, "paramSet");
  groups <- colnames(paramSet$dataSet$meta.info);
  return(filterDatasetColumn(groups));
}

GetMetaDataStatus <- function(dataName){
  paramSet <- readSet(paramSet, "paramSet");
  res <- unname(paramSet$dataSet$meta.status);
  return(res);
}

# Per-dataset processing state, read from the persisted R dataSet object.
# Java's getDataSets() rebuilds OmicsModel objects from R when its in-memory list
# is empty (session-bean recreation, project restore, example switch); those
# rebuilt objects default every processing flag to FALSE, so the meta Integrity
# Check saw allDone=false and blocked with "specify or confirm group comparison".
# The R dataSet is the source of truth — expose it so the reconstruction can
# restore annotated / normalized / compared instead of guessing.
# Returns c(annotated, normalized, compared) as 0/1 integers.
GetDataProcStatus <- function(dataName){
  ds <- readDataset(dataName);
  if (is.null(ds)) return(as.integer(c(0, 0, 0)));
  # The boolean ds$annotated flag is unreliable (not always set), so key off the
  # annotated matrix / resolved id type. Normalization is marked by norm.opt being
  # set (data.norm alone exists straight after read, so it can't distinguish).
  annotated  <- !is.null(ds$data.annotated) || (!is.null(ds$id.type) && nzchar(ds$id.type));
  normalized <- !is.null(ds$norm.opt) || (!is.null(ds$data.filtered) && !is.null(ds$data.norm));
  compared   <- !is.null(ds$comp.res) || (!is.null(ds$analSet) && !is.null(ds$analSet$cov$sig.mat));
  return(as.integer(c(annotated, normalized, compared)));
}


GetMetaTypes <- function(colNm="NA"){
  paramSet <- readSet(paramSet, "paramSet");
  if(colNm=="NA"){
    meta.types <- paramSet$dataSet$meta.types
  }else{
    meta.types <- paramSet$dataSet$meta.types[colNm]
  }
  return(unname(meta.types));
}

GetPrimaryType <- function(analysis.var){
  paramSet <- readSet(paramSet, "paramSet");
  primary.type <- unname(paramSet$dataSet$meta.types[analysis.var]);
  return(primary.type);
}
