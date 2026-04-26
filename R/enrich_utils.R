##################################################
## R script for ExpressAnalyst
## Description: Functions for enrichment analysis (GSEA and ORA)
## Authors:
## G. Zhou, guangyan.zhou@mail.mcgill.ca
## Jeff Xia, jeff.xia@mcgill.ca
###################################################  
.performEnrichAnalysis <- function(dataSet, file.nm, fun.type, ora.vec, vis.type){
  dataSet <<- dataSet;

  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  require(dplyr)
  # prepare lib
  setres <- .loadEnrichLib(fun.type, paramSet)
  current.geneset <- setres$current.geneset;
  current.setids <<- setres$current.setids
  # prepare query
  ora.nms <- names(ora.vec);
  
  if(is.null(ora.nms)){
    ora.nms <- ora.vec;
    names(ora.vec) <- ora.vec;
  }
  
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
  
  # also make sure pathways only contain genes measured in experiment
  #if(!is.null(dataSet$data.anot)){
   if(file.exists("data.anot.qs")){
    current.geneset <- lapply(current.geneset, function(x){x[x %in% current.universe]})
    inds <- lapply(current.geneset, length) > 0
    current.geneset <- current.geneset[inds]
  }

  # Bail out cleanly if no pathways survived filtering — otherwise set.size below
  # becomes NULL and phyper fires "Non-numeric argument to mathematical function".
  # Save an empty enr.mat.qs so downstream steps (ridgeline, enrichnet) that
  # ov_qs_read it see a valid 0-row matrix instead of "file not found".
  if(length(current.geneset) == 0 || length(ora.vec) == 0){
    msgSet$current.msg <- "No matching pathways found between the input gene list and the enrichment library.";
    saveSet(msgSet, "msgSet");
    empty.mat <- matrix(numeric(0), nrow=0, ncol=5);
    colnames(empty.mat) <- c("Total", "Expected", "Hits", "Pval", "FDR");
    ov_qs_save(empty.mat, "enr.mat.qs");
    return(0);
  }

  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval", "FDR");
  
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
  
  gene.vec <- current.universe;
  sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
  gene.nms <- sym.vec;

  current.geneset.symb <- lapply(current.geneset, 
                       function(x) {
                         gene.nms[gene.vec%in%unlist(x)];
  }
  );

  # total unique gene number
  uniq.count <- length(current.universe);
  
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
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    res.mat.all <- as.data.frame(res.mat);
    res.mat.all$Pathway <- rownames(res.mat);
    res.mat.all$Genes <- rep("NA",nrow(res.mat))
    # Iterate through the list and add comma-separated values to the data frame
    for (name in names(hits.query)) {
      if (name %in% res.mat.all$Pathway) {
        res.mat.all[which(res.mat.all$Pathway == name), "Genes"] <- paste(hits.query[[name]], collapse = ",")
      }
    }
    
    res.mat.all <- res.mat.all[which(res.mat.all$Genes != "NA"), ];
    res.mat.all$Pathway <- NULL;
    pws <- rownames(res.mat[which(res.mat.all$Genes != "NA"), ]) 
    fun.ids2 <- as.vector(setres$current.setids[pws]) 
    resTable.all <- data.frame(Pathway = pws, ID = fun.ids2, res.mat.all)

    csv.nm <- paste(file.nm, ".csv", sep="");    
    #print(paste(csv.nm, "=====", "enrichAnalysis"));
    write.csv(resTable.all, file=csv.nm, row.names=F);
    
    imp.inx <- res.mat[,4] <= 0.05;
    imp.inx[is.na(imp.inx)] <- F
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
    msgSet$current.msg <- "No overlap between queried genes and pathway library!"
    # Same as the pre-pipeline guard: persist an empty enr.mat.qs so downstream
    # readers (ridgeline, enrichnet) don't fail with "file not found".
    empty.mat <- matrix(numeric(0), nrow=0, ncol=5);
    colnames(empty.mat) <- c("Total", "Expected", "Hits", "Pval", "FDR");
    ov_qs_save(empty.mat, "enr.mat.qs");
    return(0);
  }
  
  # Check for and handle duplicate row names in enr.mat
  if(any(duplicated(rownames(res.mat)))) {
    res.mat <- res.mat[!duplicated(rownames(res.mat)), ]
    hits.query <- hits.query[match(rownames(res.mat), names(hits.query))]
    print("Duplicates in enr.mat were removed.")
  } else {
    res.mat <- res.mat
  }
  
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  ov_qs_save(res.mat, "enr.mat.qs");
  msgSet$current.msg <- "Functional enrichment analysis was completed";
  
  # write json
  fun.anot <- hits.query;

  # also include original IDs per pathway (for s2f/ko data where symbols differ from fc.log keys)
  hits.ids <- lapply(current.geneset,
                     function(x) {
                       ora.vec[ora.vec %in% unlist(x)];
                     }
  );
  names(hits.ids) <- names(current.geneset);
  hits.ids <- hits.ids[names(fun.anot)];

  total <- resTable$Total; if(length(total) ==1) { total <- matrix(total) };
  fun.pval <- resTable$Pval; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj <- resTable$FDR; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  #print(resTable$Hits);
  hit.num <- paste0(resTable$Hits,"/",resTable$Total); if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(setres$current.setids[resTable$Pathway]);

  resTable$IDs <- fun.ids;
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = setres$current.setlink[1],
    fun.anot = fun.anot,
    fun.anot.ids = hits.ids,
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
  fun.pval <<- fun.pval;
  hit.num <<- resTable$Hits;
  #csv.nm <- paste(file.nm, ".csv", sep="");  
  #print(csv.nm)  
  #fast.write(resTable, file=csv.nm, row.names=F);
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(json.nm))
  
  imgSet <- readSet(imgSet, "imgSet");
  rownames(resTable) <- NULL;
  imgSet$enrTables[[vis.type]] <- list()
  imgSet$enrTables[[vis.type]]$table <- resTable;
  imgSet$enrTables[[vis.type]]$library <- fun.type
  imgSet$enrTables[[vis.type]]$algo <- "Overrepresentation Analysis"

    imgSet$enrTables[[vis.type]]$current.geneset <- current.geneset;
    imgSet$enrTables[[vis.type]]$hits.query <- hits.query;
    imgSet$enrTables[[vis.type]]$current.setids <- current.setids;
    imgSet$enrTables[[vis.type]]$res.mat<- res.mat;
    imgSet$enrTables[[vis.type]]$current.geneset.symb <- current.geneset.symb;
  
    saveSet(imgSet, "imgSet");
  saveSet(paramSet, "paramSet");
  
  saveSet(msgSet, "msgSet");
  return(1);
}

.loadEnrichLib <- function(fun.type, paramSet){
  #if custom return here.
  if(fun.type == "custom"){
    return(.loadCustomEnrichLib(fun.type, paramSet));
  }

  if(paramSet$data.org == "generic"){
    folderNm <- paramSet$data.idType;
  }else{
    folderNm <- paramSet$data.org;
  }

  if(exists("api.lib.path")){
    lib.path <- api.lib.path;
  }else{
    lib.path <- paramSet$lib.path;
  }

  my.path <- paste(lib.path, folderNm, "/", fun.type, ".rds", sep="");
  if(!paramSet$on.public.web && !file.exists(platform.path)){
    nmdb <- basename(my.path);
    download.file(my.path, destfile = nmdb, method="libcurl", mode = "wb");
    my.path <- nmdb;
  }
  
  my.lib <- readRDS(my.path);
  
  if(substr(fun.type, 0, 2)=="go"){  
    if(is.null(names(my.lib))){ # some go lib does not give names
      names(my.lib) <- c("link", "term", "sets");
    }
  }
  
  current.geneset <- my.lib$sets;
  #remove empty pathways
  keep.inx <- lapply(current.geneset,length)>0
  current.geneset <- current.geneset[keep.inx]
  my.lib$term <- my.lib$term[keep.inx]
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- my.lib$term;
  
  if(substr(fun.type, 0, 2)=="go"){
    names(current.geneset) = firstup(names(current.geneset))
    names(current.geneset) = gsub("-", "_", names(current.geneset))
    names(set.ids) = firstup(names(set.ids));
    names(set.ids) = gsub("-", "_", names(set.ids));
  }
  ov_qs_save(current.geneset, "current_geneset.qs");
  res <- list();
  res$current.setlink <- my.lib$link;
  res$current.setids <- set.ids;
  res$current.geneset <- current.geneset;
  return(res);
}

GetRidgePlot <- function(dataName, imgNm = "abc", dpi=default.dpi, format="png", fun.type = "kegg", ridgeType = "ora", ridgeColor = "teal", gseaRankOpt="", sigLevel = 0.05, pwNum=20, inx = 1){
    dataSet <- readDataset(dataName);
    return(compute.ridgeline(dataSet, imgNm, dpi, format, fun.type, ridgeType, ridgeColor,gseaRankOpt, sigLevel, pwNum, inx));
}

PerformUpsetORA <- function(dataName="", file.nm, fun.type, IDs){
  paramSet <- readSet(paramSet, "paramSet");
  dataSet <- readDataset(dataName);
  gene.vec <- unlist(strsplit(IDs, "; "));
  sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
  names(gene.vec) <- sym.vec;
  res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, gene.vec, "upset");
  return(res);
}

PerformGSEA<- function(dataName, file.nm, fun.type, netNm, mType, selectedFactorInx=1, mode = "multi",rankOpt=""){
    return(my.perform.gsea(dataName, file.nm, fun.type, netNm, mType, selectedFactorInx, mode,rankOpt));
}

ComputeRankedVec <- function(data, opt, inx = 1){
   return(my.compute.ranked.vec(data, opt, inx));
}

PlotGSView <-function(cmpdNm, format="png", dpi=default.dpi, width=NA){
   return(plot.gs.view(cmpdNm, format, dpi, width));
}

PlotGSViewNew <-function(cmpdNm, format="png", dpi=default.dpi, imgName){
   return(plot.gs.view(cmpdNm, format, dpi, NA, imgName));
}


.loadCustomEnrichLib <- function(fun.type, paramSet){
  
  # Determine folder name based on paramSet information
  if(paramSet$data.org == "generic"){
    folderNm <- paramSet$data.idType;
  }else{
    folderNm <- paramSet$data.org;
  }

  
  # Load the custom gene set library
  my.lib <- ov_qs_read("custom_lib.qs")

  # Extract the specific gene set based on the function type (e.g., cell line)
  current.geneset <- my.lib
  
  # Remove any empty pathways (or cell lines with no genes)
  keep.inx <- lapply(current.geneset, length) > 0
  current.geneset <- current.geneset[keep.inx]

  # Get the names of the gene sets (e.g., cell line names)
  set.ids <- names(current.geneset)
  names(set.ids) <- names(current.geneset)

  # If the function type pertains to GO terms (or similar), normalize the names
  if(substr(fun.type, 0, 2) == "go") {
    names(current.geneset) <- firstup(names(current.geneset)) # Capitalize first letter
    names(current.geneset) <- gsub("-", "_", names(current.geneset)) # Replace hyphen with underscore
    names(set.ids) <- firstup(names(set.ids))
    names(set.ids) <- gsub("-", "_", names(set.ids))
  }

  # Save the processed gene set to a new file
  ov_qs_save(current.geneset, "current_geneset.qs")
  
  # Create the result object to return
  res <- list()
  res$current.setlink <- "" # Empty placeholder for set link
  res$current.setids <- set.ids # Names (IDs) of the gene sets (cell lines)
  res$current.geneset <- current.geneset # The actual gene set data

  return(res)
}

PerformDefaultEnrichment <- function(dataName="", file.nm, fun.type){
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");
  #save.image("defaultenr.RData");
  anal.type <- paramSet$anal.type;
    if(anal.type=="onedata"){
      dataSet <- readDataset(dataName); #instead of function parameter
      if(nrow(dataSet$sig.mat) == 0){
        print("No DE genes were identified, can not perform ORA analysis!");
        return(0);
      }

      gene.vec <- rownames(dataSet$sig.mat);
    }else if(anal.type=="metadata"){
      gene.vec <- rownames(analSet$meta.mat);
    }else{
      gene.vec <- rownames(paramSet$all.ent.mat);
    }
  sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
  names(gene.vec) <- sym.vec;
  res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, gene.vec, "default");

  
  return(1);
}

GetSigSetCount <- function(enrType, type, pval=0.05){
  pval <- as.numeric(pval);
  imgSet <- readSet(imgSet, "imgSet");
    
  tbl <- imgSet$enrTables[[enrType]]$table
  count <- 0;
  if(type == "raw"){
   count<-sum(tbl$Pval<pval);
  }else{
    count<-sum(tbl$FDR<pval);
  }
  return(count);
}


GetSetIDLinks <- function(type=""){
  imgSet <- readSet(imgSet, "imgSet");
  fun.type <- imgSet$enrTables[[type]]$library;

  paramSet <- readSet(paramSet, "paramSet");
  ids <- imgSet$enrTables[[type]]$table$IDs
  pathways <- imgSet$enrTables[[type]]$table$Pathway
    if(fun.type %in% c("go_bp", "go_mf", "go_cc")){
        annots <- paste("<a href='https://www.ebi.ac.uk/QuickGO/term/", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type %in% c("go_panthbp", "go_panthmf", "go_panthcc")){
        annots <- paste("<a href='https://www.pantherdb.org/panther/categoryList.do?searchType=basic&fieldName=all&organism=all&fieldValue=", ids, "&listType=5' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "kegg"){
        annots <- paste("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway+", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "reactome"){
        annots <- paste("<a href='https://reactome.org/content/query?q=", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else{
        annots <- ids;
    }
  
  return(annots);
}

GetHTMLPathSet <- function(type, setNm){
  imgSet <- readSet(imgSet, "imgSet");
  current.geneset <- imgSet$enrTables[[type]]$current.geneset.symb;
  hits.query <- imgSet$enrTables[[type]]$hits.query;
  set <- current.geneset[[setNm]]; 
  
  #set <- cur.setids[[setNm]];
  
  hits <- hits.query
  
  # highlighting with different colors
  red.inx <- which(set %in% hits[[setNm]]);
  
  # use actual cmpd names
  #nms <- names(set);
  nms <- doEntrez2SymbolMapping(set);
  nms[red.inx] <- paste("<font color=\"red\">", "<b>", nms[red.inx], "</b>", "</font>",sep="");

  return(cbind(setNm, paste(unique(nms), collapse="; ")));
}


GetEnrResultMatrix <-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  #print(names(imgSet$enrTables[[type]]));
  res <- imgSet$enrTables[[type]]$res.mat
  return(signif(as.matrix(res), 5));
}

GetEnrResultColNames<-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[[type]]$res.mat
  colnames(res);
}

GetEnrResSetIDs<-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[[type]]$table;
  return(res$IDs);
}

GetEnrResSetNames<-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[[type]]$table;
  return(res$Pathway);
}


GetGseaSetCount <- function(type, pval=0.05){
  pval <- as.numeric(pval);
  imgSet <- readSet(imgSet, "imgSet");
    
  tbl <- imgSet$enrTables[["gsea"]]$table
  count <- 0;
  if(type == "raw"){
   count<-sum(tbl$Pval<pval);
  }else{
    count<-sum(tbl$Padj<pval);
  }
  return(count);
}


GetGseaIDLinks <- function(dataName=""){
  dataSet <- readDataset(dataName);
  imgSet <- readSet(imgSet, "imgSet");
  fun.type <- imgSet$enrTables[["gsea"]]$library;

  paramSet <- readSet(paramSet, "paramSet");
  ids <- imgSet$enrTables[["gsea"]]$table$IDs
  pathways <- imgSet$enrTables[["gsea"]]$table$Name

    if(fun.type %in% c("go_bp", "go_mf", "go_cc")){
        annots <- paste("<a href='https://www.ebi.ac.uk/QuickGO/term/", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type %in% c("go_panthbp", "go_panthmf", "go_panthcc")){
        annots <- paste("<a href='https://www.pantherdb.org/panther/categoryList.do?searchType=basic&fieldName=all&organism=all&fieldValue=", ids, "&listType=5' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "kegg"){
        annots <- paste("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway+", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "reactome"){
        annots <- paste("<a href='https://reactome.org/content/query?q=", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else{
        annots <- ids;
    }
  
  return(annots);
}

GetGseaHTMLPathSet <- function(setNm){
  imgSet <- readSet(imgSet, "imgSet");
  current.geneset <- imgSet$enrTables[["gsea"]]$current.geneset.symb;
  hits.query <- imgSet$enrTables[["gsea"]]$hits.query;
  set <- current.geneset[[setNm]]; 
    
  hits <- hits.query
  
  # highlighting with different colors
  red.inx <- which(set %in% hits[[setNm]]);
  
  # use actual cmpd names
  #nms <- names(set);
  nms <- set;
  nms[red.inx] <- paste("<font color=\"red\">", "<b>", nms[red.inx], "</b>", "</font>",sep="");

  return(cbind(setNm, paste(unique(nms), collapse="; ")));
}


GetGseaResultMatrix <-function(){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[["gsea"]]$res.mat
  # Return an empty numeric matrix when GSEA has not yet been run (or when the
  # enrichment table is empty). Avoids as.matrix(NULL) throwing, which Rserve
  # surfaces as an opaque error code 127 with no usable message on the caller side.
  if (is.null(res) || (is.null(dim(res)) && length(res) == 0)) {
    return(matrix(numeric(0), 0, 0));
  }
  return(signif(as.matrix(res), 5));
}

GetGseaResultColNames <-function(){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[["gsea"]]$res.mat
  colnames(res);
}

GetGseaResSetIDs <-function(){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[["gsea"]]$table;

  return(res$IDs);
}

GetGseaResSetNames<-function(){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[["gsea"]]$table;
  return(res$Name);
}

GetEnrResTypes<-function(){
  imgSet <- readSet(imgSet, "imgSet");
  nms <- names(imgSet$enrTables);
  nms <- setdiff(nms, "default");
  nms <- setdiff(nms, "gsea")
  return(nms);
}

GetEnrResLibrary<-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  summary <- imgSet$enrTables[[type]];
  res <- summary$library
  return(res);
}

# ── Server-side PNG: Enrichment Network (igraph) ──
PlotEnrichNetworkPNG <- function(dataName, imgName, format="png", dpi=150, width=NA) {
  tryCatch({
    require(igraph); require(reshape)
    cat("[PlotEnrichNetworkPNG] wd=", getwd(), "\n")
    cat("[PlotEnrichNetworkPNG] enr.mat.qs exists:", file.exists("enr.mat.qs"), "\n")
    # List qs files in current directory
    qs_files <- list.files(pattern = "\\.qs$")
    cat("[PlotEnrichNetworkPNG] qs files in wd:", paste(qs_files, collapse=", "), "\n")
    enr.mat <- qs::qread("enr.mat.qs")
    hits.query <- qs::qread("hits_query.qs")
    if (is.null(enr.mat) || nrow(enr.mat) == 0) return(0)
    if ("FDR" %in% colnames(enr.mat)) {
      ord.inx <- order(enr.mat[, "FDR"])
      enr.mat <- enr.mat[ord.inx[1:min(20, nrow(enr.mat))], , drop = FALSE]
    }
    hits <- as.numeric(enr.mat[, "Hits"]); id <- rownames(enr.mat)
    names(hits) <- id; hits.query <- hits.query[id]
    n <- nrow(enr.mat); w <- matrix(0, n, n); colnames(w) <- rownames(w) <- id
    for (i in 1:n) for (j in i:n) w[i,j] <- overlap_ratio(hits.query[id[i]], hits.query[id[j]], "mixed")
    w[lower.tri(w)] <- t(w)[lower.tri(w)]
    wd <- reshape::melt(w); wd <- wd[wd[,1]!=wd[,2] & !is.na(wd[,3]),]
    g <- graph_from_data_frame(wd[,-3], directed=FALSE)
    g <- delete_edges(g, E(g)[wd[,3] < 0.3])
    if (vcount(g) == 0) return(0)
    cnt2 <- hits[V(g)$name]
    V(g)$size <- if (all(cnt2==cnt2[1])) rep(12, length(cnt2)) else scales::rescale(log(cnt2+1,10), to=c(8,24))
    V(g)$color <- "orange"; V(g)$frame.color <- "white"
    V(g)$label <- V(g)$name; V(g)$label.cex <- 0.6; V(g)$label.color <- "black"
    E(g)$arrow.mode <- 0; E(g)$color <- "#cccccc"
    l <- layout_with_graphopt(g)
    imgPath <- paste0(imgName, ".", format)
    w.val <- if (is.na(width)) 8 else width/dpi
    png(imgPath, width=w.val, height=w.val*0.75, units="in", res=dpi)
    par(mar=c(1,1,2,1)); plot(g, layout=l, main="Enrichment Network (KEGG)"); dev.off()
    return(1)
  }, error = function(e) { message("PlotEnrichNetworkPNG error: ", e$message); return(0) })
}

# ── Server-side PNG: Gene-Pathway Enrichment Heatmap ──
PlotEnrichHeatmapPNG <- function(dataName, imgName, format="png", dpi=150, width=NA) {
  tryCatch({
    enr.mat <- qs::qread("enr.mat.qs")
    current.geneset <- if (file.exists("current_geneset.qs")) qs::qread("current_geneset.qs") else NULL
    if (is.null(enr.mat) || nrow(enr.mat) < 2 || is.null(current.geneset)) return(0)
    dataSet <- readDataset(dataName)
    if (is.null(dataSet) || is.null(dataSet$prot.mat)) return(0)
    user.entrez <- rownames(dataSet$prot.mat)  # EA gene lists: prot.mat has Entrez IDs
    if ("FDR" %in% colnames(enr.mat)) {
      enr.mat <- enr.mat[order(as.numeric(enr.mat[,"FDR"]))[1:min(15,nrow(enr.mat))], , drop=FALSE]
    }
    matched.pw <- intersect(rownames(enr.mat), names(current.geneset))
    if (length(matched.pw) < 2) return(0)
    hit.genes.per.pw <- lapply(matched.pw, function(pw) intersect(user.entrez, current.geneset[[pw]]))
    names(hit.genes.per.pw) <- matched.pw
    all.hit.genes <- unique(unlist(hit.genes.per.pw))
    if (length(all.hit.genes) < 2) return(0)
    gp.mat <- matrix(0, length(all.hit.genes), length(matched.pw))
    rownames(gp.mat) <- all.hit.genes; colnames(gp.mat) <- matched.pw
    for (pw in matched.pw) { g <- hit.genes.per.pw[[pw]]; if (length(g)>0) gp.mat[g,pw] <- 1 }
    tryCatch({
      ps <- readSet(paramSet, "paramSet")
      syms <- doEntrez2SymbolMapping(rownames(gp.mat), ps$data.org, ps$data.idType)
      if (!is.null(syms) && length(syms)==nrow(gp.mat)) rownames(gp.mat) <- make.unique(syms)
    }, error=function(e) {})
    if (nrow(gp.mat) > 30) gp.mat <- gp.mat[order(rowSums(gp.mat), decreasing=TRUE)[1:30], , drop=FALSE]
    gp.mat <- gp.mat[rowSums(gp.mat)>0, colSums(gp.mat)>0, drop=FALSE]
    if (nrow(gp.mat)<2 || ncol(gp.mat)<2) return(0)
    colnames(gp.mat) <- make.unique(substr(colnames(gp.mat),1,40))
    imgPath <- paste0(imgName, ".", format)
    w.val <- if (is.na(width)) max(7, ncol(gp.mat)*0.6+3) else width/dpi
    h.val <- max(5, nrow(gp.mat)*0.25+2)
    png(imgPath, width=w.val, height=h.val, units="in", res=dpi, bg="white")
    par(mar=c(1, 8, max(4, max(nchar(colnames(gp.mat)))*0.3), 1))
    nr <- nrow(gp.mat); nc <- ncol(gp.mat)
    plot(NA, xlim=c(0,nc), ylim=c(0,nr), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", asp=NA)
    for (i in 1:nr) for (j in 1:nc) if (gp.mat[i,j]==1) rect(j-1,nr-i,j-0.05,nr-i+0.95, col="#F4837D", border="white", lwd=0.5)
    for (i in 0:nr) abline(h=i, col="#e0e0e0", lwd=0.3); for (j in 0:nc) abline(v=j, col="#e0e0e0", lwd=0.3)
    mtext(rownames(gp.mat), side=2, at=nr:1-0.5, las=1, line=0.5, cex=0.7, adj=1)
    text(x=1:nc-0.5, y=nr+0.3, labels=colnames(gp.mat), srt=45, adj=0, xpd=TRUE, cex=0.65)
    mtext("Input Genes", side=2, line=6.5, cex=0.9, font=2)
    mtext("Enriched Terms", side=3, line=max(2, max(nchar(colnames(gp.mat)))*0.15), cex=0.9, font=2)
    dev.off(); return(1)
  }, error = function(e) { message("PlotEnrichHeatmapPNG error: ", e$message); tryCatch(dev.off(), error=function(x){}); return(0) })
}