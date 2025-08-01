##################################################
## R script for ExpressAnalyst
## Description: Compute ORA heatmap
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## G. Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

# use by GSEA heatmap module as well to generate all the heatmap for the GSEA pathways to select from

# prepare data for heatmap plotting include 
# 1. meta info (top bars)
# 2. expression matrix
# 3. function annotation (left bars)
# all matrix will be combined, 
# 1 and 2 separated by a row of 'null' 
# 3 and 1+2 separated by a column of 'null'
my.prepare.heatmap.json <- function(dataSet, displayOpt="sig"){
  #save.image("prephm.RData")
  data.stat <- qs::qread("data.stat.qs");
  paramSet <- readSet(paramSet, "paramSet");
  res.tbl <- dataSet$comp.res; #dataSet$sig.mat for sig only
  res.tbl <- res.tbl[which(rownames(res.tbl) %in% rownames(data.stat)),]
  sig.ids <- rownames(dataSet$sig.mat);
  if("P.Value" %in% colnames(res.tbl)){
    stat.pvals <- res.tbl$P.Value; 
 }else if("Pval" %in% colnames(res.tbl)){
    stat.pvals <- res.tbl$Pval; 
  }else{
    stat.pvals <- res.tbl$PValue; 
  }

  # Selecting significant gene ids based on data.stat
  all.ids <- rownames(res.tbl)
  all.ids <- all.ids[all.ids %in% rownames(data.stat)]
  # Calculate variance for each gene
  if(displayOpt == "all"){
    #get all
  }else if(displayOpt == "top5000"){
    gene.variances <- apply(data.stat[all.ids, , drop = FALSE], 1, var)
    # Sorting genes by variance and selecting top 5000 if there are more than that
    if (length(gene.variances) > 5000) {
      all.ids <- names(sort(gene.variances, decreasing = TRUE)[1:5000])
    } else {
      all.ids <- all.ids
    }  }else{
    all.ids <- all.ids[all.ids %in% sig.ids];
  }
  
  if("logFC" %in% colnames(res.tbl)){
    stat.fc <- res.tbl$logFC; 
  }else{
    stat.fc <- res.tbl[,paramSet$selectedFactorInx]; 
  }
  
  # scale each gene 
  all.ids <- all.ids[all.ids %in% rownames(data.stat)];
  hit.inz <- sig.ids %in% rownames(data.stat);
  sig.ids <- sig.ids[hit.inz];

  
  #sig only
  #dat <- t(scale(t(data.stat[sig.ids, , drop=F])));
  
  #all genes
  dat <- t(scale(t(data.stat[all.ids, , drop=F])));
  
  # now pearson and euclidean will be the same after scaling
  dat.dist <- dist(dat); 
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  # do clustering and save cluster info
  # convert order to rank (score that can used to sort) 
  if(nrow(dat)> 1){
    dat.dist <- dist(dat);
    gene.ward.ord <- fastcluster::hclust(dat.dist, "ward.D")$order;
    gene.ward.rk <- match(orig.gene.nms, orig.gene.nms[gene.ward.ord]);
    gene.ave.ord <- fastcluster::hclust(dat.dist, "ave")$order;
    gene.ave.rk <- match(orig.gene.nms, orig.gene.nms[gene.ave.ord]);
    gene.single.ord <- fastcluster::hclust(dat.dist, "single")$order;
    gene.single.rk <- match(orig.gene.nms, orig.gene.nms[gene.single.ord]);
    gene.complete.ord <- fastcluster::hclust(dat.dist, "complete")$order;
    gene.complete.rk <- match(orig.gene.nms, orig.gene.nms[gene.complete.ord]);
    
    dat.dist <- dist(t(dat));
    smpl.ward.ord <- fastcluster::hclust(dat.dist, "ward.D")$order;
    smpl.ward.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ward.ord])
    smpl.ave.ord <- fastcluster::hclust(dat.dist, "ave")$order;
    smpl.ave.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ave.ord])
    smpl.single.ord <- fastcluster::hclust(dat.dist, "single")$order;
    smpl.single.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.single.ord])
    smpl.complete.ord <- fastcluster::hclust(dat.dist, "complete")$order;
    smpl.complete.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.complete.ord])
  }else{
    # force not to be single element vector which will be scaler
    stat.pvals <- matrix(stat.pvals);
    stat.fc <- matrix(stat.fc);
    gene.ward.rk <- gene.ave.rk <- gene.single.rk <- gene.complete.rk <- matrix(1);
    smpl.ward.rk <- smpl.ave.rk <- smpl.single.rk <- smpl.complete.rk <- 1:ncol(dat);
}
  
  gene.cluster <- list(
    pval = stat.pvals, 
    fc = stat.fc,
    ward = gene.ward.rk,
    average = gene.ave.rk,
    single = gene.single.rk,
    complete = gene.complete.rk
  );
  
  sample.cluster <- list(
    ward = smpl.ward.rk,
    average = smpl.ave.rk,
    single = smpl.single.rk,
    complete = smpl.complete.rk
  );
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  meta <- data.frame(dataSet$meta.stat);
  grps <- colnames(meta);
  nmeta <- meta.vec <- NULL;
  uniq.num <- 0;
  meta.grps <- vector();
  disc.inx <- rep(F, ncol(meta)*nrow(meta));  
  
  for (i in 1:ncol(meta)){
    cls <- meta[,i];
    grp.nm <- grps[i];
    meta.vec <- c(meta.vec, as.character(cls))
    # make sure each label are unqiue across multiple meta data

    if(dataSet$disc.inx[grp.nm]){
      ncls <- paste(grp.nm, as.numeric(cls)+99); # note, here to retain ordered factor
      disc.inx[c((nrow(meta)*(i-1)+1): (nrow(meta)*i))] <- T;
      sample.cluster[[grps[i]]] <- order(cls);
      
    }else{
      ncls <- rep("NA", length(cls));  
      ncls[!is.na(cls) & cls != "NA"] <- as.numeric(cut(rank(as.numeric(as.character(cls[!is.na(cls) & cls != "NA"]))), breaks=30)); # note, here to retain ordered factor
      ord <- match(orig.smpl.nms, orig.smpl.nms[order(cls)]);
      sample.cluster[[grps[i]]] <- ord;
    }
    meta.grps <- c(meta.grps, paste(grp.nm, rownames(meta))); 
    nmeta <- c(nmeta, ncls);
  }

  # convert back to numeric
  nmeta[disc.inx] <- as.numeric(as.factor(nmeta[disc.inx]))+99;
  unik.inx <- !duplicated(nmeta)     
  
  # get corresponding names
  #meta_anot <- meta.vec[unik.inx]; 
  #names(meta_anot) <- nmeta[unik.inx]; # name annotatation by their numbers
    
  meta_anot <- meta.vec; 
  names(meta_anot) <- meta.grps;

  nmeta <- matrix(nmeta, ncol=ncol(meta), byrow=F);
  colnames(nmeta) <- grps;
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks 
  res <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
  
  # note, use {} will lose order; use [[],[]] to retain the order
  
  sig.inx <- which(rownames(res) %in% sig.ids) -1
  
  if(dataSet$annotated){
    anot.id <- rownames(res);
    gene.map <- readDataQs("symbol.map.qs", paramSet$anal.type, dataSet$name);
    anot.res <- doIdMappingGeneric(anot.id, gene.map, "gene_id", "symbol", "matrix")
    # single element vector will be converted to scalar, not array, need to prevent that
    gene.id <- anot.res$symbol; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
    gene.entrez <- anot.res$gene_id; if(length(gene.entrez) ==1) { gene.entrez <- matrix(gene.entrez) };        
    gene.name <- anot.res$name; if(length(gene.name) ==1) { gene.name <- matrix(gene.name) };
    json.res <- list(
      data.type = dataSet$type, 
      gene.id = anot.res$symbol,
      gene.entrez = gene.entrez,
      gene.name = anot.res$name,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    );
  }else if(file.exists("annotation.qs")){
    # special gene.id and new gene.symbol
    anot.id <- rownames(res);
    anot.res <- doEntrezIDAnot(anot.id, paramSet$data.org, paramSet$data.idType);
    gene.id <- rownames(anot.res); if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
    gene.entrez <- anot.res$gene_id; if(length(gene.entrez) ==1) { gene.entrez <- matrix(gene.entrez) };  
    gene.name <- paste(anot.res$symbol, anot.res$name, sep=" | "); if(length(gene.name) ==1) { gene.name <- matrix(gene.name) };
    
    json.res <- list(
      data.type = dataSet$type, 
      gene.id = gene.id,
      gene.entrez = gene.entrez,
      gene.name = gene.name,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    );
  }else{          
    gene.id <- orig.gene.nms; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
    json.res <- list(
      data.type = dataSet$type,
      gene.id = gene.id,
      gene.entrez = gene.id,
      gene.name = gene.id,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    );
  }
  json.res[["sigInx"]] <- sig.inx;
  return(json.res);
}

# use by GSEA heatmap module as well to generate all the heatmap for the GSEA pathways to select from

# prepare data for heatmap plotting include
# 1. meta info (top bars)
# 2. expression matrix
# 3. function annotation (left bars)
# all matrix will be combined, 
# 1 and 2 separated by a row of 'null' 
# 3 and 1+2 separated by a column of 'null'
my.prepare.metadata.heatmap.json <- function(dataSet,displayOpt="sig"){
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  meta.avgFC <- analSet$meta.avgFC;
  meta.mat <- analSet$meta.mat;
  mdata.all <- paramSet$mdata.all;
  all.meta.mat <- qs::qread("allMeta.mat.qs");
  datanm.vec <- names(mdata.all)[mdata.all==1];
  inmex.meta <- qs::qread("inmex_meta.qs");
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  gene.variances <- apply(inmex.meta$plot.data, 1, var)
  # Sorting genes by variance and selecting top 5000 if there are more than that

  if(displayOpt == "top5000"){
    gene.variances <- apply(inmex.meta$plot.data, 1, var)
    # Sorting genes by variance and selecting top 5000 if there are more than that
    if (length(gene.variances) > 5000) {
      gene.vec <- names(sort(gene.variances, decreasing = TRUE)[1:5000])
    } else{
      gene.vec <- rownames(all.meta.mat);
    }
  } else{
    gene.vec <- rownames(meta.mat);
  }
  dat <- inmex.meta$plot.data[gene.vec, dat.inx, drop=F]; 

  
  # scale each gene for each dataset
  dat <- t(scale(t(dat)));
  
  # now need to remove na or constant rows
  dat <- na.omit(dat);
  # check for columns with all constant (var =0)
  varCol <- apply(dat, 1, var, na.rm=T);
  constCol <- (varCol == 0 | is.na(varCol));
  dat <- dat[!constCol, ];
  
  anot.res <- list();
  ids <- rownames(dat);
  if(inmex.meta$id.type == "entrez"){   
    anot.res <- doEntrezIDAnot(ids, paramSet$data.org, paramSet$data.idType);
  }else{ # no annotation, then use the default feature ID
    anot.res$gene_id <- anot.res$symbol <- anot.res$name <- ids; 
  }
  
  data.nms <- as.factor(inmex.meta$data.lbl[dat.inx]);
  cls.lbls <- as.factor(inmex.meta$cls.lbl[dat.inx]);
  
  # setup annotation info
  data.nms <- as.character(data.nms);
  datasets <- substr(as.character(data.nms), 0, nchar(data.nms)-4);
  
  sel.nms <- names(mdata.all)[mdata.all==1];
  #if(length(sel.nms) > 1){
  #  annotation <- data.frame(class= cls.lbls, dataset = as.factor(datasets));
  #}else{ # single data
  #  annotation <- data.frame(class= cls.lbls);
  #}
  annotation <- paramSet$dataSet$meta.info[dat.inx,];
  rownames(annotation) <- colnames(dat);

  ####
  sig.inx <-which( rownames(dat) %in% rownames(meta.mat)) -1 ;
  if(paramSet$inmex.method != "votecount"){
    stat.pvals <- as.numeric(all.meta.mat[,2]);
    stat.fc <- as.numeric(meta.avgFC[rownames(all.meta.mat)]);
  }else{
    stat.pvals <- as.numeric(all.meta.mat[,1]);
    stat.fc <- as.numeric(meta.avgFC[rownames(all.meta.mat)]);
  }
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  # do clustering and save cluster info
  if(nrow(dat)> 1){
    dat.dist <- dist(dat);
    gene.ward.ord <- hclust(dat.dist, "ward.D")$order;
    gene.ward.rk <- match(orig.gene.nms, orig.gene.nms[gene.ward.ord]);
    gene.ave.ord <- hclust(dat.dist, "ave")$order;
    gene.ave.rk <- match(orig.gene.nms, orig.gene.nms[gene.ave.ord]);
    gene.single.ord <- hclust(dat.dist, "single")$order;
    gene.single.rk <- match(orig.gene.nms, orig.gene.nms[gene.single.ord]);
    gene.complete.ord <- hclust(dat.dist, "complete")$order;
    gene.complete.rk <- match(orig.gene.nms, orig.gene.nms[gene.complete.ord]);
    
    dat.dist <- dist(t(dat));
    smpl.ward.ord <- hclust(dat.dist, "ward.D")$order;
    smpl.ward.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ward.ord])
    smpl.ave.ord <- hclust(dat.dist, "ave")$order;
    smpl.ave.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ave.ord])
    smpl.single.ord <- hclust(dat.dist, "single")$order;
    smpl.single.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.single.ord])
    smpl.complete.ord <- hclust(dat.dist, "complete")$order;
    smpl.complete.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.complete.ord])
  }else{
    # force not to be single element vector which will be scaler
    stat.pvals <- matrix(stat.pvals);
    gene.ward.rk <- gene.ave.rk <- gene.single.rk <- gene.complete.rk <- matrix(1);
    smpl.ward.rk <- smpl.ave.rk <- smpl.single.rk <- smpl.complete.rk <- 1:ncol(dat);
  }
  
  gene.cluster <- list(
    pval = stat.pvals,
    fc = stat.fc, 
    ward = gene.ward.rk,
    average = gene.ave.rk,
    single = gene.single.rk,
    complete = gene.complete.rk
  );
  
  sample.cluster <- list(
    ward = smpl.ward.rk,
    average = smpl.ave.rk,
    single = smpl.single.rk,
    complete = smpl.complete.rk
  );
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  meta <- annotation;
  grps <- colnames(meta);
  nmeta <- meta.vec <- NULL;
  uniq.num <- 0;
  meta.grps <- vector();
  meta.types <- paramSet$dataSet$meta.types;
  disc.inx <- rep(F, ncol(meta)*nrow(meta));  
  
  for (i in 1:ncol(meta)){
    cls <- as.factor(meta[,i]);
    grp.nm <- grps[i];
    meta.vec <- c(meta.vec, as.character(cls))
    # make sure each label are unqiue across multiple meta data

    if( meta.types[grp.nm] == "disc"){
      ncls <- paste(grp.nm, as.numeric(cls)+99); # note, here to retain ordered factor
      disc.inx[c((nrow(meta)*(i-1)+1): (nrow(meta)*i))] <- T;
      sample.cluster[[grps[i]]] <- order(cls);
      
    }else{
      ncls <- as.numeric(cut(rank(as.numeric(as.character(cls))), breaks=30)); # note, here to retain ordered factor
      ord <- match(orig.smpl.nms, orig.smpl.nms[order(cls)]);
      sample.cluster[[grps[i]]] <- ord;
      
    }
    meta.grps <- c(meta.grps, paste(grp.nm, rownames(meta))); 
    nmeta <- c(nmeta, ncls);
  }
  
  # convert back to numeric 
  nmeta <- as.numeric(as.factor(nmeta))+99;
  unik.inx <- !duplicated(nmeta)   
  
  # get corresponding names
  meta_anot <- meta.vec[unik.inx]; 
  names(meta_anot) <- nmeta[unik.inx]; # name annotatation by their numbers
  
  nmeta <- matrix(nmeta, ncol=ncol(meta), byrow=F);
  colnames(nmeta) <- grps;
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks 
  res <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
  
  # note, use {} will lose order; use [[],[]] to retain the order
  # single element vector will be converted to scalar, not array, need to prevent that
  gene.id <- anot.res$symbol; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
  gene.entrez <- anot.res$gene_id; if(length(gene.entrez) ==1) { gene.entrez <- matrix(gene.entrez) };        
  gene.name <- anot.res$name; if(length(gene.name) ==1) { gene.name <- matrix(gene.name) };
  json.res <- list(
    data.type = "array",
    gene.id = as.character(anot.res$symbol),
    gene.entrez = gene.entrez,
    gene.name = anot.res$name,
    gene.cluster = gene.cluster,
    sample.cluster = sample.cluster,
    sample.names = orig.smpl.nms,
    meta = data.frame(nmeta),
    meta.anot = meta_anot,
    data.lbl = datanm.vec,
    data = res,
    sigInx = sig.inx
  );

  return(json.res);
}