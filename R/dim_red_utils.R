##################################################
## R script for ExpressAnalyst
## Description: Computing PCA coordinates
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

SaveClusterJSONLoading <- function(dataName="", fileNm, clustOpt, nb){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  if(anal.type == "onedata"){
    .saveExpressClusterLoadingJSON(dataSet, fileNm, clustOpt,paramSet, nb);
  }else{
    .saveMetaClusterLoadingJSON(dataSet, fileNm, clustOpt,paramSet, nb);
  }
}

SaveClusterJSON <- function(dataName="", fileNm, clustOpt, opt){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  if(anal.type == "onedata"){
    .saveExpressClusterJSON(dataSet, fileNm, clustOpt, paramSet ,opt);
  }else{
    .saveMetaClusterJSON(dataSet, fileNm, clustOpt,paramSet , opt);
  }
}

.saveMetaClusterJSON <- function(dataSet, fileName, clustOpt,paramSet, opt){
    
    msgSet <- readSet(msgSet, "msgSet");
    paramSet <- readSet(paramSet, "paramSet");
    analSet <- readSet(analSet, "analSet");

    mdata.all <- paramSet$mdata.all;

    inmex.meta <- ov_qs_read("inmex_meta.qs");
    datanm.vec <- names(mdata.all)[mdata.all==1];

    dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
    dat <- inmex.meta$data[, dat.inx, drop=F]; 

    # need to deal with missing values 
    dat <- na.omit(dat);

    pca3d <- list();
    if(clustOpt == "pca"){
        if(opt == "all"){
            pca <- prcomp(t(dat), center=T, scale=T);
            }else{
            dat <- dat[which(rownames(dat) %in% analSet$loadEntrez),]
            pca <- prcomp(t(dat), center=T, scale=T);
            }

        imp.pca<-summary(pca)$importance;
        pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
        coords <- data.frame(t(signif(pca$x[,1:3], 5)));
    }else if(clustOpt == "umap"){
        require('uwot');
        # OPTIMIZED: Calculate ncol once
        n_cols <- ncol(dat)
        if(n_cols < 100){
            neighbor_num <- n_cols
        }else{
            neighbor_num <- 100;
        }

        ndat <- as.matrix(t(dat));
        res <- umap(ndat, n_components=3, n_neighbors=neighbor_num);
        pca3d$score$axis <- paste("UMAP dim ", 1:3, sep="");
        coords <- data.frame(t(signif(res, 5)));
    }else{
        require('Rtsne');
        ndat <- as.matrix(t(dat));
        max.perx <- floor((nrow(ndat)-1)/3);
        if(max.perx > 30){
            max.perx <- 30;
        }
        res <- Rtsne(ndat, dims = 3, perplexity=max.perx);
        pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
        coords <- data.frame(t(signif(res$Y, 5)));
    }

    colnames(coords) <- NULL; 
    pca3d$score$xyz <- coords;
    pca3d$score$name <- colnames(dat);

    facA <- as.character(inmex.meta$cls.lbl[dat.inx]);
    if(all.numeric(facA)){
        facA <- paste("Group", facA);
    }
    pca3d$score$facA <- facA;

    facB <-  as.character(inmex.meta$data.lbl[dat.inx]);
    if(all.numeric(facB)){
        facB <- paste("Group", facB);
    }
    pca3d$score$facB <- facB;

    # now set color for each group
    cols <- unique(GetColorSchema(facB));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
    pca3d$score$colors <- cols;

    # add shape sphere, triangles, square, pentagon (first two)
    pca3d$score$shapes <- c("sphere", "triangle");

    mypos <- t(coords);
    colnames(mypos) <- paste("Dim", 1:3, sep="");
    coords <- data.frame(Class=facA, Data=facB, mypos);

    pos.xyz <- mypos;
    pos.xyz <- unitAutoScale(pos.xyz);
    rownames(pos.xyz) = pca3d$score$name;
    ov_qs_save(pos.xyz, "score_pos_xyz.qs");

    fast.write(coords, file="expressanalyst_3d_pos.csv");

    pca3d$org <- paramSet$data.org
    pca3d$analType <- paramSet$anal.type
    pca3d$naviString <- "Scatter 3D"
    ov_qs_save(pca3d, "pca3d.qs");

    paramSet$jsonNms$pcascore <- fileName
    paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(fileName))
    # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
    jsonlite::write_json(pca3d, fileName, auto_unbox = TRUE, pretty = FALSE);
    msgSet$current.msg <- "Annotated data is now ready for 3D visualization!";
    saveSet(msgSet, "msgSet");
    saveSet(paramSet, "paramSet");

    return(1);
}


.saveMetaClusterLoadingJSON <- function(dataSet, fileName, clustOpt, paramSet, nb){
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  mdata.all <- paramSet$mdata.all;

  inmex.meta <- ov_qs_read("inmex_meta.qs");
  datanm.vec <- names(mdata.all)[mdata.all==1];
  nb <- as.numeric(5000) # set to max 5000 datapoints
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  dat <- inmex.meta$data[, dat.inx, drop=F]; 
  
  # need to deal with missing values 
  dat <- na.omit(dat);
  variances <- apply(dat,1, function(x){var(x)})
  df <- data.frame(var = variances, inx = seq.int(1,length(variances)))
  df <- df[order(-df$var),];

  #do not take subset of loading data points now
  if(nb < length(df$inx)){
    inx <- df$inx[c(1:nb)];
  }else{
    inx <- df$inx;
  }
  dat <- dat[inx,];
  
  pca3d <- list();
  
  pca <- prcomp(t(dat), center=T, scale=T);    
  imp.pca<-summary(pca)$importance;
  pca3d$score$axis <- paste("PC", 1:3, sep="");
  coords <- data.frame(t(signif(pca$rotation[,1:3], 5)));
  
  colnames(coords) <- NULL; 
  pca3d$score$xyz <- coords;
  pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation), paramSet$data.org, paramSet$data.idType);
  pca3d$score$entrez <- rownames(pca$rotation);
  
  analSet$loadEntrez <- pca3d$score$entrez
  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  rownames(mypos) <- analSet$loadEntrez;
  mypos <- unitAutoScale(mypos);
  ov_qs_save(mypos, "loading_pos_xyz.qs");
  
  coords <- data.frame(mypos);
  fast.write(coords, file="expressanalyst_loadings_3d_pos.csv");
  
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(fileName))
  paramSet$jsonNms$pcaload <- fileName;
  ov_qs_save(pca3d, "pca3d.qs");
  # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
  jsonlite::write_json(pca3d, fileName, auto_unbox = TRUE, pretty = FALSE);
  msgSet$current.msg <- "Annotated data is now ready for 3D visualization!";
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  saveSet(analSet, "analSet");

  return(1);
}


.saveExpressClusterLoadingJSON <- function(dataSet, fileName, clustOpt, paramSet, nb){  
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  dat <- dataSet$data.norm;
  pca3d <- list();
  dat <- na.omit(dat);
  nb <- as.numeric(5000) # set to max 5000 datapoints
  if(clustOpt == "pca"){
    pca <- prcomp(t(dat), center=T, scale=T);
    imp.pca<-summary(pca)$importance;
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(pca$rotation[,1:3], 5)));
    
    colnames(coords) <- NULL; 
    pca3d$score$xyz <- coords;
    pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation), paramSet$data.org, paramSet$data.idType);
    pca3d$score$entrez <-rownames(pca$rotation);
    weights <- imp.pca[2,][1:3]
    mypos <- t(coords);
    meanpos <- apply(abs(mypos),1, function(x){weighted.mean(x, weights)})
    df <- data.frame(pos = meanpos, inx = seq.int(1,length(meanpos)))
    df <- df[order(-df$pos),]
    
    if(nrow(df) > nb){
      inx <- df$inx[c(1:nb)]
      mypos <- mypos[inx,];
      pca3d$score$xyz <- coords[inx]
      pca3d$score$name <- pca3d$score$name[inx]
      pca3d$score$entrez <- pca3d$score$entrez[inx]
    }
  }

  pca3d$cls <- dataSet$meta.info;
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  # see if there is secondary
  analSet$loadEntrez <- pca3d$score$entrez
  rownames(mypos) <- pca3d$score$name;
  rownames(mypos) <- analSet$loadEntrez;
  mypos <- unitAutoScale(mypos);
  ov_qs_save(mypos, "loading_pos_xyz.qs");
  
  fast.write(mypos, file="expressanalyst_3d_load_pos.csv");
  ov_qs_save(pca3d, "pca3d.qs");
  paramSet$jsonNms$pcaload <- fileName
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(fileName))
  # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
  jsonlite::write_json(pca3d, fileName, auto_unbox = TRUE, pretty = FALSE);
  msgSet$current.msg <- "Annotated data is now ready for PCA 3D visualization!";
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  saveSet(analSet, "analSet");

  return(1);
}

# single expression data
.saveExpressClusterJSON <- function(dataSet, fileName, clustOpt,paramSet, opt){
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  dat <- dataSet$data.norm;
  pca3d <- list();
  dat <- na.omit(dat);
  
  if(clustOpt == "pca"){
    #if(opt == "all"){
      pca <- prcomp(t(dat));
   #}else{
    #  dat <- dat[which(rownames(dat) %in% analSet$loadEntrez),]
    #  pca <- prcomp(t(dat), center=T, scale=T);
    #}
    imp.pca<-summary(pca)$importance;
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(pca$x[,1:3], 5)));
  }else if(clustOpt == "umap"){
    require('uwot');
    # OPTIMIZED: Calculate ncol once
    n_cols <- ncol(dat)
    if(n_cols < 100){
      neighbor_num <- n_cols
    }else{
      neighbor_num <- 100;
    }
    dat <- as.matrix(t(dat));
    res <- umap(dat, n_components=3, n_neighbors=neighbor_num);
    pca3d$score$axis <- paste("UMAP dim ", 1:3, sep="");
    coords <- data.frame(t(signif(res, 5)));
    
  }else{ # tsne
    require('Rtsne');
    dat <- as.matrix(t(dat));
    max.perx <- floor((nrow(dat)-1)/3);
    if(max.perx > 30){
      max.perx <- 30;
    }
    res <- Rtsne(dat, dims = 3, perplexity=max.perx);
    pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
    coords <- data.frame(t(signif(res$Y, 5)));
  }
  
  colnames(coords) <- NULL; 
  pca3d$score$xyz <- coords;
  pca3d$score$name <- colnames(dataSet$data.norm);
  
  pos.xyz <- data.frame(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3]);
  pos.xyz <- as.data.frame(pos.xyz);
  pos.xyz <- unitAutoScale(pos.xyz);
  rownames(pos.xyz) = colnames(dataSet$data.norm);
  ov_qs_save(pos.xyz, "score_pos_xyz.qs");
  
  facA <- as.character(dataSet$fst.cls);
  if(all.numeric(facA)){
    facA <- paste("Group", facA);
  }
  pca3d$score$facA <- facA;
  
  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  # see if there is secondary
  if(length(dataSet$sec.cls) > 1){
    facB <- as.character(dataSet$sec.cls);
    if(all.numeric(facB)){
      facB <- paste("Group", facB);
    }
    pca3d$score$facB <- facB;
    
    # set shape based on the first group
    pca3d$score$shapes <- c("sphere", "triangle");
    
    # now set color based on 2nd group
    cols <- unique(GetColorSchema(facB));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
    pca3d$score$colors <- cols;
    
    mypos <- data.frame(factorA=facA, factorB=facB, mypos);
  }else{
    # now set color based on first group
    cols <- unique(GetColorSchema(facA));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
    pca3d$score$colors <- cols;
    mypos <- data.frame(factorA=facA, mypos);
  }
  
  pca3d$cls <- dataSet$meta.info;
  pca3d$org <- paramSet$data.org
  pca3d$analType <- paramSet$anal.type
  pca3d$naviString <- "Scatter 3D"
  ov_qs_save(pca3d, "pca3d.qs");
  paramSet$jsonNms$pcascore <- fileName
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(fileName))
  rownames(mypos) <- colnames(dataSet$data.norm);
  
  fast.write(mypos, file="expressanalyst_3d_pos.csv");
  # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
  jsonlite::write_json(pca3d, fileName, auto_unbox = TRUE, pretty = FALSE);
  msgSet$current.msg <- "Annotated data is now ready for PCA 3D visualization!";
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  return(1);
}


#' Compute ellipsoid encasings for multiple groups in one subprocess call.
#'
#' @param filenm       Output JSON filename.
#' @param type         Encasing type (only "ellipse" supported).
#' @param groups.json  JSON array of groups; each element has fields
#'                     {grpName: "...", names: "id1; id2; id3"}.
#' @param level        Confidence level (default 0.95).
#' @param omics        Reserved for multi-omics callers (ignored here).
#' @return filenm
#' @export
ComputeEncasingBatch <- function(filenm, type, groups.json, level = 0.95, omics = "NA") {
  tryCatch({
    level <- as.numeric(level)

    if (!file.exists("score_pos_xyz.qs") && !file.exists("score_pos_xyz.qs2")) {
      sink(filenm); cat("{}"); sink()
      return(filenm)
    }
    pos.xyz <- ov_qs_read("score_pos_xyz.qs")

    groups_list <- RJSONIO::fromJSON(groups.json)
    if (is.data.frame(groups_list)) {
      groups_list <- split(groups_list, seq_len(nrow(groups_list)))
    }

    # Parse groups (expects fields grpName + names, matching Express's JS callers)
    # and collect per-group xyz coord matrices in master.
    coords_per_group <- vector("list", length(groups_list))
    group_names      <- character(length(groups_list))

    for (i in seq_along(groups_list)) {
      g <- groups_list[[i]]
      if (is.character(g)) {
        group_names[i] <- unname(g["grpName"])
        names_vec      <- unname(g["names"])
      } else if (is.data.frame(g)) {
        group_names[i] <- g$grpName[1]
        names_vec      <- g$names[1]
      } else {
        group_names[i] <- g$grpName
        names_vec      <- g$names
      }
      nms <- strsplit(names_vec, "; ")[[1]]
      inx <- rownames(pos.xyz) %in% nms
      coords_per_group[[i]] <- as.matrix(pos.xyz[inx, c(1:3)])
    }

    # Offload all ellipsoid computations to a single rsclient subprocess round-trip.
    bridge_in  <- ov_bridge_file("in")
    bridge_out <- sub("_in.qs2", "_out.qs2", bridge_in)
    ov_qs_save(list(coords_per_group = coords_per_group, level = level,
                    group_names = group_names), bridge_in)
    on.exit(unlink(c(bridge_in, bridge_out)), add = TRUE)

    run_func_via_rsclient(
      func = function(wd, bridge_in, bridge_out) {
        setwd(wd)
        Sys.setenv(RGL_USE_NULL = TRUE)
        require(rgl)
        input <- ov_qs_read(bridge_in)
        coords_list <- input$coords_per_group
        level <- input$level
        group_names <- input$group_names

        result_list <- vector("list", length(coords_list))
        for (i in seq_along(coords_list)) {
          coords <- coords_list[[i]]
          if (nrow(coords) < 4) {
            result_list[[i]] <- list(grpName = group_names[i], mesh = list(),
                                     error = "Insufficient points")
            next
          }
          tryCatch({
            pos    <- cov(coords, y = NULL, use = "everything")
            center <- colMeans(coords)
            t_val  <- sqrt(qchisq(level, 3))
            mesh   <- list()
            mesh[[1]] <- rgl::ellipse3d(x = as.matrix(pos), centre = center, t = t_val)
            result_list[[i]] <- list(grpName = group_names[i], mesh = mesh, error = NULL)
          }, error = function(e) {
            result_list[[i]] <<- list(grpName = group_names[i], mesh = list(),
                                      error = e$message)
          })
        }
        ov_qs_save(result_list, bridge_out)
      },
      args = list(wd = getwd(), bridge_in = bridge_in, bridge_out = bridge_out),
      timeout_sec = 120
    )

    mesh_results <- if (file.exists(bridge_out)) ov_qs_read(bridge_out) else NULL
    if (!is.null(mesh_results)) {
      sink(filenm); cat(RJSONIO::toJSON(mesh_results)); sink()
    } else {
      sink(filenm); cat("{}"); sink()
    }
  }, error = function(e) {
    message("[ComputeEncasingBatch] ", e$message)
    sink(filenm); cat("{}"); sink()
  })
  return(filenm)
}


ComputeEncasing <- function(filenm, type, names.vec, level=0.95, omics="NA"){
  tryCatch({
    level <- as.numeric(level)
    names <- strsplit(names.vec, "; ")[[1]]

    if (!file.exists("score_pos_xyz.qs")) {
      sink(filenm); cat("{}"); sink()
      return(filenm)
    }
    pos.xyz <- ov_qs_read("score_pos_xyz.qs")
    inx <- rownames(pos.xyz) %in% names
    coords <- as.matrix(pos.xyz[inx, c(1:3)])

    if (nrow(coords) < 4) {
      sink(filenm); cat(RJSONIO::toJSON(list())); sink()
      return(filenm)
    }

    bridge_in <- ov_bridge_file("in")
    bridge_out <- sub("_in.qs2", "_out.qs2", bridge_in)
    ov_qs_save(list(coords = coords, level = level), bridge_in, preset = "fast")
    on.exit(unlink(c(bridge_in, bridge_out)), add = TRUE)

    run_func_via_rsclient(
      func = function(wd, bridge_in, bridge_out) {
        setwd(wd)
        Sys.setenv(RGL_USE_NULL = TRUE)
        require(rgl)
        input <- ov_qs_read(bridge_in)
        pos <- cov(input$coords, y = NULL, use = "everything")
        center <- colMeans(input$coords)
        t_val <- sqrt(qchisq(input$level, 3))
        mesh <- list()
        mesh[[1]] <- rgl::ellipse3d(x = as.matrix(pos), centre = center, t = t_val)
        ov_qs_save(mesh, bridge_out, preset = "fast")
      },
      args = list(wd = getwd(), bridge_in = bridge_in, bridge_out = bridge_out),
      timeout_sec = 120
    )

    mesh <- if (file.exists(bridge_out)) ov_qs_read(bridge_out) else NULL
    if (!is.null(mesh)) {
      sink(filenm); cat(RJSONIO::toJSON(mesh)); sink()
    }
  }, error = function(e) {
    message("[ComputeEncasing] ", e$message)
    sink(filenm); cat("{}"); sink()
  })
  return(filenm)
}

unitAutoScale <- function(df){
    df <- as.data.frame(df)
    row.nms <- rownames(df);
    col.nms <- colnames(df);
    df<-apply(df, 2, AutoNorm);
    rownames(df) <- row.nms;
    colnames(df) <- col.nms;
    maxVal <- max(abs(df))
    df<- df/maxVal
    return(df)
}
