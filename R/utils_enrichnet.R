

my.enrich.net<-function(dataSet, netNm="abc", type="list", overlapType="mixed", analSet){
  enr.mat <- ov_qs_read("enr.mat.qs");

  # Filter by adjusted p-value (FDR) < 0.05 and limit to max 50 pathways
  # If fewer than 20 are significant (or none), show top 20 by FDR
  if("FDR" %in% colnames(enr.mat)){
    sig.inx <- enr.mat[,"FDR"] < 0.05;
    n.sig <- sum(sig.inx);

    if(n.sig >= 20){
      # At least 20 significant pathways exist
      enr.mat <- enr.mat[sig.inx, , drop=FALSE];

      # Limit to top 50 pathways (sorted by FDR)
      if(nrow(enr.mat) > 50){
        ord.inx <- order(enr.mat[,"FDR"]);
        enr.mat <- enr.mat[ord.inx[1:50], , drop=FALSE];
      }
    } else {
      # Fewer than 20 significant pathways (or none), show top 20 by FDR
      ord.inx <- order(enr.mat[,"FDR"]);
      max.show <- min(20, nrow(enr.mat));
      enr.mat <- enr.mat[ord.inx[1:max.show], , drop=FALSE];
    }
  }

  hits <-  enr.mat[,"Hits"];
  pvals <- enr.mat[,"Pval"];

  pvalue <- pvals;
  id <- names(pvalue);

  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;

  if(is.null(enr.mat) || nrow(enr.mat) == 0){
    return(0);
  }
  
  require(igraph);
  require(reshape);

  current.geneset <- ov_qs_read("current_geneset.qs");
  hits.query <- ov_qs_read("hits_query.qs")
  # Filter hits.query to match the filtered pathways (FDR < 0.05, max 50)
  hits.query <- hits.query[rownames(enr.mat)];
  geneSets <- hits.query;
  n <- nrow(enr.mat);
  w  <- matrix(NA_real_, n, n);   # overlap ratio (edge filter)
  cw <- matrix(0L, n, n);         # shared-gene COUNT (edge weight -> thickness)
  colnames(w) <- rownames(w) <- id;
  # Upper triangle only — compute each pathway PAIR once.
  for (i in 1:n) for (j in i:n) {
    if (i == j) next
    gi <- unlist(geneSets[id[i]]); gj <- unlist(geneSets[id[j]]);
    w[i, j]  <- overlap_ratio(gi, gj, overlapType);
    cw[i, j] <- length(intersect(gi, gj));
  }
  # ONE edge per pair (the symmetric melt used to create duplicate i->j AND j->i
  # edges, crowding the map with parallel edges). Keep pairs with overlap >= 0.3;
  # edge weight = number of shared genes so the static map can scale thickness.
  ut <- which(upper.tri(w) & !is.na(w) & w >= 0.3, arr.ind = TRUE);
  edf <- if (nrow(ut) > 0)
           data.frame(from = id[ut[, 1]], to = id[ut[, 2]],
                      weight = as.numeric(cw[ut]), stringsAsFactors = FALSE)
         else data.frame(from = character(0), to = character(0), weight = numeric(0));
  g <- graph_from_data_frame(edf, directed = FALSE,
                             vertices = data.frame(name = id, stringsAsFactors = FALSE));
  idx <- unlist(sapply(V(g)$name, function(x) which(x == id)));
  
  # define local function
  my.normalize <- function(x){
    return((x- min(x)) /(max(x)-min(x)))
  }
  my.rescale <- function(x, from, to){
    (x - min(x)) / max(x - min(x)) * (to - from) + from
  }
  
  
  # Compute normalized p-values for color gradient
  normalized_pvalues <- -log(my.normalize(pvalue) + min(pvalue/2))
  
  # Ensure that you compute colors only for existing vertices in the graph
  existing_vertices <- V(g)$name
  vertex_colors <- ComputeColorGradient(normalized_pvalues[names(normalized_pvalues) %in% existing_vertices], "black", F, F)
  vertex_colorsw <- ComputeColorGradient(normalized_pvalues[names(normalized_pvalues) %in% existing_vertices], "white", F, F)
  
  # Assign colors only to existing vertices
  V(g)$color <- vertex_colors
  V(g)$colorw <- vertex_colorsw
  
  cnt <- hits;
  names(cnt) <- id;
  cnt2 <- cnt[V(g)$name];
  
  if (all(cnt2 == cnt2[1])){
    V(g)$size <- rep(16, length(cnt2))
  }else{
    V(g)$size <- my.rescale(log(cnt2+1, base=10), 8, 24);
  }
  
  # layout using Fruchterman-Reingold with balanced parameters
  # Compact overall layout with good within-component separation
  n_nodes <- vcount(g);
  # Updated for igraph >= 0.8.0: removed deprecated parameters (area, repulserad, start.temp, grid)
  pos.xy <- layout_with_fr(g, niter=500)
  # tighten layout to reduce spacing between disconnected components
  pos.xy <- sweep(pos.xy, 2, colMeans(pos.xy), "-") * 0.6

  # Cache the plot inputs (graph, layout, counts, significance, labels) to disk
  # under netNm so a static enrichment-map image matching this interactive
  # network can be rendered later without recomputing the layout.
  tryCatch(saveRDS(list(g = g, layout = pos.xy, counts = cnt2,
                        signif = -log10(pmax(pvalue[V(g)$name], 1e-300)),
                        labels = V(g)$name),
                   paste0(netNm, "_plot_inputs.rds")),
           error = function(e) message("[enrichnet] caching plot inputs for ",
                                       netNm, " failed: ", conditionMessage(e)));

  # now create the json object
  nodes <- vector(mode="list");
  node.nms <- V(g)$name;
  node.sizes <- V(g)$size;
  node.cols <- V(g)$color;
  node.colsw <- V(g)$colorw;
  
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id = node.nms[i],
      label=node.nms[i],
      size = node.sizes[i],
      true_size=node.sizes[i], 
      molType="set",
      colorb=node.cols[i],
      colorw=node.colsw[i],
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    );
  }
  
  edge.mat <- as_edgelist(g);
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);
  
  # covert to json
  bedges <- stack(hits.query);
  b.mat <- matrix(NA, nrow=nrow(bedges), ncol=2);
  b.mat[,1] <- bedges[,"values"];
  b.mat[,2] <- as.character(bedges[,"ind"]);
  b.mat <- b.mat[complete.cases(b.mat),]
  colnames(b.mat) <- c("source", "target");
  bg <- graph_from_data_frame(b.mat, directed=F);
  idx <- unlist(sapply(V(bg)$name, function(x) which(x == id)));
  cols <- color_scale("red", "#E5C494");
  
  V(bg)$color[V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(-log(pvalue), "black", F, F);
  V(bg)$colorw[V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(-log(pvalue), "white", F, F);
  node.nms <- V(bg)$name;
  if(anal.type == "onedata"){
    tbl <- dataSet$comp.res
    tbl <- tbl[which(doEntrez2SymbolMapping(rownames(tbl), paramSet$data.org, paramSet$data.idType) %in% V(bg)$name),]
    expr.val <- tbl[,paramSet$selectedFactorInx];
    expvals <- expr.val;
    names(expvals) <- doEntrez2SymbolMapping(rownames(tbl), paramSet$data.org, paramSet$data.idType)
    expvals <- expvals[node.nms]
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T,T);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, T);
  }else if(anal.type == "genelist" && sum(as.numeric(paramSet$all.prot.mat[,1])) != 0){
    tbl <- paramSet$all.prot.mat
    gene.nms <- V(bg)$name[which(!V(bg)$name %in% rownames(enr.mat))]
    tbl <- tbl[which(tbl[,2] %in% gene.nms),]
    expr.val <- tbl[,1];
    names(expr.val) <- tbl[,2]
    expvals <- expr.val
    expvals <- expvals[node.nms]
    expvals <- expvals[!is.na(expvals)]
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, T);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, T);
    
  }else if(anal.type =="metadata"){
    if(paramSet$selDataNm == "meta_default"){
      tbl <- analSet$meta.mat.all
      tbl <- tbl[which(doEntrez2SymbolMapping(rownames(tbl), paramSet$data.org, paramSet$data.idType) %in% V(bg)$name),]
      expvals <- analSet$meta.avgFC[rownames(tbl)]
      names(expvals) <- doEntrez2SymbolMapping(rownames(tbl), paramSet$data.org, paramSet$data.idType)
      expvals <- expvals[node.nms]
      V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, T);
      V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, T);
    }else{
      dataSet <- readDataset(paramSet$selDataNm);
      tbl <- dataSet$comp.res;
      tbl <- tbl[which(doEntrez2SymbolMapping(rownames(tbl), paramSet$data.org, paramSet$data.idType) %in% V(bg)$name),]
      expvals <- tbl[,"logFC"];
      names(expvals) <- doEntrez2SymbolMapping(rownames(tbl), paramSet$data.org, paramSet$data.idType);
      expvals <- expvals[node.nms]
      expvals <- expvals[!is.na(expvals)]
      inx <- !V(bg)$name %in% rownames(enr.mat);
      V(bg)$color[inx] <- ComputeColorGradient(unname(expvals), "black", T, T);
      V(bg)$colorw[inx] <- ComputeColorGradient(unname(expvals), "black", T, T);
    }
  }else{
    expvals <- rep(0,length(V(bg)$color)); 
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
  }
  node.dgr2 <- as.numeric(degree(bg));
  V(bg)$size <- my.rescale(log(node.dgr2, base=10), 8, 24); 
  
  # layout using Fruchterman-Reingold with balanced parameters
  # Compact overall layout with good within-component separation
  n_nodes_bg <- vcount(bg);
  # Updated for igraph >= 0.8.0: removed deprecated parameters (area, repulserad, start.temp, grid)
  pos.xy <- layout_with_fr(bg, niter=500)
  # tighten layout to reduce spacing between disconnected components
  pos.xy <- sweep(pos.xy, 2, colMeans(pos.xy), "-") * 0.6
  
  # now create the json object
  bnodes <- vector(mode="list");
  node.sizes <- V(bg)$size;
  node.cols <- V(bg)$color;
  node.colsw <- V(bg)$colorw;
  
  shapes <- rep("circle", length(node.nms));
  hit.inx <- node.nms %in% b.mat[,"source"];
  shapes[hit.inx] <- "gene";
  node.lbls <- doEntrez2SymbolMapping(node.nms, paramSet$data.org, paramSet$data.idType)
  for(i in 1:length(node.sizes)){
    bnodes[[i]] <- list(
      id = node.nms[i],
      label=node.lbls[i], 
      size=node.sizes[i], 
      colorb=node.cols[i],
      colorw=node.colsw[i],
      true_size=node.sizes[i], 
      molType=shapes[i],
      exp= unname(expvals[node.nms[i]]),
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    );
  }
  
  ppi.comps <- vector(mode="list");
  paramSet$current.net.nm <- netNm
  ppi.comps[[netNm]] <- bg;
  analSet$ppi.comps <- ppi.comps
  
  bedge.mat <- as_edgelist(bg);
  bedge.mat <- cbind(id=paste0("b", 1:nrow(bedge.mat)), source=bedge.mat[,1], target=bedge.mat[,2]);
  initsbls <- doEntrez2SymbolMapping(analSet$list.genes, paramSet$data.org, paramSet$data.idType)
  names(initsbls) <- analSet$list.genes
  
  #for rjson generation
  edge.mat <- apply(edge.mat, 1, as.list)
  bedge.mat <- apply(bedge.mat, 1, as.list)
  enr.mat <- apply(enr.mat, 1, as.list)
  
  #paramSet$current.sigmajsR.nm <- paste0(netNm, ".rda");
  #save(nodes, edge.mat, file = paramSet$current.sigmajsR.nm);

  netData <- list(nodes=nodes, 
                  edges=edge.mat, 
                  bnodes=bnodes, 
                  bedges=bedge.mat, 
                  enr=unname(enr.mat), 
                  id=names(enr.mat), 
                  sizes=analSet$listSizes, 
                  hits=hits.query, 
                  genelist=initsbls, 
                  analType=anal.type, 
                  org=paramSet$data.org, 
                  backgroundColor=list("#514F6A", "#222222"),
                  dat.opt = paramSet$selDataNm,
                  naviString = "Enrichment Network");
  
  netName <- paste0(netNm, ".json");
  paramSet$partialToBeSaved <- c( paramSet$partialToBeSaved, c(netName));
  paramSet$jsonNms$network <- netName;
  saveSet(paramSet, "paramSet");
  
  analSet$enrichNet <- netData;
  saveSet(analSet, "analSet");
  sink(netName);
  cat(rjson::toJSON(netData));
  sink();
  return(analSet);
}
