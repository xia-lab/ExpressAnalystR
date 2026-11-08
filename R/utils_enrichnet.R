

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
  # ONE edge per pair; edge weight = number of shared genes so the static map can scale thickness.
  # the edges are BuildEnrichNetGraph's, the one computation every enrichment network of this package uses
  edf <- BuildEnrichNetGraph(geneSets[id], stats::setNames(seq_len(n), id), top = n)$edges[, c("from", "to", "weight")];
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

# THE EDGE RULE of the enrichment networks (owner, 7 Nov 2026): two sets are linked when their Jaccard index |A n B| / |A u B| is
# at least 0.25 (the standard cutoff is 0.25-0.375) and they share at least 2 members. Jaccard is the adaptive part: it scales with
# the sets' sizes and penalizes a size difference, so a small set inside a large one is not linked by a few members; the floor of 2
# only stops a single shared member from linking two tiny sets (metabolite and taxon sets, small gene hits). A count floor by data
# type (5 genes) removed real links between small pathways. It was an overlap ratio of 0.2-0.3 with no minimum. The same rule in
# every tool package, so every network drawn from it links pathways the same way.
enrich_net_edge_ok <- function(a, b, min_shared = 2, min_jaccard = 0.25) {
  a <- unique(as.character(unlist(a))); b <- unique(as.character(unlist(b)));
  sh <- length(intersect(a, b)); u <- length(union(a, b));
  sh >= min_shared && u > 0 && sh / u >= min_jaccard
}
# Keep the edges of g whose two sets pass enrich_net_edge_ok; `sets` is a list of members named by the vertex names.
enrich_net_keep_edges <- function(g, sets, min_shared = 2, min_jaccard = 0.25) {
  if (igraph::ecount(g) == 0) return(g);
  el <- igraph::as_edgelist(g, names = TRUE);
  ok <- vapply(seq_len(nrow(el)), function(k) enrich_net_edge_ok(sets[[el[k, 1]]], sets[[el[k, 2]]], min_shared, min_jaccard), logical(1));
  igraph::delete_edges(g, igraph::E(g)[!ok])
}

# THE ENRICHMENT NETWORK of any set of enriched pathways (7 Nov 2026): the computing, here in the tool package; a caller
# lays it out and draws it. `hits.query` = the members of each pathway hit by the query (a list named by pathway), `pvals` = a
# p-value per pathway (named). The `top` most significant pathways are kept (a pathway named twice keeps its first entry), and two
# are linked by enrich_net_edge_ok. Returns list(ids, members, hits, pvals, edges = data.frame(from, to, weight = shared members,
# jaccard)), or NULL when no pathway has a member. Base R only.
BuildEnrichNetGraph <- function(hits.query, pvals, top = 20, min_shared = 2, min_jaccard = 0.25) {
  if (is.null(hits.query) || !length(hits.query)) return(NULL);
  hq <- lapply(hits.query, function(x) { x <- unique(as.character(unlist(x))); x[!is.na(x) & nzchar(x)] });
  hq <- hq[lengths(hq) > 0];
  hq <- hq[!duplicated(names(hq))];
  if (!length(hq)) return(NULL);
  ids <- names(hq);
  pv <- suppressWarnings(as.numeric(pvals[ids])); pv[!is.finite(pv)] <- 1;
  ord <- order(pv);
  if (length(ord) > top) ord <- ord[seq_len(top)];
  ids <- ids[ord]; pv <- stats::setNames(pv[ord], ids); hq <- hq[ids];
  n <- length(ids);
  fr <- character(0); to <- character(0); wt <- numeric(0); jc <- numeric(0);
  if (n >= 2) for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    sh <- length(intersect(hq[[i]], hq[[j]]));
    if (sh > 0 && enrich_net_edge_ok(hq[[i]], hq[[j]], min_shared, min_jaccard)) {
      fr <- c(fr, ids[i]); to <- c(to, ids[j]); wt <- c(wt, sh);
      jc <- c(jc, sh / length(union(hq[[i]], hq[[j]])));
    }
  }
  list(ids = ids, members = hq, hits = lengths(hq), pvals = pv,
       edges = data.frame(from = fr, to = to, weight = wt, jaccard = jc, stringsAsFactors = FALSE))
}
