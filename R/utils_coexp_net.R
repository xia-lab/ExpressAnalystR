
# ====================================================================
# BuildIgraphFromCEM  —  make a data-derived gene–gene network
# --------------------------------------------------------------------
# file        : path to the qs-saved CEMiTool object (“cem.qs”)
# thresh      : numeric, keep edges with weight > thresh
# return      : igraph object with vertex/edge attributes
# ====================================================================
BuildIgraphFromCEM <- function(thresh    = 0.05,
                               layoutFun = igraph::layout_nicely) {
  
  ## ── 0 · packages ─────────────────────────────────────────────────
  library(CEMiTool)
  library(igraph)
  library(reshape2)
  
  ## ── 1 · read the CEMiTool object ────────────────────────────────
  cem <- qs::qread("cem.qs")
  
  ## ── 2 · make sure we have an adjacency matrix -------------------
  ##     (CEMiTool stores β only if you explicitly asked for it.)
  get_beta <- function(cem) {
    # 1) try stored β
    if (!is.null(cem@parameters$beta))
      return(as.numeric(cem@parameters$beta))
    
    # 2) try default scale-free heuristic (≥ v1.29)
    beta <- tryCatch({
      get_cemitool_r2_beta(cem)[2]   # returns c(R2, β)
    }, error = function(e) NA)
    
    # 3) last-resort: WGCNA pickSoftThreshold on the expression matrix
    if (is.na(beta)) {
      expr <- CEMiTool:::get_expression(cem)   # matrix genes × samples
      sft  <- WGCNA::pickSoftThreshold(t(expr), verbose = 0)
      beta <- sft$powerEstimate
      if (is.na(beta)) beta <- 6              # fallback default
    }
    beta
  }
  
  if (is.null(cem@adjacency)) {
    beta <- get_beta(cem)
    cem  <- get_adj(cem, beta = beta)
  }
  adj <- adj_data(cem)                        # square genes × genes
  
  ## ── 3 · build edge list above threshold -------------------------
  edge.df <- melt(adj)
  edge.df <- subset(edge.df, value > thresh & Var1 != Var2)
  
  ## keep at most 2 000 heaviest edges
  if (nrow(edge.df) > 2000) {
    edge.df <- edge.df[order(edge.df$value, decreasing = TRUE), ]
    edge.df <- edge.df[1:2000, ]
  }
  
  g       <- graph_from_data_frame(edge.df, directed = FALSE)
  E(g)$weight <- edge.df$value
  
  
  
  ## ── 4 · add vertex-level annotations ----------------------------
  mod.df <- cem@module                       # cols: genes, modules
  idx    <- match(V(g)$name, mod.df$genes)
    V(g)$module <- mod.df$modules[idx]

deg   <- igraph::degree(g)
ldeg  <- log10(deg + 1)                       # stabilise high degrees

# rescale helper
resc <- function(x) (x - min(x)) / (max(x) - min(x))

# colour ramp: yellow → dark red (works on white & black)
pal <- colorRampPalette(c("#FFD54F", "#FFA726", "#EF5350", "#B71C1C"))(10)

V(g)$color  <- pal[ ceiling( resc(ldeg) * 9 ) + 1 ]   # for light bg
V(g)$colorw <- V(g)$color                             # same for dark bg
  
  ## size by degree
  rescale <- function(x, from = 8, to = 20)
    (x - min(x)) / (max(x) - min(x)) * (to - from) + from
  V(g)$size <- rescale(log10(degree(g) + 1))
  #V(g)$size <- 8;

  ## ── 5 · 2-D layout coordinates ----------------------------------
  xy <- layoutFun(g)
  V(g)$posx <- xy[, 1]
  V(g)$posy <- xy[, 2]
  analSet <- qs::qread("analSet.qs");
  analSet$overall.graph <- g;
  analSet$overall.graph.orig <- g;

  saveSet(analSet, "analSet"):

  return(1)
}

CorrIgraph2SigmaJS <- function(g,
                               netNm     = "coexp_net",
                               paramSet,
                               analSet) {
  
  symVec <- doEntrez2SymbolMapping(V(g)$name,
                                   paramSet$data.org,
                                   paramSet$data.idType)
  
  nodes <- lapply(seq_len(vcount(g)), function(i) {
    v   <- V(g)[i]
    lbl <- if (!is.na(symVec[i]) && nzchar(symVec[i])) symVec[i] else v$name
    
    list(
      id        = as.character(v$name),   # still Entrez as key
      label     = lbl,                    # SYMBOL shown in SigmaJS
      size      = unclass(v$size)[1],
      true_size = unclass(v$size)[1],
      molType   = "gene",
      colorb    = as.character(v$color),
      colorw    = as.character(v$colorw),
      exp       = if (!is.null(v$expr)) unclass(v$expr)[1] else 0,
      posx      = unclass(v$posx)[1],
      posy      = unclass(v$posy)[1]
    )
  })
  
  
  ## ── edges with rescaled size 0.5–2.5  ----------------------------
  el <- igraph::as_data_frame(g, what = "edges")       # from, to, weight
  
  wMin <- min(el$weight)
  wMax <- max(el$weight)
  rescale <- function(x, from = 0.5, to = 2.5) {
    if (wMax == wMin) return((from + to) / 2)          # avoid 0/0
    (x - wMin) / (wMax - wMin) * (to - from) + from
  }
  
  edges <- lapply(seq_len(nrow(el)), function(i) {
    w <- as.numeric(el$weight[i])
    list(
      id     = paste0("e", i),
      source = as.character(el$from[i]),
      target = as.character(el$to[i]),
      weight = w,                     # keeps the raw weight
      size   = rescale(w)             # stroke width 0.5–2.5
    )
  })
  
  
  dataSet <- readDataset(paramSet$dataName);
  nodeTable <- BuildNodeTable(g,paramSet,dataSet,analSet);
  
  initsbls <- doEntrez2SymbolMapping( V(g)$name, paramSet$data.org, paramSet$data.idType)
  names(initsbls) <- V(g)$name
  ppi.net <- list();

  ppi.net[["node.data"]] <- data.frame(Id=V(g)$name, Label=unname(initsbls));
  ppi.net <<- ppi.net;
  ## ── 3 · assemble JSON payload -----------------------------------
  netData <- list(nodes            = nodes,
                  edges            = edges,
                  backgroundColor  = list("#f5f5f5", "#0066CC"),
                  naviString       = "Correlation Network",
                  org              = paramSet$data.org,
                  genelist=initsbls, 
                  nodeTable = nodeTable);
  
  fileNm <- paste0(netNm, ".json")
  jsonlite::write_json(netData, fileNm, auto_unbox = TRUE)
  
  ## track for later download
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, fileNm)
  paramSet$jsonNms$coexpNet <- basename(fileNm)
  saveSet(paramSet, "paramSet")
  analSet$corNet <- netData
  saveSet(analSet, "analSet")
  
  invisible(netData)
}

SplitIgraphByModule <- function(g, keepXTalk = FALSE) {

  stopifnot("module" %in% vertex_attr_names(g))
  
  mods <- sort(unique(V(g)$module))
  subG <- setNames(vector("list", length(mods)), mods)
  
  for (m in mods) {
    genes.m <- V(g)[module == m]
    
    if (keepXTalk) {
      # keep all edges touching those genes
      subG[[m]] <- induced_subgraph(g, vids = genes.m)
    } else {
      # keep *only* edges whose BOTH endpoints are in module m
      eKeep <- E(g)[inc(V(g)[module == m])]
      eKeep <- eKeep[ which(ends(g, eKeep)[,1] %in% genes.m$name &
                              ends(g, eKeep)[,2] %in% genes.m$name) ]
      subG[[m]] <- subgraph.edges(g, eKeep)
    }
  }
  subG
}

GenerateCEMModuleNetworks <- function(fileName  = "coexp_network",
                                      thresh    = 0.05,
                                      keepXTalk = FALSE,
                                      minNodeNum = 3) {

  paramSet <- readSet(paramSet, "paramSet")
  analSet  <- readSet(analSet,  "analSet")
  library(igraph);

  print(names(analSet));
  g.all   <- analSet$overall.graph.orig
  g.byMod <- SplitIgraphByModule(g.all, keepXTalk = keepXTalk)
  
  comps <- g.byMod


  # first compute subnet stats
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("module", 1:length(comps), sep="");

  # note, we report stats for all nets (at least 3 nodes);
  hit.inx <- net.stats$Node >= minNodeNum;
  comps <- comps[hit.inx];


  # now record
  net.stats <<- net.stats;
  sub.stats <- unlist(lapply(comps, vcount));

  
  ## ── 3 · write JSON for *only the first* module in the list ─────
  firstMod <- names(g.byMod)[1]              # e.g. "M1"
  netNm <- fileName;
  analSet$ppi.comps <- comps;
  saveSet(analSet, "analSet")
  return(c(vcount(g.all), ecount(g.all), length(comps), sub.stats));
}


BuildNodeTable <- function(g,
                           paramSet,
                           dataSet = NULL,
                           analSet = NULL) {

  ids    <- V(g)$name
  labels <- doEntrez2SymbolMapping(ids,
                                   paramSet$data.org,
                                   paramSet$data.idType)
  anal.type <- paramSet$anal.type
  ## --- topological features --------------------------------------
  deg  <- igraph::degree(g)
  btw  <- igraph::betweenness(g)

  ## --- expression values (log2FC) --------------------------------
  expr <- rep(0, length(ids))         # default NA

  if (anal.type == "onedata") {

    tbl <- dataSet$comp.res
    inx <- match(ids, rownames(tbl))
    expr <- tbl[inx, paramSet$selectedFactorInx]

  } else if (anal.type == "metadata") {

    if (paramSet$selDataNm == "meta_default") {
      tbl  <- analSet$meta.mat.all
      sy   <- doEntrez2SymbolMapping(rownames(tbl),
                                     paramSet$data.org,
                                     paramSet$data.idType)
      inx  <- match(ids, sy)
      expr <- analSet$meta.avgFC[rownames(tbl)][inx]

    } else {                     # metadata but user-selected dataset
      ds   <- readDataset(paramSet$selDataNm)
      tbl  <- ds$comp.res
      sy   <- doEntrez2SymbolMapping(rownames(tbl),
                                     paramSet$data.org,
                                     paramSet$data.idType)
      inx  <- match(ids, sy)
      expr <- tbl[inx, "logFC"]
    }
  }
  ## for "genelist" expr stays NA  (per requirement)

  ## --- assemble data-frame ---------------------------------------
  node.df <- data.frame(
    id          = ids,
    label       = ifelse(is.na(labels) | labels == "", ids, labels),
    degree      = deg,
    betweenness = btw,
    expr        = expr,
    stringsAsFactors = FALSE
  )

  node.df <- node.df[order(node.df$degree, decreasing = TRUE), ]

  node.df
}


GetNetsName <- function(){
  rownames(net.stats);
}

GetNetsNameString <- function(){
  paste(rownames(net.stats), collapse="||");
}

GetNetsEdgeNum <- function(){
  as.numeric(net.stats$Edge);
}

GetNetsNodeNum <- function(){
  as.numeric(net.stats$Node);
}

GetNetsQueryNum <- function(){
  as.numeric(net.stats$Query);
}


ComputeSubnetStats <- function(comps){
  library(igraph);
  net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
  colnames(net.stats) <- c("Node", "Edge", "Query");
  for(i in 1:length(comps)){
    g <- comps[[i]];
    net.stats[i,] <- c(vcount(g),ecount(g),vcount(g));
  }
  return(net.stats);
}
# ====================================================================
# filterNetByThresh  —  edge-weight filtering for an igraph object
# --------------------------------------------------------------------
# g          : igraph object that already has a numeric edge attribute
#              called 'weight'
# thresh     : keep edges with weight > thresh
# maxEdges   : cap the network at this many heaviest edges (NULL = no cap)
# rmIsolated : TRUE → delete nodes that become isolated after filtering
# layoutFun  : (optional) layout recalculation if you need updated coords
# return     : list(graph = <filtered igraph>,
#                   stats = c(nodes, edges, n.components))
# ====================================================================
FilterNetByThresh <- function(thresh      = 0.05,
                                 maxEdges    = 2000,
                                 rmIsolated  = TRUE) {
  # save.image("filter.RDAta");
  analSet  <- readSet(analSet,  "analSet")
  overall.graph <- analSet$overall.graph;
  g <- overall.graph;
  if (!"weight" %in% edge_attr_names(g))
    stop("edge attribute 'weight' not found")

  # ── 1 · keep only edges above threshold ───────────────────────────
  g <- subgraph.edges(g, E(g)[weight > thresh], delete.vertices = FALSE)

  # ── 2 · cap total edges if requested ──────────────────────────────
  if (!is.null(maxEdges) && ecount(g) > maxEdges) {
    el <- igraph::as_data_frame(g, what = "edges")
    el <- el[order(el$weight, decreasing = TRUE), ][seq_len(maxEdges), ]
    g  <- graph_from_data_frame(el,
                                directed = FALSE,
                                vertices = igraph::as_data_frame(g, what = "vertices"))
    E(g)$weight <- el$weight                    # restore weights
  }

  # ── 3 · optionally drop newly-isolated nodes ──────────────────────
  if (rmIsolated) {
    iso <- which(igraph::degree(g) == 0)
    if (length(iso) > 0)
      g <- delete_vertices(g, iso)
  }

  # ── 4 · final tidy-up (loops / multi-edges) ───────────────────────
  g <- simplify(g, edge.attr.comb = list("first"))

  # ── 6 · book-keeping & return  ────────────────────────────────────
  current.msg <<-
    paste("FilterNetByThreshold:",
          vcount(g), "nodes and", ecount(g), "edges retained at thresh >", thresh);
  

    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
  outStats <- c(vcount(g), ecount(g), length(substats))
  analSet$overall.graph <- g;
    return(saveSet(analSet, "analSet", outStats));
}


FilterBipartiNet <- function(nd.type, min.dgr, min.btw){
    paramSet <- readSet(paramSet, "paramSet");
    analSet <- readSet(analSet, "analSet");
    overall.graph <- analSet$overall.graph
    all.nms <- V(overall.graph)$name;
    edge.mat <- get.edgelist(overall.graph);
    dgrs <- degree(overall.graph);
    nodes2rm.dgr <- nodes2rm.btw <- NULL;

    if(nd.type == "gene"){
        hit.inx <- all.nms %in% edge.mat[,1];
    }else if(nd.type=="other"){
        hit.inx <- all.nms %in% edge.mat[,2];
    }else{ # all
        hit.inx <- rep(TRUE, length(all.nms));
    }

    if(min.dgr > 0){
        rm.inx <- dgrs <= min.dgr & hit.inx;
        nodes2rm.dgr <- V(overall.graph)$name[rm.inx];
    }
    if(min.btw > 0){
        btws <- betweenness(overall.graph);
        rm.inx <- btws <= min.btw & hit.inx;
        nodes2rm.btw <- V(overall.graph)$name[rm.inx];
    }

    nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
    overall.graph <- simplify(delete.vertices(overall.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
    if(!is.null(substats)){
        output <- c(vcount(overall.graph), ecount(overall.graph), length(analSet$ppi.comps), substats);
    }else{
        output <- 0;
    }
    analSet$overall.graph <- overall.graph;

    return(saveSet(analSet, "analSet", output));
}

PrepareNetwork <- function(net.nm, jsonNm){
   analSet <- readSet(analSet, "analSet");
   paramSet <- readSet(paramSet, "paramSet");
   #print(analSet$ppi.comps);

   my.ppi <- analSet$ppi.comps[[net.nm]];
   nd.nms <- V(my.ppi)$name;
  print(names(analSet))
  CorrIgraph2SigmaJS(my.ppi,
                     netNm    = jsonNm,
                     paramSet = paramSet,
                     analSet  = analSet)

   current.net.nm <<- net.nm;
   paramSet$current.net.nm <- net.nm;
   saveSet(paramSet, "paramSet");

   return(saveSet(analSet, "analSet", 1));

}

PerformNetEnrichment <- function(dataName="", file.nm, fun.type, IDs){
  #dataName <<- dataName;
  #file.nm <<- file.nm;
  #fun.type <<- fun.type;
  #IDs <<- IDs;
  #save.image("PerformNetEnrichment.RData");
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  data.org <- paramSet$data.org;
  # prepare query
  ora.vec <- NULL;
  idtype <- "entrez";
  
    ora.vec <- unlist(strsplit(IDs, "; "));
    names(ora.vec) <- as.character(ora.vec);
  
  
  if(fun.type %in% c("trrust", "encode", "jaspar", "mirnet", "met", "drugbank", "disease")){
    res <- PerformRegEnrichAnalysis(dataSet, file.nm, fun.type, ora.vec, "inverse");
  }else{
    res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, ora.vec, "coexp");
  }

  return(res);
}

PerformRegEnrichAnalysis <- function(dataSet, file.nm, fun.type, ora.vec, netInv){
    if(!exists("my.reg.enrich")){ # public web on same user dir
        compiler::loadcmp(paste0(resource.dir, "rscripts/ExpressAnalystR/R/_utils_regenrich.Rc"));    
    }
    return(my.reg.enrich(dataSet, file.nm, fun.type, ora.vec, netInv));
}


FindCommunities <- function(method = "walktrap",
                            use.weight = FALSE,
                            component = "largest",#c("largest", "all", "seed"),
                            min_genes = 5,
                            require_query_hit = TRUE) {
  # save.image("find.Rdata");
  paramSet    <- readSet(paramSet, "paramSet")
  analSet    <- readSet(analSet, "analSet")

  seed.expr   <- paramSet$seed.expr
  ppi.comps <- analSet$ppi.comps
  current.net <- ppi.comps[[current.net.nm]]
  
  if (igraph::vcount(current.net) < 2L) return("NA||Graph too small")
  
  # ---- choose component(s) ----------------------------------------------------
  pick_largest <- function(g) {
    comp <- igraph::components(g)
    igraph::induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  }
  
  graph_list <- switch(
    component,
    "largest" = list(pick_largest(current.net)),
    "all" = {
      comp <- igraph::components(current.net)
      split(seq_len(igraph::vcount(current.net)), comp$membership) |>
        lapply(function(vs) igraph::induced_subgraph(current.net, vs))
    },
    "seed" = {
      seeds <- intersect(seed.proteins, igraph::V(current.net)$name)
      if (length(seeds) == 0) return("NA||No seeds in network!")
      comp <- igraph::components(current.net)
      # keep components that contain at least one seed
      keep_comp <- unique(comp$membership[match(seeds, igraph::V(current.net)$name)])
      vs <- which(comp$membership %in% keep_comp)
      list(igraph::induced_subgraph(current.net, vs))
    }
  )
  
  # ---- helper: run detection on one graph ------------------------------------
  run_one <- function(g) {
    if (!igraph::is_connected(g)) g <- pick_largest(g)
    if (igraph::vcount(g) < 2L) return(list(vec = character(0), tbl = NULL))
    
    # -- FIX 1: ensure vertex names exist; fallback to other attrs/indices -----
    vnames <- igraph::V(g)$name
    if (is.null(vnames)) vnames <- rep(NA_character_, igraph::vcount(g))
    if (all(is.na(vnames))) {
      alt <- NULL
      if (!is.null(igraph::V(g)$Id))    alt <- as.character(igraph::V(g)$Id)
      if (!is.null(igraph::V(g)$Label)) alt <- as.character(igraph::V(g)$Label)
      vnames <- if (!is.null(alt)) alt else as.character(seq_len(igraph::vcount(g)))
      igraph::V(g)$name <- vnames
    }
    
    # symbol mapping with fallback to name
    hit.x <- match(vnames, ppi.net$node.data[, 1])
    sybls <- ppi.net$node.data[hit.x, 2]
    sybls[is.na(sybls)] <- vnames[is.na(sybls)]
    names(sybls) <- vnames
    
    # -- FIX 2: robust weight assignment (no invalid indexing) -----------------
    weights <- igraph::E(g)$weight
    # community detection (weights used when supported)
    fc <- switch(
      method,
      "walktrap"  = igraph::cluster_walktrap(g, weights = weights),
      "infomap"   = igraph::cluster_infomap(g, e.weights = weights),
      "labelprop" = igraph::cluster_label_prop(g, weights = weights),
      { return(list(err = "NA||Unknown method!")) }
    )
    
    comm_list <- igraph::communities(fc)
    if (length(comm_list) == 0 || igraph::modularity(fc, weights = weights) == 0) {
      return(list(vec = character(0), tbl = NULL))
    }
    
    # iterate communities
    community.vec  <- character(0)
    gene.community <- NULL
    rowcount <- 0L
    
    for (i in seq_along(comm_list)) {
      vids       <- comm_list[[i]]     # vertex indices
      comm_size  <- length(vids)
      if (comm_size < min_genes) next
      
      comm_names  <- vnames[vids]
      # ensure no NA names in printable path
      comm_names[is.na(comm_names)] <- as.character(vids[is.na(comm_names)])
      
      # label strings (symbols if available)
      comm_labels <- sybls[comm_names]
      comm_labels[is.na(comm_labels)] <- comm_names
      
      # query hits
      qnums <- comm_size

      if (require_query_hit && qnums == 0) next
      
      # in/out degree test
      subg    <- igraph::induced_subgraph(g, vids)
      in.deg  <- igraph::degree(subg)
      out.deg <- igraph::degree(g, vids) - in.deg
      ppval   <- suppressWarnings(
        wilcox.test(in.deg, out.deg, exact = FALSE, alternative = "two.sided")$p.value
      )
      ppval   <- signif(ppval, 3)
      
      # record
      rowcount <- rowcount + 1L
      pids     <- paste(comm_labels, collapse = "->")
      
      com.mat <- cbind(Id     = comm_names,
                       Label  = comm_labels,
                       Module = as.character(i))  # keep as before
      gene.community <- rbind(gene.community, com.mat)
      
      community.vec[rowcount] <- paste(c(comm_size, qnums, ppval, pids), collapse = ";")
    }
    
    if (length(community.vec) == 0) return(list(vec = character(0), tbl = NULL))
    
    # order: size desc, p-value asc (same as before)
    community_data <- do.call(rbind, lapply(community.vec, function(x) {
      parts <- strsplit(x, ";")[[1]]
      data.frame(size = as.numeric(parts[1]), p_value = as.numeric(parts[3]))
    }))
    ord <- with(community_data, order(-size, p_value))
    
    list(vec = community.vec[ord], tbl = gene.community)
  }
  # run and combine
  out_list <- lapply(graph_list, run_one)
  if (any(vapply(out_list, function(x) !is.null(x$err), logical(1)))) {
    return(out_list[[which(vapply(out_list, function(x) !is.null(x$err), logical(1)))] ]$err)
  }
  
  vecs <- unlist(lapply(out_list, `[[`, "vec"), use.names = FALSE)
  tbls <- do.call(rbind, lapply(out_list, `[[`, "tbl"))
  
  if (length(vecs) == 0) return("NA||No communities were detected!")
  
  all.communities <- paste(vecs, collapse = "||")
  if (!is.null(tbls) && nrow(tbls) > 0) {
    colnames(tbls) <- c("Id", "Label", "Module")
    fast.write(tbls, file = "module_table.csv", row.names = FALSE)
    return(all.communities)
  } else {
    return("NA")
  }
}


community.significance.test <- function(graph, vs, ...) {
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}


# from to should be valid nodeIDs
GetShortestPaths <- function(from, to){
  current.net <- ppi.comps[[current.net.nm]];
  paths <- get.all.shortest.paths(current.net, from, to)$res;
  if(length(paths) == 0){
    return (paste("No connection between the two nodes!"));
  }
  
  path.vec <- vector(mode="character", length=length(paths));
  for(i in 1:length(paths)){
    path.inx <- paths[[i]]; 
    path.ids <- V(current.net)$name[path.inx];
    path.sybls <- path.ids;
    pids <- paste(path.ids, collapse="->");
    psbls <- paste(path.sybls, collapse="->");
    path.vec[i] <- paste(c(pids, psbls), collapse=";")
  }
  
  if(length(path.vec) > 50){
    path.vec <- path.vec[1:50];
  }
  
  all.paths <- paste(path.vec, collapse="||");
  return(all.paths);
}

DecomposeGraph <- function(gObj,analSet, minNodeNum = 3, jsonBool = F){
  # now decompose to individual connected subnetworks
    if(jsonBool == "netjson"){
        comps <-list(gObj)
    }else{
        comps <-decompose.graph(gObj, min.vertices=minNodeNum);
    }
  if(length(comps) == 0){
    current.msg <<- paste("No subnetwork was identified with at least", minNodeNum, "nodes!");
    return(NULL);
  }
  
  # first compute subnet stats
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("module", 1:length(comps), sep="");
  
  # note, we report stats for all nets (at least 3 nodes);
  hit.inx <- net.stats$Node >= minNodeNum;
  comps <- comps[hit.inx];
  #overall <- list();
  #overall[["overall"]] <- g
  #ppi.comps <- append(overall, ppi.comps, after=1);
  
  # now record
  net.stats <<- net.stats;
  sub.stats <- unlist(lapply(comps, vcount)); 
  analSet$ppi.comps <- comps;
  analSet$net.stats <- net.stats;
  analSet$substats <- sub.stats;
  return(analSet);
}