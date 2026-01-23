# Helper function to suppress ALL graphics device popups on macOS
.suppress_quartz <- function(expr) {
  # Save current state
  old_device <- getOption("device")
  old_bitmapType <- getOption("bitmapType")

  # Override device options to prevent Quartz
  options(device = function(...) pdf(file = NULL))
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    options(bitmapType = "cairo")
  }

  # Close any stray devices
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)

  # Ensure cleanup happens
  on.exit({
    while (dev.cur() > 1) try(dev.off(), silent = TRUE)
    options(device = old_device, bitmapType = old_bitmapType)
  })

  # Execute the expression
  eval(expr, envir = parent.frame())
}

.trim_cem_object_for_save <- function(cem) {
  # Reset slots that are not used downstream so the saved object stays small.
  slot(cem, "fit_indices") <- data.frame()
  slot(cem, "selected_genes") <- character(0)
  slot(cem, "enrichment") <- list()
  slot(cem, "ora") <- data.frame()
  slot(cem, "interactions") <- list()
  slot(cem, "interaction_plot") <- list()
  slot(cem, "profile_plot") <- list()
  slot(cem, "enrichment_plot") <- list()
  slot(cem, "mean_k_plot") <- list()
  slot(cem, "barplot_ora") <- list()
  slot(cem, "mean_var_plot") <- ggplot2::ggplot()
  slot(cem, "hist_plot") <- ggplot2::ggplot()
  slot(cem, "qq_plot") <- ggplot2::ggplot()
  slot(cem, "sample_tree_plot") <- gtable::gtable(
    widths = grid::unit(0, "cm"),
    heights = grid::unit(0, "cm"))
  slot(cem, "mod_colors") <- character(0)
  slot(cem, "input_params") <- list()
  slot(cem, "calls") <- list()
  slot(cem, "adjacency") <- matrix(numeric(0), 0, 0)
  cem
}

my.build.cemi.net <- function(dataName,
                              filter      = TRUE,
                              min_ngen    = 30,
                              cor_method  = "pearson",
                              verbose     = FALSE,
                              classCol    = NULL) {   # <-- optional argument
  tryCatch({

    ## 1 · load dataset -------------------------------------------------
    dataSet  <- readDataset(dataName)
    expr_mat <- as.data.frame(dataSet$data.norm)    # genes × samples

    ## metadata: keep *all* columns, coerce factors -> character
    meta_df <- dataSet$meta.info
    meta_df[] <- lapply(meta_df, \(x) if (is.factor(x)) as.character(x) else x)

    ## decide which column is the class
    if (is.null(classCol)) {
      classCol <- colnames(meta_df)[1]              # first column by default
    }
    if (!classCol %in% colnames(meta_df)) {
      stop("classCol '", classCol, "' not found in meta.info")
    }

    ## build annotation table (SampleName + all meta)
    annot_df <- data.frame(SampleName = rownames(meta_df),
                           meta_df,
                           check.names = FALSE,
                           stringsAsFactors = FALSE)

    ## 2 · run CEMiTool -------------------------------------------------
    suppressPackageStartupMessages({
      library(CEMiTool)
      library(WGCNA)
    })

print("buildingceminet");

    # FIX: Suppress Quartz popup on macOS - completely disable plotting during cemitool
    # We'll generate plots separately using the other functions
    cem <- cemitool(expr              = expr_mat,
                    annot             = annot_df,
                    filter            = filter,
                    min_ngen          = min_ngen,
                    cor_method        = match.arg(cor_method),
                    class_column      = classCol,
                    verbose           = verbose,
                    plot              = FALSE,           # Disable all plotting
                    plot_diagnostics  = FALSE)

    ## 3 · save & return -----------------------------------------------
    cem <- .trim_cem_object_for_save(cem)
    qs::qsave(cem, "cem.qs")

  mod <- attr(cem, "module")

   if (is.null(mod) || !is.data.frame(mod) || nrow(mod) == 0 || !("modules" %in% colnames(mod))) {
      if(nrow(cem@sample_annotation) < 25){
      return("ERROR: No modules found. Beta selection likely failed during network construction due to low number of samples. Please have look at the Scale-free fit plot.");
      }else{
      return(paste0("ERROR: No modules found. Beta selection likely failed during network construction. ",
           "Consider relaxing filtering, lowering minimum module size or ",
           "increasing sample size. Please have look at the Scale-free fit plot."));
      }
    }else{
      return("OK")
    }
  }, error = function(e) {
    return(paste("Error:", conditionMessage(e)))
  })
}


PlotCEMiDendro <- function(mode      = c("sample", "module"),
                           metaClass = "NA",
                           imgName   = "cem_dendro",
                           dpi       = 72,
                           format    = "png") {

  library(Cairo); library(WGCNA)

  # FIX: Suppress Quartz popup on macOS - override device at function start
  old_device <- getOption("device")
  old_bitmapType <- getOption("bitmapType")
  options(device = function(...) pdf(file = NULL))
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    options(bitmapType = "cairo")
  }
  on.exit({
    options(device = old_device, bitmapType = old_bitmapType)
  }, add = TRUE)

  cem <- qs::qread("cem.qs")
  if (!inherits(cem, "CEMiTool"))
    stop("'cem.qs' does not contain a valid CEMiTool object.")

  mode <- match.arg(mode)
  expr <- attr(cem, "expression")       # genes × samples
  mod <- attr(cem, "module")

  ## helper ----------------------------------------------------------
  plotDendroColoured <- function(hc, colNamed, label, file, legendPal) {

    leaves <- labels(as.dendrogram(hc))
    if (!all(leaves %in% names(colNamed)))
      stop("Colour vector missing leaves:\n  ",
           paste(setdiff(leaves, names(colNamed)), collapse = ", "))

    colMat <- matrix(colNamed[leaves], nrow = length(leaves), ncol = 1,
                     dimnames = list(leaves, NULL))

    if (dpi == 72) dpi <- 96
    width_in <- 10
    height_in <- 6

    # FIX: Suppress Quartz popup on macOS - close any existing devices first
    while (dev.cur() > 1) dev.off()

    Cairo(file, width = width_in, height = height_in, dpi = dpi,
          bg = "white", type = format, units = "in")

    oldMar <- par("mar"); par(mar = oldMar + c(0, 0, 0, 4))
    plotDendroAndColors(
      dendro       = hc,
      colors       = colMat,
      groupLabels  = label,
      dendroLabels = FALSE,
      addGuide     = TRUE,
      guideHang    = 0.05,
      main         = paste("Module dendrogram"),
      ylab         = "1 − Pearson correlation")

    ## legend
    par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), xpd = NA)
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("topright",
           legend = names(legendPal),
           fill   = legendPal,
           border = NA,
           cex    = 0.8,
           bty    = "n")

    # FIX: Suppress Quartz popup on macOS
    par(oldMar)
    invisible(dev.off())
  }

  # ── MODULE dendrogram ────────────────────────────────────────────
  if (mode == "module") {
    ME  <- moduleEigengenes(t(expr[mod$genes, ]),
                            colors = mod$modules)$eigengenes
    hc  <- hclust(as.dist(1 - cor(ME)), method = "average")

    ids   <- colnames(ME)
    pal   <- setNames(rainbow(length(ids)), ids)
    idVec <- setNames(ids, ids)

    file <- sprintf("%s_module_dendro_dpi%d.%s", imgName, dpi, format)
    plotDendroColoured(hc, idVec, "modules", file, pal)
    return(1)
  }

  # ── SAMPLE dendrogram ────────────────────────────────────────────
  sa <- cem@sample_annotation

  ## choose metadata column sensibly
  if (is.na(metaClass) || metaClass == "NA") metaClass <- 2  # 1 is SampleName
  if (is.numeric(metaClass)) {
    if (metaClass < 1 || metaClass > ncol(sa))
      stop("'metaClass' index out of range.")
    if (metaClass == 1)
      stop("metaClass 1 is 'SampleName'; choose a metadata column (>=2).")
    classes <- sa[[metaClass]]
  } else {
    if (!metaClass %in% colnames(sa))
      stop("metaClass '", metaClass, "' not found in sample_annotation.")
    if (metaClass == "SampleName")
      stop("metaClass 'SampleName' is invalid; choose real metadata.")
    classes <- sa[[metaClass]]
  }

  names(classes) <- sa$SampleName

  hc  <- hclust(as.dist(1 - cor(expr)), method = "average")
  pal <- setNames(rainbow(length(unique(classes))), unique(classes))
  colNamed <- setNames(pal[classes], names(classes))

  file <- sprintf("%sdpi%d.%s", imgName, dpi, format)
  plotDendroColoured(hc, colNamed, "sample class", file, pal)
    imgSet <- readSet(imgSet, "imgSet");
    imgSet$coexp_dendrogram <- file;
    saveSet(imgSet, "imgSet"):
  return("OK")
}

# =======================================================================
# PlotCEMiTreatmentHeatmap
# -----------------------------------------------------------------------
# factorName : name (character) or index (numeric) of the categorical
#              column in cem@sample_annotation to expand into dummies
# imgName    : file stem for the output
# dpi        : resolution
# format     : "png" | "pdf"
# -----------------------------------------------------------------------
# returns 1 on success, 0 on failure
# =======================================================================
PlotCEMiTreatmentHeatmap <- function(factorName,
                                     imgName = "cem_treatment_heatmap",
                                     dpi     = 96,
                                     format  = c("png", "pdf")) {

  tryCatch({
    dbg_file <- "coexp_plot_debug.log"
    try(writeLines(paste0("[PlotCEMiTreatmentHeatmap] wd=", getwd(),
                          " imgName=", imgName, " dpi=", dpi,
                          " format=", paste(format, collapse = ",")),
                   con = dbg_file, sep = "\n", useBytes = TRUE),
        silent = TRUE)

    library(CEMiTool); library(WGCNA); library(Cairo)

    # FIX: Suppress Quartz popup on macOS - override device at function start
    old_device <- getOption("device")
    old_bitmapType <- getOption("bitmapType")
    options(device = function(...) pdf(file = NULL))
    if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
      options(bitmapType = "cairo")
    }
    on.exit({
      options(device = old_device, bitmapType = old_bitmapType)
    }, add = TRUE)

    cem <- qs::qread("cem.qs")
    stopifnot(inherits(cem, "CEMiTool"))

    sa <- cem@sample_annotation
    ## ── 1 · validate factorName  ----------------------------------
    if (is.numeric(factorName)) {
      stopifnot(factorName >= 2, factorName <= ncol(sa))
      fac      <- sa[[factorName]]
      colLabel <- colnames(sa)[factorName]
    } else {
      stopifnot(factorName %in% colnames(sa), factorName != "SampleName")
      fac      <- sa[[factorName]]
      colLabel <- factorName
    }
    fac <- as.factor(fac)

    ## ── 2 · dummy matrix  -----------------------------------------
    mm <- model.matrix(~ 0 + fac)
    colnames(mm) <- levels(fac)
    rownames(mm) <- sa$SampleName

    ## ── 3 · module eigengenes  ------------------------------------
    expr   <- attr(cem, "expression")
    modTbl <- attr(cem, "module")
    g      <- intersect(rownames(expr), modTbl$genes)
    colors <- modTbl$modules[match(g, modTbl$genes)]
    MEs    <- moduleEigengenes(t(expr[g, , drop = FALSE]), colors)$eigengenes

    mm <- mm[rownames(MEs), , drop = FALSE]

    ## ── 4 · correlations  -----------------------------------------
    corMat <- cor(MEs, mm, use = "p")
    pMat   <- corPvalueStudent(corMat, nSamples = nrow(MEs))
    textMat <- paste0(formatC(corMat, 2), "\n(",
                      formatC(pMat , 1, format = "e"), ")")

    ## ── 5 · device  (min 8 × 6 in)  -------------------------------
    colfun <- colorRampPalette(c("royalblue4", "white", "tomato"))

    outFile <- sprintf("%sdpi%d.%s", imgName, dpi, match.arg(format))
    if(dpi == 72){
        dpi <- 96;
    }
    width_px  <- max(50 + 40 * ncol(corMat), 8 * dpi)
    height_px <- max(50 + 20 * nrow(corMat), 6 * dpi)
    width_in  <- width_px  / dpi
    height_in <- height_px / dpi

    # FIX: Suppress Quartz popup on macOS - close any existing devices first
    while (dev.cur() > 1) dev.off()

    if (tolower(format) == "png") {
      Cairo(file   = outFile,
            width  = width_in,
            height = height_in,
            dpi    = dpi,
            bg     = "white",
            type   = "png",
            units  = "in")
    } else {
      Cairo(file   = outFile,
            width  = width_in,
            height = height_in,
            bg     = "white",
            type   = "pdf")
    }

    ## optional: tighten default margins a little
    oldMar <- par("mar"); on.exit(par(oldMar), add = TRUE)
    par(mar = c(5, 9, 4, 2) + 0.1)

    labeledHeatmap(Matrix          = corMat,
                   xLabels         = colnames(mm),
                   yLabels         = rownames(corMat),
                   ySymbols        = rownames(corMat),
                   colorLabels     = FALSE,
                   colors          = colfun(50),
                   textMatrix      = textMat,
                   setStdMargins   = FALSE,
                   cex.text        = 0.7,
                   zlim            = c(-1, 1),
                   main            = paste("Module ×", colLabel, "Levels"))

    # FIX: Suppress Quartz popup on macOS
    invisible(dev.off())
    #message("Heat-map written to: ", outFile)
    try(writeLines(paste0("[PlotCEMiTreatmentHeatmap] outFile=", outFile,
                          " exists=", file.exists(outFile)),
                   con = dbg_file, sep = "\n", useBytes = TRUE, append = TRUE),
        silent = TRUE)
    imgSet <- readSet(imgSet, "imgSet");
    imgSet$coexp_traitheat <- outFile;
    saveSet(imgSet, "imgSet");
    1

  }, error = function(e) {
    try(writeLines(paste0("[PlotCEMiTreatmentHeatmap] error=", conditionMessage(e)),
                   con = dbg_file, sep = "\n", useBytes = TRUE, append = TRUE),
        silent = TRUE)
    #message("PlotCEMiTreatmentHeatmap: ", e$message)
    0
  })
}

PlotCemiScaleFree <- function(imgName = "coexp_scalefree",
                                     dpi = 72,
                                     format = "png") {
  library(Cairo); library(CEMiTool)

  # FIX: Suppress Quartz popup on macOS - override device at function start
  old_device <- getOption("device")
  old_bitmapType <- getOption("bitmapType")
  options(device = function(...) pdf(file = NULL))
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    options(bitmapType = "cairo")
  }
  on.exit({
    options(device = old_device, bitmapType = old_bitmapType)
  }, add = TRUE)

  cem <- qs::qread("cem.qs")
  stopifnot(inherits(cem, "CEMiTool"))

  # Ensure the plot exists (some versions only populate it after calling plot_beta)
  if (is.null(cem@beta_r2_plot)) {
    plot_beta <- try(getFromNamespace("plot_beta", "CEMiTool"), silent = TRUE)
    if (!inherits(plot_beta, "try-error") && is.function(plot_beta)) {
      cem <- plot_beta(cem)  # fills cem@beta_r2_plot
    }
  }

  # Extract ggplot object from the slot
  g <- NULL
  if (!is.null(cem@beta_r2_plot)) {
    if (is.list(cem@beta_r2_plot) && "beta_r2_plot" %in% names(cem@beta_r2_plot)) {
      g <- cem@beta_r2_plot$beta_r2_plot
    } else {
      g <- cem@beta_r2_plot
    }
  }
  if (is.null(g)) return("Error: beta_r2_plot not available on the CEMiTool object.")

  # Save
  file <- sprintf("%sdpi%d.%s", imgName, dpi, format)
    if (dpi == 72) dpi <- 96
  width_in <- 10
  height_in <- 6

    # FIX: Suppress Quartz popup on macOS - close any existing devices first
    while (dev.cur() > 1) dev.off()

    Cairo(file, width = width_in, height = height_in, dpi = dpi, bg = "white", type = format, units = "in")
  invisible(print(g))    # ggplot draw
  invisible(dev.off())
    imgSet <- readSet(imgSet, "imgSet");
    imgSet$coexp_scalefree <- file;
    saveSet(imgSet, "imgSet");
  return(1);
}
