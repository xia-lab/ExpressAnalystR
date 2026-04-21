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
                              classCol    = NULL,
                              force_beta  = NULL) {   # <-- optional argument to force beta value
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

    # PRACTICAL LIMIT: Cap at top 5000 features by IQR to ensure reasonable computation time
    # Co-expression analysis is computationally expensive (O(n²) for correlation matrix)
    MAX_FEATURES <- 5000
    if (nrow(expr_mat) > MAX_FEATURES) {
      message("Dataset contains ", nrow(expr_mat), " features, selecting top ", MAX_FEATURES, " by IQR...");

      # Compute IQR for each feature
      feature_iqr <- apply(expr_mat, 1, IQR, na.rm = TRUE)

      # Select top features by IQR
      top_features <- order(feature_iqr, decreasing = TRUE)[1:MAX_FEATURES]
      expr_mat <- expr_mat[top_features, , drop = FALSE]

      message("Selected top ", MAX_FEATURES, " features by IQR for co-expression analysis");
    }

    ## 2 · run CEMiTool -------------------------------------------------
    suppressPackageStartupMessages({
      library(CEMiTool)
      library(WGCNA)
    })

    # Log package versions for debugging
    message("=== Environment Info ===")
    message("CEMiTool version: ", as.character(packageVersion("CEMiTool")))
    message("WGCNA version: ", as.character(packageVersion("WGCNA")))
    message("R version: ", R.version.string)

print("buildingceminet");

    # Log data dimensions and quality for debugging
    message("=== Data Info ===")
    message("Data dimensions: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")
    message("Number of samples: ", nrow(annot_df))
    message("Filter: ", filter, ", min_ngen: ", min_ngen)
    message("Class column: ", classCol)

    # Data quality checks
    message("Data range: [", min(expr_mat, na.rm = TRUE), ", ", max(expr_mat, na.rm = TRUE), "]")
    message("NA values: ", sum(is.na(expr_mat)))
    message("Infinite values: ", sum(is.infinite(as.matrix(expr_mat))))
    message("Zero variance genes: ", sum(apply(expr_mat, 1, var, na.rm = TRUE) == 0))

    # Check class column
    class_vals <- annot_df[[classCol]]
    message("Class column '", classCol, "' unique values: ", paste(unique(class_vals), collapse = ", "))
    message("Class counts: ", paste(table(class_vals), collapse = ", "))

    # FIX: Suppress Quartz popup on macOS - completely disable plotting during cemitool
    # We'll generate plots separately using the other functions
    message("=== Running CEMiTool ===")

    # Capture warnings during cemitool execution
    # If force_beta is specified, use it; otherwise let CEMiTool auto-select
    if (!is.null(force_beta)) {
      message("Using forced beta value: ", force_beta)
      cem <- withCallingHandlers(
        cemitool(expr              = expr_mat,
                 annot             = annot_df,
                 filter            = filter,
                 min_ngen          = min_ngen,
                 cor_method        = match.arg(cor_method),
                 class_column      = classCol,
                 verbose           = verbose,
                 plot              = FALSE,           # Disable all plotting
                 plot_diagnostics  = TRUE,            # Enable diagnostics to debug beta selection
                 force_beta        = force_beta),     # Force specific beta value
        warning = function(w) {
          message("WARNING during cemitool: ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    } else {
      message("Auto-selecting beta value")
      cem <- withCallingHandlers(
        cemitool(expr              = expr_mat,
                 annot             = annot_df,
                 filter            = filter,
                 min_ngen          = min_ngen,
                 cor_method        = match.arg(cor_method),
                 class_column      = classCol,
                 verbose           = verbose,
                 plot              = FALSE,           # Disable all plotting
                 plot_diagnostics  = TRUE),           # Enable diagnostics to debug beta selection
        warning = function(w) {
          message("WARNING during cemitool: ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    }

    ## 3 · save & return -----------------------------------------------

    # Log diagnostic information before trimming
    message("CEMiTool analysis completed")

    # Check if fit_indices exist and extract beta information
    if (!is.null(cem@fit_indices) && nrow(cem@fit_indices) > 0) {
      message("Fit indices available: ", nrow(cem@fit_indices), " rows")
      message("Beta values tested: ", paste(head(cem@fit_indices$beta, 10), collapse = ", "))

      # Try to find selected beta from fit_indices
      if ("selected" %in% colnames(cem@fit_indices)) {
        selected_row <- cem@fit_indices[cem@fit_indices$selected == TRUE, ]
        if (nrow(selected_row) > 0) {
          message("Selected beta: ", selected_row$beta[1])
        } else {
          message("WARNING: No beta marked as selected in fit_indices")
        }
      } else if ("beta" %in% colnames(cem@fit_indices)) {
        message("Beta range: ", min(cem@fit_indices$beta), " to ", max(cem@fit_indices$beta))
      }
    } else {
      message("WARNING: No fit indices found - beta selection may have failed")
    }

    # Check parameters slot if it exists
    if (length(slotNames(cem)) > 0 && "parameters" %in% slotNames(cem)) {
      params <- slot(cem, "parameters")
      if (is.list(params) && "beta" %in% names(params)) {
        message("Beta from parameters: ", params$beta)
      }
    }

    mod <- attr(cem, "module")
    if (!is.null(mod) && is.data.frame(mod) && nrow(mod) > 0) {
      message("Modules found: ", length(unique(mod$modules)))
      message("Genes in modules: ", nrow(mod))
    } else {
      message("ERROR: No modules detected")

      # If auto-selection failed and no forced beta was used, try with a reasonable beta value
      if (is.null(force_beta) && !is.null(cem@fit_indices) && nrow(cem@fit_indices) > 0) {
        message("=== Attempting retry with best available beta ===")

        # Find beta with highest SFT.R.sq that's >= 0.75 (slightly relaxed threshold)
        fit_df <- cem@fit_indices
        rsq_col <- NULL

        # Check for both possible column names
        if ("SFT.R.sq" %in% colnames(fit_df)) {
          rsq_col <- "SFT.R.sq"
        } else if ("R.sq" %in% colnames(fit_df)) {
          rsq_col <- "R.sq"
        }

        if (!is.null(rsq_col)) {
          message("Using R-squared column: ", rsq_col)
          good_fits <- fit_df[fit_df[[rsq_col]] >= 0.75, ]

          if (nrow(good_fits) > 0) {
            # Choose lowest beta among good fits (more conservative, larger modules)
            best_idx <- which.min(good_fits$Power)
            best_beta <- good_fits$Power[best_idx]
            best_rsq <- good_fits[[rsq_col]][best_idx]

            message("Retrying with beta = ", best_beta, " (", rsq_col, " = ", best_rsq, ")")

            # Recursive call with forced beta
            return(my.build.cemi.net(dataName, filter, min_ngen, cor_method, verbose, classCol, force_beta = best_beta))
          } else {
            message("No beta values with ", rsq_col, " >= 0.75 found")
            message("Available R.sq values: ", paste(fit_df[[rsq_col]], collapse = ", "))
          }
        } else {
          message("Could not find R.sq column in fit_indices")
          message("Available columns: ", paste(colnames(fit_df), collapse = ", "))
        }
      }
    }

    cem <- .trim_cem_object_for_save(cem)
    ov_qs_save(cem, "cem.qs")

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
    error_msg <- paste("Error:", conditionMessage(e))
    message("=== ERROR DETAILS ===")
    message(error_msg)
    message("Stack trace:")
    message(paste(capture.output(traceback()), collapse = "\n"))

    # Try to get more context
    if (exists("expr_mat")) {
      message("expr_mat exists: dim = ", paste(dim(expr_mat), collapse = " x "))
    }
    if (exists("annot_df")) {
      message("annot_df exists: dim = ", paste(dim(annot_df), collapse = " x "))
    }
    if (exists("cem")) {
      message("cem object exists: class = ", class(cem))
      if (inherits(cem, "CEMiTool")) {
        message("cem slots: ", paste(slotNames(cem), collapse = ", "))
      }
    }

    return(error_msg)
  })
}


PlotCEMiDendro <- function(mode      = c("sample", "module"),
                           metaClass = "NA",
                           imgName   = "cem_dendro",
                           dpi = default.dpi,
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

  cem <- ov_qs_read("cem.qs")
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
    # Detect feature column name (genes vs features)
    feature_col_dendro <- if ("genes" %in% colnames(mod)) "genes" else "features";

    ME  <- moduleEigengenes(t(expr[mod[[feature_col_dendro]], ]),
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

    cem <- ov_qs_read("cem.qs")
    stopifnot(inherits(cem, "CEMiTool"))

    sa <- cem@sample_annotation
    ## ── 1 · validate factorName  ----------------------------------
    # Default: use first metadata column (column 2, since column 1 is SampleName)
    if (is.na(factorName) || is.null(factorName) || factorName == "NA" || factorName == "") {
      message("factorName is NA/NULL/empty, using first metadata column");
      if (ncol(sa) < 2) {
        stop("Sample annotation has no metadata columns (only SampleName)");
      }
      colLabel <- colnames(sa)[2]
      fac <- sa[[2]]
      message("Auto-selected column: '", colLabel, "'");
    } else if (is.numeric(factorName)) {
      message("factorName is numeric: ", factorName);
      if (factorName < 2 || factorName > ncol(sa)) {
        stop("Numeric factorName (", factorName, ") is out of range. Must be between 2 and ", ncol(sa),
             ". Available columns: ", paste(colnames(sa), collapse=", "));
      }
      fac      <- sa[[factorName]]
      colLabel <- colnames(sa)[factorName]
    } else {
      message("factorName is character: '", factorName, "'");

      # Validate column exists - try case-insensitive match
      available_cols <- colnames(sa)

      if (factorName %in% available_cols) {
        # Exact match
        fac      <- sa[[factorName]]
        colLabel <- factorName
      } else {
        # Try case-insensitive match
        matched_col <- available_cols[tolower(available_cols) == tolower(factorName)]

        if (length(matched_col) > 0) {
          message("Found case-insensitive match: '", matched_col[1], "' for requested '", factorName, "'");
          fac <- sa[[matched_col[1]]]
          colLabel <- matched_col[1]
        } else {
          # No match - provide helpful error
          stop("factorName '", factorName, "' not found in sample annotation (case-insensitive search failed). ",
               "Available columns: ", paste(available_cols, collapse=", "), ". ",
               "Please check that your metadata column name matches exactly.");
        }
      }

      if (colLabel == "SampleName") {
        stop("factorName cannot be 'SampleName'. Available columns: ",
               paste(available_cols[available_cols != "SampleName"], collapse=", "));
      }
    }
    fac <- as.factor(fac)

    ## ── 2 · dummy matrix  -----------------------------------------
    mm <- model.matrix(~ 0 + fac)
    colnames(mm) <- levels(fac)
    rownames(mm) <- sa$SampleName

    ## ── 3 · module eigengenes  ------------------------------------
    expr   <- attr(cem, "expression")
    modTbl <- attr(cem, "module")

    # CRITICAL FIX: CEMiTool uses different column names in different versions
    feature_col <- if ("genes" %in% colnames(modTbl)) {
      "genes"
    } else if ("features" %in% colnames(modTbl)) {
      "features"
    } else {
      stop("Module table has no 'genes' or 'features' column. Available: ",
           paste(colnames(modTbl), collapse=", "));
    }

    g      <- intersect(rownames(expr), modTbl[[feature_col]])
    colors <- modTbl$modules[match(g, modTbl[[feature_col]])]
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
            unit   = "in",
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
                                     dpi = default.dpi,
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

  cem <- ov_qs_read("cem.qs")
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
