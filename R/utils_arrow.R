#' Arrow Shadow Save Utility Functions
#'
#' This file provides helper functions for the Arrow migration.
#' The shadow_save() function writes data in both .qs and .arrow formats
#' to enable gradual migration from qs to Apache Arrow.
#'
#' Special feature: Preserves R matrix/data.frame rownames by converting
#' them to a column named 'row_names_id' in the Arrow file.
#'
#' @author ExpressAnalyst Team

#' Shadow save: writes both .qs and .arrow formats
#'
#' This function saves data in the original .qs format and additionally
#' creates an .arrow (Feather) file for Java-side reading. The .qs file
#' maintains backward compatibility while the .arrow file enables
#' direct Java reading without RServe.
#'
#' R matrix rownames are preserved by converting them to a column named
#' 'row_names_id' before writing to Feather format.
#'
#' @param obj The R object to save (data.frame, matrix, or list)
#' @param file The file path (should end with .qs)
#' @return Invisible NULL. Called for side effects.
#' @export
shadow_save <- function(obj, file) {
  # Original qs save (keep for backward compatibility)
  qs::qsave(obj, file)

  # Arrow shadow save
  arrow_path <- sub("\\.qs$", ".arrow", file)
  if (!grepl("\\.arrow$", arrow_path)) {
    arrow_path <- paste0(file, ".arrow")
  }

  # Convert to data.frame if possible, preserving rownames
  tryCatch({
    if (is.data.frame(obj) || is.matrix(obj)) {
      df <- preserve_rownames(obj)
      arrow::write_feather(df, arrow_path, compression = "uncompressed")
    } else if (is.list(obj) && !inherits(obj, "phyloseq") && !is.null(names(obj))) {
      # For simple lists with data.frame components, save each component
      for (nm in names(obj)) {
        if (is.data.frame(obj[[nm]]) || is.matrix(obj[[nm]])) {
          component_path <- sub("\\.arrow$", paste0("_", nm, ".arrow"), arrow_path)
          df <- preserve_rownames(obj[[nm]])
          arrow::write_feather(df, component_path, compression = "uncompressed")
        }
      }
    }
    # Skip phyloseq and other complex S4 objects for now
  }, error = function(e) {
    warning(paste("Arrow shadow save failed for", file, ":", e$message))
  })

  invisible(NULL)
}

#' Preserve rownames by converting to a column
#'
#' Converts R matrix/data.frame rownames to a column named 'row_names_id'.
#' This allows Java Arrow readers to access the original R rownames.
#'
#' @param obj A matrix or data.frame
#' @return A data.frame with row_names_id as the first column
#' @keywords internal
preserve_rownames <- function(obj) {
  df <- as.data.frame(obj)

  # Check if rownames are meaningful (not just default 1, 2, 3...)
  rn <- rownames(df)
  if (!is.null(rn) && length(rn) > 0) {
    # Check if rownames are not just sequential integers
    default_names <- as.character(seq_len(nrow(df)))
    if (!identical(rn, default_names)) {
      # Prepend row_names_id column
      df <- cbind(row_names_id = rn, df)
    }
  }

  return(df)
}

#' Direct Arrow save (no qs, for new code paths)
#'
#' Saves data directly to Arrow format without qs backup.
#' Use this for new code paths where qs compatibility is not needed.
#'
#' @param obj The R object to save (data.frame or matrix)
#' @param file The file path (should end with .arrow)
#' @return Invisible NULL. Called for side effects.
#' @export
arrow_save <- function(obj, file) {
  if (!grepl("\\.arrow$", file)) {
    file <- paste0(file, ".arrow")
  }

  tryCatch({
    if (is.data.frame(obj) || is.matrix(obj)) {
      df <- preserve_rownames(obj)
      arrow::write_feather(df, file, compression = "uncompressed")
    } else {
      stop("arrow_save only supports data.frame and matrix objects")
    }
  }, error = function(e) {
    stop(paste("Arrow save failed for", file, ":", e$message))
  })

  invisible(NULL)
}
