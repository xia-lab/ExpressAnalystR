##################################################
## qs/qs2 compatibility shims
## Loaded first so the rest of the package can use these helpers.
##################################################

.expressanalyst_qsave <- function(object, file, ...) {
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qs_save(object, file = file)
  } else if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(object, file = file, preset = "fast")
  } else {
    stop("Neither qs2 nor qs is available")
  }
}

.expressanalyst_qread <- function(file, ...) {
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qs_read(file)
  } else if (requireNamespace("qs", quietly = TRUE)) {
    qs::qread(file)
  } else {
    stop("Neither qs2 nor qs is available")
  }
}

.expressanalyst_qd_read <- function(raw_or_file, ...) {
  if (is.raw(raw_or_file)) {
    tmp <- tempfile(pattern = "load_qs_", fileext = ".qs", tmpdir = getwd())
    on.exit(unlink(tmp), add = TRUE)
    writeBin(raw_or_file, tmp)
    if (requireNamespace("qs2", quietly = TRUE)) {
      qs2::qd_read(tmp, ...)
    } else if (requireNamespace("qs", quietly = TRUE)) {
      qs::qdeserialize(raw_or_file)
    } else {
      stop("Neither qs2 nor qs is available")
    }
  } else {
    if (requireNamespace("qs2", quietly = TRUE)) {
      qs2::qd_read(raw_or_file, ...)
    } else if (requireNamespace("qs", quietly = TRUE)) {
      qs::qdeserialize(raw_or_file)
    } else {
      stop("Neither qs2 nor qs is available")
    }
  }
}
