# Standalone check for ExpressAnalystR.
#
# The package is loaded by the web application as plain script files, so a
# function that is defined somewhere else in the deployment but not here would
# only fail at run time, for the user. This script makes that a build-time
# failure: every R/*.R file must parse, every top-level function must be
# definable, and every call to a shared helper (ov_*, .ov_*, .ai_*) must resolve
# to a function defined in this package. No package dependencies are needed.
#
#   Rscript tools/check_standalone.R
#
# Exit status is non-zero on any problem.

r.dir <- if (dir.exists("R")) "R" else file.path(dirname(sys.frame(1)$ofile), "..", "R")
files <- list.files(r.dir, pattern = "[.]R$", full.names = TRUE)
if (length(files) == 0) stop("no R files found under ", r.dir)

env <- new.env(parent = globalenv())
problems <- character(0)

# 1. parse and define
for (f in files) {
  exprs <- tryCatch(parse(f, keep.source = FALSE), error = function(e) {
    problems <<- c(problems, sprintf("%s: parse error: %s", basename(f), conditionMessage(e)))
    NULL
  })
  if (is.null(exprs)) next
  for (e in exprs) {
    is.fun.def <- is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
    if (is.fun.def) {
      tryCatch(eval(e, env), error = function(err) {
        problems <<- c(problems, sprintf("%s: cannot define %s: %s", basename(f),
                                         as.character(e[[2]]), conditionMessage(err)))
      })
    }
  }
}
cat(sprintf("parsed %d files, defined %d functions\n", length(files), length(ls(env, all.names = TRUE))))

# 2. shared helpers must be defined here
defined <- ls(env, all.names = TRUE)
shared.pattern <- "^[.]?(ov|ai)_"
missing <- list()
for (fn in defined) {
  obj <- get(fn, env)
  if (!is.function(obj)) next
  used <- tryCatch(codetools::findGlobals(obj, merge = FALSE)$functions, error = function(e) character(0))
  used <- used[grepl(shared.pattern, used) & !(used %in% defined)]
  for (u in used) missing[[u]] <- c(missing[[u]], fn)
}
for (u in names(missing)) {
  problems <- c(problems, sprintf("undefined shared helper %s() called from: %s",
                                  u, paste(unique(missing[[u]]), collapse = ", ")))
}

if (length(problems)) {
  cat("FAIL\n"); cat(paste0("  - ", problems, "\n"), sep = "")
  quit(status = 1)
}
cat("OK: ExpressAnalystR is self-contained\n")
