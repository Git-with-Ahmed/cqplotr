script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
app_file <- file.path(script_dir, "app.R")
docs_dir <- file.path(script_dir, "docs")
build_dir <- file.path(script_dir, ".shinylive-build")

if (!file.exists(app_file)) {
  stop(sprintf("app.R not found at %s", app_file))
}

if (!requireNamespace("shinylive", quietly = TRUE)) {
  stop("Package 'shinylive' is not installed.")
}

unlink(docs_dir, recursive = TRUE, force = TRUE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)
unlink(build_dir, recursive = TRUE, force = TRUE)
dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
copied <- file.copy(app_file, file.path(build_dir, "app.R"), overwrite = TRUE)
if (!isTRUE(copied)) {
  stop("Failed to copy app.R into the Shinylive build directory.")
}

if (dir.exists(file.path(script_dir, "www"))) {
  dir.create(file.path(build_dir, "www"), recursive = TRUE, showWarnings = FALSE)
  file.copy(
    file.path(script_dir, "www"),
    build_dir,
    recursive = TRUE,
    overwrite = TRUE
  )
}

shinylive::export(
  appdir = build_dir,
  destdir = docs_dir
)

app_hash <- substr(unname(tools::md5sum(app_file)), 1, 12)
build_stamp <- list(
  source = "app.R",
  hash = app_hash,
  built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)

replace_once <- function(path, pattern, replacement, fixed = FALSE) {
  size <- file.info(path)$size
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  text <- rawToChar(readBin(con, what = "raw", n = size))
  matched <- grepl(pattern, text, perl = !fixed, fixed = fixed)
  if (!matched) {
    stop(sprintf("Pattern not found in %s", path))
  }
  updated <- gsub(pattern, replacement, text, perl = !fixed, fixed = fixed)
  out <- file(path, open = "wb")
  on.exit(close(out), add = TRUE)
  writeBin(charToRaw(updated), out)
}

replace_once(
  file.path(docs_dir, "index.html"),
  '\\./shinylive/load-shinylive-sw\\.js"',
  sprintf('./shinylive/load-shinylive-sw.js?v=%s"', app_hash)
)

replace_once(
  file.path(docs_dir, "index.html"),
  '\\./shinylive/shinylive\\.js"',
  sprintf('./shinylive/shinylive.js?v=%s"', app_hash)
)

replace_once(
  file.path(docs_dir, "shinylive", "shinylive.js"),
  'fetch("./app.json")',
  sprintf('fetch("./app.json?v=%s")', app_hash),
  fixed = TRUE
)

replace_once(
  file.path(docs_dir, "shinylive-sw.js"),
  'var version = "v9";',
  sprintf('var version = "v9-%s";', app_hash)
)

writeLines(
  c(
    "{",
    sprintf('  "source": "%s",', build_stamp$source),
    sprintf('  "hash": "%s",', build_stamp$hash),
    sprintf('  "built_at": "%s"', build_stamp$built_at),
    "}"
  ),
  file.path(docs_dir, "build.json"),
  useBytes = TRUE
)

message(sprintf("Shinylive build complete: %s", docs_dir))
message(sprintf("Build hash: %s", app_hash))

unlink(build_dir, recursive = TRUE, force = TRUE)
