# Internal helpers -----------------------------------------------------
make_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d_%H%M")
}

.check_project <- function(path) {
  if (
    !dir.exists(file.path(path, "output")) ||
      !dir.exists(file.path(path, "data"))
  ) {
    stop("Project folders not found. Run initialiseProject() first.")
  }
}
