#' Initialise a maudr project folder structure
#'
#' Creates a self-contained top-level project folder (default: "maudr_assignments")
#' containing the subfolders `data/`, `code/`, and `output/` (with
#' `assignments_output/` and `answers_output/` inside).
#' Works in both RStudio and Positron; no .Rproj file is created.
#'
#' @param path Base directory in which to create the project folder (default ".").
#' @param name Name of the top-level project folder (default "maudr_assignments").
#' @param verbose Logical; show progress messages (default TRUE).
#'
#' @return The full path to the newly created project folder (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#' initialiseProject()
#' initialiseProject(name = "enzyme_project", path = "~/Documents")
#' }
initialiseProject <- function(
  path = ".",
  name = "maudr_assignments",
  verbose = TRUE
) {
  # Define top-level project folder
  project_path <- file.path(path, name)

  # Define subfolders
  dirs <- file.path(
    project_path,
    c("data", "output/assignments_output", "output/answers_output")
  )

  # Create directories recursively
  purrr::walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

  # # Add a small README with metadata
  # readme_path <- file.path(project_path, "README.txt")
  # if (!file.exists(readme_path)) {
  #   info <- c(
  #     "# maudr project",
  #     paste("Created:", Sys.time()),
  #     paste("maudr package version:", utils::packageVersion("maudr"))
  #   )
  #   writeLines(info, readme_path)
  # }

  # Optional message
  if (verbose) {
    message("Project folder initialised at ", normalizePath(project_path))
  }

  invisible(project_path)
}
