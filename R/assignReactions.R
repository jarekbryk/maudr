#' Assign reaction conditions to students
#'
#' Randomly assigns each student one substrate and one inhibition type
#' (excluding "no_inhibition"), then creates a paired condition for
#' "no_inhibition". Each student therefore has two rows:
#' one inhibited and one uninhibited reaction. The function writes
#' a metadata Excel file to the timestamped `assignments_output` folder.
#'
#' If no input files are provided, demo datasets bundled with the package
#' are used automatically.
#'
#' @param student_file Optional path to Excel file containing student list
#'   with columns `student_no`, `first_name`, and `surname`.
#' @param enzyme_file Optional path to Excel file containing enzyme
#'   parameters (e.g. `reaction_parameters.xlsx`).
#' @param project_path Path to the top-level project folder created by
#'   [initialiseProject()]. Default is current working directory.
#' @param seed Random seed for reproducibility (default = 1234).
#' @param verbose Logical; if TRUE (default), prints progress messages.
#'
#' @return Invisibly returns a list with:
#'   \describe{
#'     \item{timestamp}{Timestamp string identifying the output folder.}
#'     \item{metadata}{A tibble with all student–reaction assignments.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' proj <- initialiseProject()
#' assignReactions(project_path = proj, seed = 4321)
#' }
assignReactions <- function(
  student_file = NULL,
  enzyme_file = NULL,
  project_path = ".",
  seed = 1234,
  verbose = TRUE
) {
  # Ensure project folder structure exists
  # initialiseProject(path = project_path, verbose = FALSE)

  # Create timestamped subfolder for this run
  timestamp <- make_timestamp()
  assignment_dir <- file.path(
    project_path,
    "output",
    "assignments_output",
    timestamp
  )
  dir.create(assignment_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Load student file -----------------------------------------------------
  if (is.null(student_file)) {
    demo_student <- system.file(
      "extdata",
      "student_names.xlsx",
      package = "maudr"
    )
    if (demo_student == "") {
      stop("Demo student file missing in package inst/extdata.")
    }
    students <- readxl::read_excel(demo_student)
    if (verbose) message("Using bundled demo student file.")
  } else {
    students <- readxl::read_excel(student_file)
    file.copy(
      student_file,
      file.path(project_path, "data", basename(student_file)),
      overwrite = TRUE
    )
    if (verbose) message("Copied student file to data/ folder.")
  }

  # --- Load enzyme parameters file ------------------------------------------
  if (is.null(enzyme_file)) {
    demo_enzyme <- system.file(
      "extdata",
      "reaction_parameters.xlsx",
      package = "maudr"
    )
    if (demo_enzyme == "") {
      stop("Demo enzyme file missing in package inst/extdata.")
    }
    enzymes <- readxl::read_excel(demo_enzyme)
    if (verbose) message("Using bundled demo enzyme parameters.")
  } else {
    enzymes <- readxl::read_excel(enzyme_file)
    file.copy(
      enzyme_file,
      file.path(project_path, "data", basename(enzyme_file)),
      overwrite = TRUE
    )
    if (verbose) message("Copied enzyme parameters file to data/ folder.")
  }

  # --- Prepare student IDs ---------------------------------------------------
  students <- students |>
    dplyr::mutate(
      student_id = stringr::str_trim(
        toupper(paste(student_no, first_name, sep = "_"))
      )
    )

  # --- Random assignment of reaction conditions ----------------------------------------------
  set.seed(seed)

  rxn_substr <- unique(enzymes$rxn_substrate)
  inh_types <- setdiff(unique(enzymes$inhibition_actual), "no_inhibition")

  # Each student gets one substrate + one inhibition type
  metadata <- students |>
    dplyr::mutate(
      rxn_substrate = sample(rxn_substr, dplyr::n(), replace = TRUE),
      inhibition_actual = sample(inh_types, dplyr::n(), replace = TRUE)
    ) |>
    tidyr::uncount(2, .id = "id") |>
    dplyr::mutate(
      inhibition_actual = ifelse(id == 2, "no_inhibition", inhibition_actual)
    ) |>
    dplyr::select(-student_no, -first_name, -surname, -id) |>
    dplyr::left_join(enzymes, by = c("rxn_substrate", "inhibition_actual"))

  # --- Add run metadata ------------------------------------------------------
  meta_out <- metadata |>
    dplyr::mutate(seed = seed, timestamp = timestamp)

  # --- Write metadata Excel file --------------------------------------------
  meta_path <- file.path(assignment_dir, "setup_metadata.xlsx")
  writexl::write_xlsx(meta_out, meta_path)

  if (verbose) {
    message("Reaction assignments saved to: ", normalizePath(meta_path))
  }

  # Return summary invisibly
  invisible(list(timestamp = timestamp, metadata = meta_out))
}
