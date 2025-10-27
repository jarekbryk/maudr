#' Generate assignment Excel files for each student
#'
#' Reads the setup metadata (from [assignReactions()]) and produces per-student
#' Excel files with simulated absorbance vs. time data. Each student receives
#' both inhibited and uninhibited conditions in one sheet. The true inhibition
#' type remains in the instructor metadata.
#'
#' Random measurement noise ("jitter") is optionally applied per student in a
#' reproducible way based on the run seed and student ID.
#'
#' @param run_timestamp Timestamp string (YYYY-MM-DD_HHMM) from [assignReactions()].
#' @param project_path Path to the top-level project folder (default ".").
#' @param verbose Logical; show progress messages (default TRUE).
#' @param cuv_vol Cuvette volume in litres (default 0.003).
#' @param eps Extinction coefficient (default 6220).
#' @param enz_vol Enzyme volume in millilitres (default 0.1).
#' @param substr_conc Vector of substrate concentrations in mM.
#' @param time_vec Vector of reaction time points in minutes.
#' @param use_jitter Logical; add random noise to absorbances (default TRUE).
#' @param jitter_sd Numeric; SD of multiplicative noise if `use_jitter = TRUE`
#'   (default 0.02, i.e. ±2% variability).
#'
#' @return Invisibly returns TRUE after writing all student assignment files
#'   and updating metadata.
#' @export
generateAssignments <- function(
  run_timestamp,
  project_path = ".",
  verbose = TRUE,
  cuv_vol = 0.003,
  eps = 6220,
  enz_vol = 0.1,
  substr_conc = c(0, 10, 20, 40, 80, 160),
  time_vec = c(0.17, 0.33, 0.5, 0.66, 0.83, 1),
  use_jitter = TRUE,
  jitter_sd = 0.02
) {
  # --- Locate metadata -------------------------------------------------------
  assignment_dir <- file.path(
    project_path,
    "output",
    "assignments_output",
    run_timestamp
  )
  meta_file <- file.path(assignment_dir, "setup_metadata.xlsx")
  if (!file.exists(meta_file)) {
    stop("Metadata file not found for timestamp: ", run_timestamp)
  }
  meta <- readxl::read_excel(meta_file)

  # --- Helper functions ------------------------------------------------------
  calculateV <- function(Vmax, Km, substrate_conc) {
    Vmax * substrate_conc / (Km + substrate_conc)
  }

  calculateGradient <- function(Vmax, Km, substrate_conc) {
    V <- calculateV(Vmax, Km, substrate_conc)
    V * (enz_vol / 1e6) / cuv_vol * eps
  }

  # internal helper: generate Excel for one student
  generateAbsVsTimeData <- function(student_id, df, output_path) {
    df |>
      dplyr::mutate(
        student_id = student_id,
        substrate_conc_char = paste0(substrate_conc, "_mM")
      ) |>
      dplyr::select(-inhibition_actual, -substrate_conc) |>
      tidyr::pivot_wider(
        names_from = substrate_conc_char,
        values_from = absorbance
      ) |>
      dplyr::select(
        student_id,
        rxn_substrate,
        rxn_condition,
        rxn_time,
        tidyselect::everything()
      ) |>
      writexl::write_xlsx(file.path(
        output_path,
        paste0(student_id, "_data.xlsx")
      ))
  }

  # --- Simulate absorbance data ---------------------------------------------
  # Prepare full dataset (both conditions)
  all_data <- meta |>
    dplyr::mutate(
      substrate_conc = list(substr_conc),
      rxn_time = list(time_vec),
      rxn_condition = dplyr::case_when(
        inhibition_actual == "no_inhibition" ~ "without_inhibitor",
        TRUE ~ "with_inhibitor"
      )
    ) |>
    tidyr::unnest(substrate_conc) |>
    tidyr::unnest(rxn_time) |>
    dplyr::mutate(
      gradient = calculateGradient(Vmax, Km, substrate_conc),
      absorbance = gradient * rxn_time
    )

  # --- Apply deterministic jitter per student -------------------------------
  if (use_jitter) {
    # Recover seed from metadata for reproducibility
    seed <- unique(meta$seed)
    if (length(seed) != 1 || is.na(seed)) {
      seed <- 1234
    } # fallback

    all_data <- all_data |>
      dplyr::group_by(student_id) |>
      dplyr::mutate(
        # Set seed to enable individual students' jitter recovery
        seed_offset = match(student_id, unique(all_data$student_id)),
        absorbance = purrr::map2_dbl(
          absorbance,
          seed_offset,
          ~ {
            set.seed(seed + .y)
            .x * rnorm(1, 1, jitter_sd)
          }
        )
      ) |>
      dplyr::ungroup()
  }

  all_data <- all_data |>
    dplyr::mutate(absorbance = round(absorbance, 3)) |>
    dplyr::select(
      student_id,
      rxn_substrate,
      inhibition_actual,
      rxn_condition,
      substrate_conc,
      rxn_time,
      absorbance
    )

  # --- Write individual student files ---------------------------------------
  all_data |>
    dplyr::group_by(student_id) |>
    tidyr::nest() |>
    dplyr::mutate(output_path = assignment_dir) |>
    purrr::pwalk(~ generateAbsVsTimeData(..1, ..2, ..3))

  # --- Update metadata -------------------------------------------------------
  meta2 <- meta |>
    dplyr::mutate(
      assignment_file = file.path(
        assignment_dir,
        paste0(student_id, "_data.xlsx")
      ),
      timestamp = run_timestamp
    )
  writexl::write_xlsx(meta2, meta_file)

  if (verbose) {
    message("Updated metadata with assignment filenames: ", meta_file)
  }

  invisible(TRUE)
}
