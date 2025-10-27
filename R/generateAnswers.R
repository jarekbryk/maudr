#' Generate answer PDFs from assignment Excel files
#'
#' Reads student assignment files, estimates rates and (optionally) Km/Vmax,
#' and renders per-student answer pages with:
#'  - Absorbance vs time (faceted by rxn_condition)
#'  - Michaelis-Menten (both conditions on one plot; optional nls fit)
#'  - Lineweaver-Burke (both conditions on one plot; fitted lm lines)
#'
#' @param run_timestamp Timestamp string (YYYY-MM-DD_HHMM) from assignReactions().
#' @param output_files One of "separate", "single", or "both". Default "single".
#' @param project_path Path to the top-level project folder (default ".").
#' @param verbose Logical; show progress messages (default TRUE).
#' @return Invisibly TRUE.
#' @export
generateAnswers <- function(
  run_timestamp,
  output_files = c("single", "separate", "both"),
  project_path = ".",
  verbose = TRUE
) {
  output_files <- match.arg(output_files)

  # --- Paths & metadata ------------------------------------------------------
  assignment_dir <- file.path(
    project_path,
    "output",
    "assignments_output",
    run_timestamp
  )
  answers_dir <- file.path(
    project_path,
    "output",
    "answers_output",
    run_timestamp
  )
  dir.create(answers_dir, recursive = TRUE, showWarnings = FALSE)

  meta_file <- file.path(assignment_dir, "setup_metadata.xlsx")
  if (!file.exists(meta_file)) {
    stop("No setup metadata found for timestamp: ", run_timestamp)
  }
  meta <- readxl::read_excel(meta_file)

  # Ensure assignment_file column exists
  if (!"assignment_file" %in% names(meta)) {
    meta <- meta |>
      dplyr::mutate(
        assignment_file = file.path(
          assignment_dir,
          paste0(student_id, "_data.xlsx")
        )
      )
  }

  # --- Read all student files (map_dfr) -------------------------------------
  # Each file has: student_id | rxn_substrate | rxn_condition | rxn_time | 0_mM ... 160_mM
  have_files <- meta |> dplyr::filter(file.exists(assignment_file))
  if (nrow(have_files) == 0) {
    stop("No assignment files found to read.")
  }

  student_raw <- purrr::map_df(
    unique(have_files$assignment_file),
    readxl::read_xlsx
  )

  # --- Long format for analysis ---------------------------------------------
  # Parse numeric substrate concentration from column names like "40_mM"
  long_dat <- student_raw |>
    tidyr::pivot_longer(
      cols = tidyselect::matches("^[0-9]+_mM$"),
      names_to = "substrate_label",
      values_to = "absorbance"
    ) |>
    dplyr::mutate(
      substrate_conc = readr::parse_number(substrate_label)
    ) |>
    dplyr::select(
      student_id,
      rxn_substrate,
      rxn_condition,
      rxn_time,
      substrate_conc,
      absorbance
    )

  # --- Estimate rates: slope of absorbance ~ time per (student, condition, [S]) ----
  # Keep 0 mM (rate ~ 0) for MM plotting points; exclude from LB fitting later.
  rate_dat <- long_dat |>
    dplyr::group_by(
      student_id,
      rxn_substrate,
      rxn_condition,
      substrate_conc
    ) |>
    dplyr::summarise(
      rate = {
        fit <- try(stats::lm(absorbance ~ rxn_time), silent = TRUE)
        if (inherits(fit, "try-error")) {
          NA_real_
        } else {
          unname(coef(fit)[["rxn_time"]])
        }
      },
      .groups = "drop"
    )

  # --- Optional MM (nls) fit per student+condition (best-effort, safe if fails) ----
  # Not required, but nice to have for table; if it fails we return NA.
  est_params <- rate_dat |>
    dplyr::group_by(student_id, rxn_substrate, rxn_condition) |>
    dplyr::summarise(
      KM = {
        df <- dplyr::cur_data_all()
        tryCatch(
          {
            # nls over non-negative conc with a small lower bound for stability
            df2 <- df |>
              dplyr::filter(substrate_conc > 0, is.finite(rate), rate >= 0)
            if (nrow(df2) < 3) {
              stop("too few points")
            }
            start <- list(
              Vmax = max(df2$rate, na.rm = TRUE),
              Km = stats::median(df2$substrate_conc, na.rm = TRUE)
            )
            nlsfit <- stats::nls(
              rate ~ Vmax * substrate_conc / (Km + substrate_conc),
              data = df2,
              start = list(Vmax = start$Vmax, Km = start$Km),
              control = list(warnOnly = TRUE)
            )
            round(coef(nlsfit)[["Km"]], 3)
          },
          error = function(e) NA_real_
        )
      },
      VMAX = {
        df <- dplyr::cur_data_all()
        tryCatch(
          {
            df2 <- df |>
              dplyr::filter(substrate_conc > 0, is.finite(rate), rate >= 0)
            if (nrow(df2) < 3) {
              stop("too few points")
            }
            start <- list(
              Vmax = max(df2$rate, na.rm = TRUE),
              Km = stats::median(df2$substrate_conc, na.rm = TRUE)
            )
            nlsfit <- stats::nls(
              rate ~ Vmax * substrate_conc / (Km + substrate_conc),
              data = df2,
              start = list(Vmax = start$Vmax, Km = start$Km),
              control = list(warnOnly = TRUE)
            )
            round(coef(nlsfit)[["Vmax"]], 3)
          },
          error = function(e) NA_real_
        )
      },
      .groups = "drop"
    )

  # --- Plot builders ---------------------------------------------------------
  # Abs vs time: facet by condition; colour by [S]
  plotAbsVsTime <- function(df_long) {
    df_long |>
      dplyr::mutate(
        substrate_f = factor(
          substrate_conc,
          levels = sort(unique(df_long$substrate_conc))
        )
      ) |>
      ggplot2::ggplot(ggplot2::aes(
        x = rxn_time,
        y = absorbance,
        colour = substrate_f
      )) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", se = FALSE) +
      ggplot2::facet_wrap(~rxn_condition, ncol = 1) +
      ggplot2::labs(
        x = "Reaction time (min)",
        y = "Absorbance (AU)",
        colour = "[S] (mM)"
      ) +
      ggplot2::theme_minimal()
  }

  # MM: points for both conditions on one plot (optional smooth from nls estimates if available)
  plotMM <- function(df_rate, df_est) {
    p <- ggplot2::ggplot(
      df_rate,
      ggplot2::aes(x = substrate_conc, y = rate, colour = rxn_condition)
    ) +
      ggplot2::geom_point() +
      ggplot2::labs(
        x = "Substrate concentration (mM)",
        y = "Reaction rate (ΔAbs/min)",
        title = "Michaelis-Menten (points)"
      ) +
      ggplot2::theme_minimal()

    # Optional curve: if both VMAX & KM exist, overlay a curve per condition
    add_curve <- function(p, Vmax, Km, cond, xs) {
      if (is.na(Vmax) || is.na(Km)) {
        return(p)
      }
      ys <- Vmax * xs / (Km + xs)
      p +
        ggplot2::geom_line(
          data = data.frame(
            substrate_conc = xs,
            rate = ys,
            rxn_condition = cond
          ),
          ggplot2::aes(x = substrate_conc, y = rate, colour = rxn_condition)
        )
    }

    # Try to add per-condition curves
    for (cond in unique(df_rate$rxn_condition)) {
      pars <- df_est |> dplyr::filter(rxn_condition == cond)
      if (nrow(pars) == 1) {
        xs <- seq(
          0,
          max(df_rate$substrate_conc, na.rm = TRUE),
          length.out = 200
        )
        p <- add_curve(p, pars$VMAX, pars$KM, cond, xs)
      }
    }
    p
  }

  # LB: lines per condition (lm on reciprocals), dynamic x-axis to include all x-intercepts
  plotLB <- function(df_rate) {
    df_pos <- df_rate |>
      dplyr::filter(substrate_conc > 0, is.finite(rate), rate > 0)
    if (nrow(df_pos) == 0) {
      return(
        ggplot2::ggplot() +
          ggplot2::labs(title = "Lineweaver-Burke: insufficient data") +
          ggplot2::theme_minimal()
      )
    }
    df_lb <- df_pos |>
      dplyr::mutate(x = 1 / substrate_conc, y = 1 / rate)

    # Fit lm per condition and compute x-intercepts (-a/b)
    fits <- df_lb |>
      dplyr::group_by(rxn_condition) |>
      dplyr::summarise(
        a = {
          fit <- try(stats::lm(y ~ x), silent = TRUE)
          if (inherits(fit, "try-error")) NA_real_ else unname(coef(fit))
        },
        b = {
          fit <- try(stats::lm(y ~ x), silent = TRUE)
          if (inherits(fit, "try-error")) NA_real_ else unname(coef(fit))
        },
        .groups = "drop"
      ) |>
      dplyr::mutate(
        x_intercept = ifelse(is.na(a) | is.na(b) | b == 0, NA_real_, -a / b)
      )

    xmin <- min(c(0, fits$x_intercept), na.rm = TRUE)
    xmax <- max(df_lb$x, na.rm = TRUE)
    pad <- 0.1 * (xmax - xmin)
    xlim <- c(xmin - pad, xmax + pad)

    ggplot2::ggplot(df_lb, ggplot2::aes(x = x, y = y, colour = rxn_condition)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
      ggplot2::coord_cartesian(xlim = xlim) +
      ggplot2::labs(
        x = "1 / [S] (1/mM)",
        y = "1 / rate (min/ΔAbs)",
        title = "Lineweaver-Burke (fitted lines)"
      ) +
      ggplot2::theme_minimal()
  }

  # --- Build per-student objects (nest/map) ---------------------------------
  # Join rates + estimates back for convenience
  answers_tbl <- long_dat |>
    dplyr::group_by(student_id, rxn_substrate, rxn_condition) |>
    tidyr::nest(data = c()) |>
    # tidyr::nest(data = dplyr::everything()) |>
    dplyr::left_join(
      rate_dat |>
        dplyr::select(student_id, rxn_substrate, rxn_condition, rate) |>
        dplyr::group_by(student_id, rxn_substrate, rxn_condition) |>
        tidyr::nest(rates = rate),
      # tidyr::nest(rates = dplyr::everything()),
      by = c("student_id", "rxn_substrate", "rxn_condition")
    ) |>
    dplyr::left_join(
      est_params |>
        dplyr::group_by(student_id, rxn_substrate, rxn_condition) |>
        tidyr::nest(est = c(KM, VMAX)),
      by = c("student_id", "rxn_substrate", "rxn_condition")
    ) |>
    dplyr::mutate(
      abs_vs_time_plot = purrr::map(data, plotAbsVsTime),
      mm_plot = purrr::map2(
        rates,
        est,
        ~ {
          df_rate <- .x
          df_est <- if (is.null(.y) || nrow(.y) == 0) {
            tibble::tibble(
              rxn_condition = character(),
              KM = NA_real_,
              VMAX = NA_real_
            )
          } else {
            .y
          }
          plotMM(df_rate, df_est)
        }
      ),
      lb_plot = purrr::map(rates, plotLB),
      # Small tables
      table1 = purrr::map2(
        student_id,
        rxn_substrate,
        ~ tibble::tibble(info = c("student", "substrate"), value = c(.x, .y))
      ),
      table2 = purrr::map(
        est,
        ~ {
          if (is.null(.x) || nrow(.x) == 0) {
            tibble::tibble(rxn_condition = NA, KM = NA, VMAX = NA)
          } else {
            .x |> dplyr::select(rxn_condition, KM, VMAX)
          }
        }
      )
    )

  # Layout helper
  createAnswerPage <- function(table1, table2, p_abs, p_mm, p_lb) {
    g1 <- gridExtra::tableGrob(table1, rows = NULL)
    g2 <- gridExtra::tableGrob(table2, rows = NULL)
    layout <- rbind(c(1, 2), c(3, 3), c(3, 3), c(4, 5), c(4, 5))
    gridExtra::grid.arrange(g1, g2, p_abs, p_mm, p_lb, layout_matrix = layout)
  }

  # --- Write PDFs ------------------------------------------------------------
  if (output_files %in% c("separate", "both")) {
    answers_tbl |>
      dplyr::mutate(
        pdf_path = file.path(answers_dir, paste0(student_id, "_answers.pdf"))
      ) |>
      purrr::pwalk(function(
        student_id,
        rxn_substrate,
        data,
        rates,
        est,
        abs_vs_time_plot,
        mm_plot,
        lb_plot,
        table1,
        table2
      ) {
        grDevices::pdf(
          file.path(answers_dir, paste0(student_id, "_answers.pdf")),
          width = 8.27,
          height = 11.69
        ) # A4 inches
        print(createAnswerPage(
          table1,
          table2,
          abs_vs_time_plot,
          mm_plot,
          lb_plot
        ))
        grDevices::dev.off()
        if (verbose) message("Wrote answers PDF for ", student_id)
      })
  }

  if (output_files %in% c("single", "both")) {
    pdf_all <- file.path(answers_dir, "answers_all_students.pdf")
    grDevices::pdf(pdf_all, width = 8.27, height = 11.69)
    purrr::pwalk(answers_tbl, function(...) {
      row <- tibble::as_tibble(list(...))
      print(createAnswerPage(
        row$table1[[1]],
        row$table2[[1]],
        row$abs_vs_time_plot[[1]],
        row$mm_plot[[1]],
        row$lb_plot[[1]]
      ))
    })
    grDevices::dev.off()
    if (verbose) message("Wrote combined PDF: ", pdf_all)
  }

  invisible(TRUE)
}
