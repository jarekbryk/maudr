utils::globalVariables(c(
  # column names resolved via tidy evaluation (dplyr, tidyr)
  "student_no",
  "first_name",
  "surname",
  "student_id",
  "rxn_substrate",
  "rxn_condition",
  "rxn_time",
  "substrate_conc",
  "substrate_label",
  "absorbance",
  "gradient",
  "rate",
  "Vmax",
  "Km",
  "inhibition_actual",
  "assignment_file",
  "timestamp",
  # assignReactions: .id column created by tidyr::uncount()
  "id",
  # generateAssignments: intermediate mutate columns
  "substrate_conc_char",
  "student_idx",
  # generateAnswers: nested list-columns selected via dplyr::select()
  "rates",
  "est",
  "abs_vs_time_plot",
  "mm_plot",
  "lb_plot",
  "table1",
  "table2",
  # generateAnswers plotAbsVsTime: aes column created in mutate
  "substrate_f",
  # generateAnswers plotLB: aes columns and lm-derived variables
  "x",
  "y",
  "intercept",
  "slope",
  # generateAnswers plotLB: ggpmisc after_stat computed variables
  "eq.label",
  "rr.label"
))
