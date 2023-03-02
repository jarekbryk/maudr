# Calculate reaction rate V
calculateV <- function(Vmax, Km, substrate_conc){
	V = Vmax  * substrate_conc / (Km + substrate_conc)
	# V = Vmax  *  ifelse(substrate_conc > 0, jitter(substrate_conc), substrate_conc)  /  (Km    +    ifelse(substrate_conc > 0, jitter(substrate_conc), substrate_conc))
	return(V)
}

# Calculate gradients for each reaction rate
calculateGradient <- function(Vmax, Km, substrate_conc){
	V = Vmax  *  substrate_conc / (Km + substrate_conc)
	# V = Vmax  *  ifelse(substrate_conc > 0, jitter(substrate_conc), substrate_conc)  /  (Km    +    ifelse(substrate_conc > 0, jitter(substrate_conc), substrate_conc))
	gradient = V * enz_vol / 1000000 / cuv_vol * eps
	return(gradient)
}

# Fit non-linear model to Michaelis-Menten curve and estimate Km and Vmax
# WARNING: The nls function won't fit the model work unless you round the rxn_rate values (floating point error?)
estimateKmVmax <- function(df){
  nls(rxn_rate ~ estimated_Vmax * substrate_conc / (estimated_Km + substrate_conc), data = df, start = list(estimated_Km = 5, estimated_Vmax = 0.05)) %>%
    tidy() %>%
    select(parameter = term, estimate) %>%
    mutate(estimate = round(estimate, digits = 2)) %>% # Important - see note above
    pivot_wider(names_from = parameter, values_from = estimate)
}

# Estimate the Km and Vmax from rxn_rate and substrate_conc
# The nls function doesn't work unless you round the rxn_rate values!!! (floating point error?)
estimated_Km_Vmax <- students_rxn_params %>%
  mutate(rxn_rate = pmap(select(., Vmax, Km, substrate_conc), calculateV)) %>%  # Calculation of V for different substrate concentrations
  unnest(c(rxn_rate, substrate_conc)) %>%
  mutate(rxn_rate = round(rxn_rate, digits = 2)) %>%
  group_by(student_id, rxn_substrate, inhibition_actual) %>%
  nest() %>%
  mutate(estimated_params = map(data, estimateKmVmax)) %>%
  select(-data) %>%
  mutate(estimated_params = bind_cols(inhibition_type = inhibition_actual, estimated_params)) %>%
  unnest(estimated_params) %>%
  ungroup() %>%
  select(-inhibition_actual) %>%
  group_by(student_id, rxn_substrate) %>%
  nest() %>%
  summarise(estimated_params = list(data)) %>% # Is this necessary?
  unnest(estimated_params)

# Generate data files for each student
generateAbsVsTimeData <- function(student_id, df) {
  df %>%
    mutate(student_id = student_id,
           substrate_conc_char = paste(substrate_conc, "mM", sep = "_")) %>%
    select(-c(inhibition_actual, substrate_conc)) %>%
    pivot_wider(names_from = substrate_conc_char, values_from = absorbance) %>%
    select(student_id, rxn_substrate, rxn_condition, rxn_time, everything()) %>%
    write_xlsx(here("output", paste(student_id, "data_new.xlsx", sep = "_")))
}

# Create the student data
students_rxn_params %>%
  mutate(gradient = pmap(select(., Vmax, Km, substrate_conc), calculateGradient),
         rxn_time = list(time)) %>%  # Calculate of gradients for different reaction rates
  unnest(c(rxn_time)) %>%
  unnest(c(substrate_conc, gradient)) %>% # Expand the dataframe to have multiple time points for each gradients
  mutate(absorbance = round(gradient * rxn_time, digits = 3),
         rxn_condition = case_when(
           inhibition_actual == "no_inhibition" ~ "without_inhibitor",
           TRUE ~ "with_inhibitor")) %>%
  relocate(rxn_condition, .after = inhibition_actual) %>%
  select(student_id, rxn_substrate, inhibition_actual, rxn_condition, substrate_conc, rxn_time, absorbance) %>%
  group_by(student_id) %>%
  nest() %$%
  map2(student_id, data, generateAbsVsTimeData)


