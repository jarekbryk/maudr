# Randomly assign substrate and inhibition type to each student
assignReactionConditions <- function(
		students = load("student_names.rda"), # list of students
		rxn_params = load("reaction_parameters.rda"), # table with reaction parameters
		time = c(.17, .33, .5, .66, .83, 1), # measurements timepoints (min)
		substr_conc = c(0,2,4,8,10,20,40,80,160), # substrate concentration (mM)
		cuv_vol = 0.003, # volume of cuvette (l)
		eps = 6220, # extinction coefficient
		enz_vol = 0.1, # volume of enzyme added (ml)
		save = TRUE){ # do you want to write the prepared set to disk?

	seed  <- sample(100:999, 1) # Record the seed to recreate the data later if needed
	set.seed(seed)
	print(paste0("The seed value for this round is: ", seed, ". If save is set to FALSE, record this number for future reference"))

	students_rxn_params <- students %>%
		mutate(student_id = paste(toupper(student_no), toupper(first_name), sep = "_"),
					 rxn_substrate = sample(rxn_substr, nrow(students), replace = TRUE),
					 inhibition_actual = sample(inh_type, nrow(students), replace = TRUE)) %>%
		uncount(., 2, .id = "id") %>%
		mutate(inhibition_actual = ifelse(id == 2, "no_inhibition", inhibition_actual)) %>%
		select(-student_no, -first_name, -surname, -id) %>%
		left_join(., rxn_params, by = c("rxn_substrate", "inhibition_actual")) %>%
		mutate(substrate_conc = list(substr_conc)) # Add substrate concentrations as list-column

	if(save == TRUE){
		write_csv(x = students_rxn_params, file = here("output", paste0(Sys.Date(), "_", students-rxn-params, "_seed_", seed, ".csv"))) # Save the file with the conditions for each student - useful for debugging and checking later on
		return(students_rxn_params)
	}
	else {
		return(students_rxn_params)
	}
}
