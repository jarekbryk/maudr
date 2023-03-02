## Prepare answer files

I generate the answers from the students' data files simulated data. I am **not** re-loading students' answer files to generate the answers.

### Functions to generate plots

Note that for the Lineweaver-Burke plot the x axis limits need to be dynamically calculated to include the x axis intercept for the non-inhibited data, but I don't know how to do it. Currently the limit is hard-coded to -0.05 but it won't be sufficient in all cases.

# Absorbance vs time plot
plotAbsVsTime <- function(df){
	df %>%
		mutate(gradient = pmap(select(., Vmax, Km, substrate_conc), calculateGradient), # Calculate of gradients for different substrate concentrations
					 rxn_time = list(time)) %>%
		unnest(c(rxn_time)) %>%
		unnest(c(substrate_conc, gradient)) %>% # Expand the dataframe for reaction rate and concentration
		mutate(absorbance = round(gradient * rxn_time, digits = 3), # Calculate absorbance
					 substrate_conc = fct_relevel(as.character(substrate_conc), levels = c("0", "10", "20", "40", "80", "160"))) %>%
		ggplot() + aes(x = rxn_time, y = absorbance, colour = substrate_conc) + geom_point() + geom_smooth(method = "lm", se = FALSE) + facet_wrap(~inhibition_actual) +
		labs(x = "Reaction time (min)",
				 y = "Absorbance (AU)",
				 title = "Absorbance vs reaction time plot",
				 colour = "Substrate concentration (mM)") +
		theme_minimal() +
		theme(legend.position = "bottom")
}

# Michaelis-Menten plot
plotMM <- function(df) {
	df %>%
		mutate(rxn_rate = pmap(select(., Vmax, Km, substrate_conc), calculateV)) %>%  # Calculate reaction rates
		unnest(c(substrate_conc, rxn_rate)) %>% # Expand the dataframe for reaction rate and concentration
		ggplot() + aes(x = substrate_conc, y = rxn_rate, colour = inhibition_actual) + geom_point() +
		labs(x = "Substrate concentration (mM)",
				 y = "Reaction rate (∆abs/min)",
				 title = "Michaelis-Menten plot",
				 colour = "Type of inhibition") +
		theme_minimal() +
		theme(legend.position = "bottom")
}

# Lineweaver-Burke plot
plotLB <- function(df) {
	df %>%
		mutate(rxn_rate = pmap(select(., Vmax, Km, substrate_conc), calculateV)) %>%  # Calculate of gradients for different reaction rates
		unnest(c(substrate_conc, rxn_rate)) %>% # Expand the dataframe for reaction rate and concentration
		filter(substrate_conc != 0) %>%
		mutate(substrate_conc_reciprocal = 1/substrate_conc,
					 rxn_rate_reciprocal = 1/rxn_rate) %>%
		ggplot() + aes(x = substrate_conc_reciprocal, y = rxn_rate_reciprocal, colour = inhibition_actual) + geom_point() + expand_limits(x = -0.05) + geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
		labs(x = "Reciprocal substrate concentration (mM)",
				 y = "Reciprocal reaction rate",
				 title = "Lineweaver-Burke plot",
				 colour = "Type of inhibition") +
		theme_minimal() +
		theme(legend.position = "bottom") +
		stat_poly_eq(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*")))
}

### Generate all the plots

...and add the table with estimated Km and Vmax for the final print

students_answers_data <- students_rxn_params %>%
	left_join(., estimated_Km_Vmax, by = c("student_id", "rxn_substrate")) %>%
	mutate(inhibition_actual = fct_relevel(inhibition_actual, "no_inhibition", after = 0L)) %>%
	group_by(student_id) %>%
	nest() %>%
	mutate(abs_vs_time_plot = map(data, plotAbsVsTime),
				 mm_plot = map(data, plotMM),
				 lb_plot = map(data, plotLB)) %>%
	unnest(data) %>%
	select(student_id, rxn_substrate, inhibition_actual, abs_vs_time_plot, mm_plot, lb_plot, estimated_params) %>%
	filter(inhibition_actual != "no_inhibition")

### Combine all the elements and print the answer file for each student

Answers for each students are saved on a single page PDF and then combined into a multi-page document for ease of printing.

createAnswerPDF <- function(student_id, rxn_substrate, inhibition_actual, abs_vs_time_plot, mm_plot, lb_plot, estimated_params){
	table1 <- tableGrob(tibble("info" = c("student", "substrate"), " " = c(student_id, rxn_substrate)), rows = NULL, cols = NULL)
	table2 <- tableGrob(unnest(estimated_params), rows = NULL)
	layout <- rbind(c(1,2), c(3,3), c(3,3), c(4,5), c(4,5))
	answers <- grid.arrange(table1, table2, abs_vs_time_plot, mm_plot, lb_plot, layout_matrix = layout)
	return(answers)
}

ggsave(plot = marrangeGrob(pmap(students_answers_data, createAnswerPDF), nrow=1, ncol=1), here("output", "SIB2004_answers_all_students.pdf"), width = 210, height = 297, units = "mm", dpi = "retina")

# Solution for separate files for each student
createAnswerPDF2 <- function(student_id, rxn_substrate, inhibition_actual, abs_vs_time_plot, mm_plot, lb_plot, estimated_params){
	table1 <- tableGrob(tibble("info" = c("student", "substrate"), " " = c(student_id, rxn_substrate)), rows = NULL, cols = NULL)
	table2 <- tableGrob(unnest(estimated_params), rows = NULL)
	layout <- rbind(c(1,2), c(3,3), c(3,3), c(4,5), c(4,5))
	answers <- grid.arrange(table1, table2, abs_vs_time_plot, mm_plot, lb_plot, layout_matrix = layout)
	ggsave(plot = answers, path = here("output"), filename = paste0(student_id, "_answers_new.pdf"), units = "mm", width = 210, height = 297, dpi = "retina")
}

pmap(students_answers_data, createAnswerPDF2)

