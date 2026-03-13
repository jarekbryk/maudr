---
title: "maudr: R package to generate and evaluate enzyme kinetics data for teaching in undergraduate laboratory courses."
tags:
  - R
  - enzyme kinetics
  - education
  - biochemistry
  - simulation
authors:
  - name: Lewis T. J. Ward
    orcid: 
    affiliation: 1
  - name: Dr Richard Bingham
    orcid: 0000-0003-3805-0364
    affiliation: 1
  - name: Wayne-hon Lau
    orcid: 
    affiliation: 1
  - name: Dr Shamus Burns
    orcid: 
    affiliation: 1
  - name: Dr Jarosław Bryk
    orcid: 0000-0002-2791-9709
    affiliation: 1
affiliations:
  - name: University of Huddersfield
    index: 1
date: "2026-03-13"
bibliography: paper.bib
---

# Summary

maudr is an R package designed to generate, distribute, and evaluate enzyme kinetics data for teaching and assessment in undergraduate biochemistry laboratory courses. Given a list of students as the only required input, it produces individualized experimental datasets based on real enzyme parameters, complete with enzyme inhibition conditions, and automatically generates corresponding answer PDFs with visual summaries. The package uses yeast alcohol dehydrogenase with five alcohol substrates as a default, but can be extended to any enzyme for which basic processing parameters are known. The name of the package is a play on the name of Dr Maud Menten, a Canadian biochemist who established fundamental principles of enzyme kinetics, along with Dr Leonor Michaelis.

# Statement of Need

Educational resources for teaching enzyme kinetics in undergraduate courses have been recently summarised and discussed by Mak and colleagues (2024)). These resources tend to focus on the calculations of enzymatic transformations (REF) or simulations of the reactions across a range of kinetic models (REF), and are less concerned with the practical logistics of providing appropriate datasets to large cohorts of students and ensuring accurate marking of work where students perform calculations on and visualise enzyme kinetics.

maudr fills this gap by generating an arbitrary number of individualised student datasets in Excel format, matching data typically acquired in an "enzyme inhibition" laboratory class, and automatically produces visual and numerical answers in PDF format to facilitate marking and feedback. The package requires only a basic ability to run R functions and by default provides widely used set of kinetic parameters for yeast alcohol dehydrogenase and its inhibitors (REF Ganzhorn).

Originally, the set of scripts that became the maudr package was developed during COVID-19 lockdowns, when students could not attend laboratory classes but still needed to practice enzyme kinetics. Since then, the scripts have been routinely used to generate datasets and solutions as formative exercises for first- and second-year undergradautes, during revision sessions and in other contexts, such as in programming classes for biology students. One of us (JB) presented the idea for the package in 2024 on the useR! conference in Salzburg, which provided enough encouragement for maudr's development.

# Functionality

In a typical use scenario, an instructors provides a list of students and maudr randomly assigns each student a substrate and an inhibitor from the default reaction parameters table, and then simulates absorbance vs time data for both an inhibited and non-inhibited reaction across substrate concentrations from 0 and 160 mM, using the enzyme's Kcat, Km and Vmax (REF Ganzhorn). A small amount of user-controllable random noise in is added to each dataset by default to simulate experimental variation (this option can be switched off). Each student receives a single Excel file with individualised data, with information of the tytpe of inhibition hidden.

Students are expected to use their file to analyse and visualise reaction rates as absorbance vs time across substrate concentrations, to estimate Km and Vmax from a Michaelis-Menten plot and identify the inhibition type from a Lineweaver-Burk double-reciprocal plot. 

As individual datasets would then require solutions to be calculated and visualised by instructors to check students' solutions, maudr generates the corresponding answers, containing all three plots and a table of estimated kinetic parameters and expected inhibition type, automatically from each students' file. This output can be saved either as a one-page PDF per student or as a single collated PDF with all students' expected solutions, or both.

# Acknowledgements

Use of the `nls` function REF to plot and estimate Michaelis-Menten curves was inspired by [this 2015 post by prof. Paul Brennan describing prof. Rob Benton's script on how to use it](https://rforbiochemists.blogspot.com/2015/05/plotting-and-fitting-enzymology-data.html). We are also grateful to [Fonti Kar](https://github.com/fontikar) for helpful suggestions that improved the package.

# References
