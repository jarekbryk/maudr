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
date: "2025-10-30"
bibliography: paper.bib
---

# Summary

maudr is an R package designed to generate, distribute, and evaluate enzyme kinetics data for teaching and assessment in undergraduate biochemistry laboratory courses. It allows instructors to produce individualized, experimental datasets based on real enzyme parameters, complete with enzyme inhibition conditions, and to produce corresponding answers with visual summaries. The package uses yeast alcohol dehydrogenase with five alcohol substrates as a default, but can be extended to any enzyme for which basic processing parameters are known. The name of the package is a play on the name of Dr Maud Menten, a Canadian biochemist who established fundamental principles of enzyme kinetics, along with Dr Leonor Michaelis.

# Statement of Need

There are several educational resources related to teaching enzyme kinetics in undergraduate courses (recently summarised and discussed by Mak and colleagues (2024)), but they often concentrate on the calculations of enzymatic transformations (REF) or simulations of the reactions for a range of kinetic models (REF), or both. These resources are less focused on the "logistics" of providing appropriate datasets to large cohorts of students and ensuring accurate marking of work where students perform calculations on and visualise enzyme kinetics.

maudr fills this gap by providing a mechanism to generate an arbitrary number of individualised datasets for students in Excel that match data typically acquired by students in a biochemistry "enzyme inhibition" laboratory class. Students can then be tasked with visualising reaction rates and deriving Km and Vmax values from a Michaelis-Menten plot, as well as identifying an unknown inhibition type using a Lineweaver-Burk plot. maudr automatically generates these visual and numerical answers from each of students' datasets in a PDF format, facilitating marking and feedback on students work. The package by default provides widely used set of kinetic parameters of yeast alcohol dehydrogenase and its inhibitors (REF) and requires only a basic ability to run functions in the statistical programming language R (REF). 

Originally, the set of scripts that became the maudr package was developed during COVID-19 lockdowns, when students could not attend laboratory classes but still needed to practice topics in enzyme kinetics. Since then, we have been using maudr to routinely generate datasets and solutions as formative exercises for undergradaute students in year 1 and 2, during revision sessions or in other contexts, such as in programming classes for biology students. One of us (JB) presented the idea for the package in 2024 on the useR! conference in Salzburg, which provided enough encouragement for maudr's development.

# Functionality

A typical use scenario of the package is for instructors to produce individual datasets that contain absorbance values for an inhibited and non-inhibited reaction of yeast alcohol dehydrogenase over a range of substrate concentrations and with different inhibitors. The enzyme parameters for each substrate are taken from REF and students' data is randomised from the reference table included as a default in the package. Information about type of inhibition is also hidden in students' datasets.

Given a list of students - the only required user-supplied input file - the package will randomly assign substrate and inhibitor to each student and then generate absorbance values for a non-inhibited and inhibited reactions over substrate concentrations between 0 and 160 mM, using the enzyme's Kcat, Km and Vmax. The individual datasets by default also include a small amount of user-controllable variation in absorbance values, to simulate experimental variation (this option can be switched off). maudr will generate a single Excel file with the absorbance data for each student.

The files can then be provided to students for analysis and to produce visualisations of enzyme reaction rates as absorbance vs time, dependent on substrate concentration, to plot Michaelis-Menten curves to estimate Km values and to identify the type of inhibition using the Lineweaver-Burk double-reciprocal plot.

As individual datasets analysed by students would then require solutions to be calculated and visualised by instructors to check students' solutions, the package also produces and plots estimated Km and Vmax values and comparisons of inhibited and non-inhibited conditions from individual students' datasets. Solutions for each student are generated as PDF files (one page per student), with an option of saving a separate file for each student or to collate all individual PDFs into a single file, or both.

# Acknowledgements

Use of the `nls` function REF to plot and estimate Michaelis-Menten curves was inspired by [this 2015 post by prof. Paul Brennan describing prof. Rob Benton's script on how to use it](https://rforbiochemists.blogspot.com/2015/05/plotting-and-fitting-enzymology-data.html). We are also grateful to [Fonti Kar](https://github.com/fontikar) for helpful suggestions that improved the package.

# References
