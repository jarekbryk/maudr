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



# Functionality

A typical use scenario of the package is for instructors to produce individual datasets that contain absorbance values for an inhibited and non-inhibited reaction of yeast alcohol dehydrogenase over a range of substrate concentrations and with different inhibitors. The enzyme parameters for each substrate are taken from REF. 

Given a list of students - the only required user-supplied input file - the package will randomly assign substrate and inhibitor to each student and then generate absorbance values for a non-inhibited and inhibited reactions over substrate concentrations between 0 and 160 mM, using the enzyme's Kcat, Km and Vmax. The individual datasets can also include a small amount of user-controllable variation in absorbance values, to simulate experimental variation. maudr will generate a single Excel file with the absorbance data for each student.

The files can then be provided to students for analysis and to produce visualisations of enzyme reaction rates as absorbance vs time, dependent on substrate concentration, to plot Michaelis-Menten curves to estimate Km values and to identify the type of inhibition using the Lineweaver-Burk double-reciprocal plot.

As individual datasets analysed by students would then require solutions to be calculated and visualised by instructors to check students' solutions, the package also produces the plots estimated Km and Vmax values and comparisons of inhibited and non-inhibited conditions from individual students' datasets. Solutions for each student are generated as PDF files (one page per student), with an option of saving a separate file for each student or to collate all individual PDFs into a single file, or both.

# Experience of Use

We have been using maudr (or earlier, the set of scripts that became maudr) to generate datasets and solutions as formative exercises for undergradaute students in year 1 and 2, during revision sessions or in other contexts, such as in programming classes for biology students.

# Project Origin

Originally, the set of scripts that became the maudr package was developed during COVID-19 lockdowns, when students could not attend laboratory classes but still needed to practice topics in enzyme kinetics. One of us (JB) presented the idea for the package in 2024 on the useR! conference in Salzburg, which provided enough encouragement for maudr's development.

# Acknowledgements

Use of the `nls` function REF to plot and estimate Michaelis-Menten curves was inspired by [this 2015 post by prof. Paul Brennan describing prof. Rob Benton's script on how to use it](https://rforbiochemists.blogspot.com/2015/05/plotting-and-fitting-enzymology-data.html).

# References
