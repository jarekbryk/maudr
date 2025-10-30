---
title: "maudr: R package to simulate, distribute, and evaluate enzyme kinetics data for teaching and assessment in undergraduate laboratory courses."
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

When it comes to enzyme kinetics, many students can struggle to grasp the mathematical models and concepts behind it. [@Silverstein2016; @AlOdat2024; @NelsonCox2021] This was further compounded upon by the COVID-19 pandemic depriving students of laboratory practicals. [@Stenson2022] **maudr** is an R package designed to generate, simulate and visualise enzyme kinetics data within an educational context. Educators are able to create realistic, unique datasets with a range of inhibition types, and verify student answers against model answers.

# Statement of Need

Although many enzyme kinetics tools exist, a majority are designed for research rather than teaching. Some can be complex, require programming skills, or are locked behind a paywall. [@Aledo2022; @MakDunn2024; @Johnson2009; @LiorZ2022]

**maudr** provides a lightweight, easy-to-use interface through one single wrapper function. It focuses on pedagogical clarity for educators and learners, and enables the generation of datasets and their visualisations; supporting conceptual learning in undergraduate biochemistry and molecular biology courses.

# Functionality

**maudr** includes the following features:

-   Modular pipeline of functions [@PengFunctions; @WickhamPackages; @WickhamBestPractices; @TaylorForecastML].

-   User-modifiable parameters, as well as built-in defaults.

-   Simulation of Michaelis-Menten kinetic models to generate student data.

-   Visualisation of absorbance-time, Lineweaver-Burk and Michaelis-Menten plots.

-   Exporting of results and figures to .pdf files, formatted using user-determined arguments, for utilisation in teaching materials.

Usage examples and tutorials are included in the package's documentation.

# Experience of Use

**maudr** has been tested on systems running Windows 11 and macOS, and is intended to be integrated into undergraduate biochemistry courses to generate enzyme kinetics datasets for large-scale classes.

# Project Origin

The project originated from the need to provide an interactive, open-source alternative to proprietary simulation software in teaching labs.

# Acknowledgements

Thank you to Wayne-hon Lau, Dr. Shamus Burns, and Dr. Richard Bingham, for their work in the early development of maudr.

# References
