# maudr

**maudr** is an R package designed to simulate, distribute, and evaluate enzyme kinetics data for teaching and assessment in biochemistry and related laboratory courses. It allows instructors to automatically generate individualized, realistic experimental datasets for each student (complete with enzyme inhibition conditions) and to produce corresponding answers with visual summaries. 

The name of the package is a play on the name of [Dr Maud Menten](https://en.wikipedia.org/wiki/Maud_Menten), a Canadian biochemist who established fundamental principles of enzyme kinetics, along with Dr Leonor Michaelis.

## ⬇️ Installation

**maudr** is only available on GitHub, there are no plans to upload it to CRAN.

``` r
install.packages("remotes")
remotes::install_github("jarekbryk/maudr")
```

## 📈 Usage

### Introduction

The package uses enzyme kinetic parameters for _S. cerevisiae_ alcohol dehydrogenase, as established in a [1987 publication by Ganzhorn _et al._](https://www.jbc.org/article/S0021-9258(18)61419-X/pdf). An Excel file with the parameters for ADH and its inhibitors is provided by default (`reaction_parameters.xlsx`), however, parameters for other enzymes and inhibitors can be provided by the user, as long as the column names remain unchanged:

- `rxn_substrate`
- `Kcat`
- `Km`
- `Vmax`
- `enzyme_conc`
- `inhibition_actual`

If the default parameters are used, the only input file required for the package is an Excel file with a list of students for whom the datasets will be generated and analysed. This file must have one row per student, with the following column headers:

- `student_no` – unique student ID (e.g. "u123456")
- `first_name` – student's first name
- `surname` – student's surname

### Process

The package has four main functions to produce students' datasets and analysis of their results; the functions should be used in order.

#### `initialiseProject()`

It sets up a top-level folder (default location - current folder) where the input data will be stored and output produced. The top-level folder will be named "maudr_assignments" unless a different name is specified by the user. Folders `data` and `output` are created, with subfolders `output/assignments_output` and `output/answers_outputa. R projects are not used.

#### `assignReactions()`

It creates a metadada table with an assignment of every student-to-reaction conditions, randomly drawn from the reaction_parameters.xlsx file. This table is deposited in the `output/assignments_output/` folder in the top-level project folder created by `initialiseProject()`.

The output of `assignReactions()` is a list with a timestamp and the metadata table with student-reaction assignments. The timestamp is used as an input for the next function.

#### `generateAssignments()`

#### `generateAnswers()`

### ▶️ Example use

```r
library(maudr)
initialiseProject(path = "~/Desktop")
setup <- assignReactions(student_file = "~/Documents/students_file.xlsx", project_path = "~/Desktop") # students_file.xlsx is copied into the data folder in the top level folder
# str(setup) # to check the list created by the assignReaction() function
generateAssignments(run_timestamp = setup$timestamp, project_path = "~/Desktop") # Assignments (one Excel file per student) will be deposited in the output/assignment_output folder
generateAnswers(run_timestamp = setup$timestamp, project_path = "~/Desktop", output_files = "both") # Assignments (one PDF file per student plus a single PDf with answers from all students' datasets) will be deposited in the output/answers_output folder
```
