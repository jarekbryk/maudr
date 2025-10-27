library(testthat)
test_that('demo workflow runs', {
  tmp <- tempdir()
  initialiseProject(path = tmp, verbose = FALSE)
  res <- generateSetup(path = tmp, seed = 4321, verbose = FALSE)
  expect_true(file.exists(file.path(tmp, 'output', 'assignments_output', res$timestamp, 'setup_metadata.xlsx')))
  # now create assignments from that setup
  generateAssignments(res$timestamp, path = tmp, verbose = FALSE)
  # check assignment files
  files <- list.files(file.path(tmp, 'output', 'assignments_output', res$timestamp), pattern = '_data\\.xlsx$', full.names = TRUE)
  expect_true(length(files) > 0)
  # generate answers as single pdf
  generateAnswers(res$timestamp, output_files = 'single', path = tmp, verbose = FALSE)
  expect_true(file.exists(file.path(tmp, 'output', 'answers_output', res$timestamp, 'answers_all_students.pdf')))
})
