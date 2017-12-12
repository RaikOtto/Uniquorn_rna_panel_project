#library("Uniquorn")
#library("testthat")
context("Correct identification results for test data")

test_that("identification of test data returns correct results", {
  HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package = "Uniquorn")
  ident_result = identify_vcf_file(HT29_vcf_file, ref_gen = "GRCH37", verbose = FALSE)
  
  expect_is(ident_result, "data.frame")
  expect_equal(dim(ident_result)[1], 60)
  expect_equal(dim(ident_result)[2], 7)
  expect_equal(as.character(unique(ident_result$Library)), "CELLMINER")
  expect_equal(as.character(ident_result[1,1]), "HT29")
  expect_true(ident_result$Identification_sig[1])
  expect_false(any(ident_result$Identification_sig[2:60]))
  expect_true(ident_result$P_value_sig[1])
  expect_false(any(ident_result$P_value_sig[2:60]))
})