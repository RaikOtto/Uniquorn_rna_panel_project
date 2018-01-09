context("Correct add and remove of custom Cancer Cell Line")

test_that("adding and removing a custom Cancer Cell Line works succesfully", {
  # Check add
  HT29_vcf_file = system.file("extdata/HT29_TEST.vcf", package = "Uniquorn")
  add_custom_vcf_to_database(vcf_input_files = HT29_vcf_file, library_name = "CELLMINER")
  ccls_all = show_contained_ccls()
  variants = show_contained_variants_for_ccl(name_ccl = "HT29_TEST", library_name = "CELLMINER")
  
  expect_message(show_contained_ccls(), "CELLMINER amount CCLs: 61")
  expect_equal(dim(ccls_all)[1], 61)
  expect_true(any(grepl("HT29_TEST", ccls_all[,1], fixed = TRUE)))
  expect_equal(length(variants), 24969)
  
  # Check remove
  remove_ccls_from_database(ccl_names = "HT29_TEST", library_name = "CELLMINER")
  ccls_all = show_contained_ccls()
  variants = show_contained_variants_for_ccl(name_ccl = "HT29_TEST", library_name = "CELLMINER")
  
  expect_message(show_contained_ccls(), "CELLMINER amount CCLs: 60")
  expect_equal(dim(ccls_all)[1], 60)
  expect_false(any(grepl("HT29_TEST", ccls_all[,1], fixed = TRUE)))
  expect_equal(length(variants), 0)
})