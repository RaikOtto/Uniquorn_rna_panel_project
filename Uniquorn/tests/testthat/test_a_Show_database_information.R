context("Correct access and structure of stored database information")

test_that("reference library contains only CELLMINER", {
  names = read_library_names(ref_gen = "GRCH37")
  
  expect_equal(names, "CELLMINER")
  expect_false(any(grepl("COSMIC", names, fixed = TRUE)))
  expect_false(any(grepl("CCLE", names, fixed = TRUE)))
  expect_false(any(grepl("EGA", names, fixed = TRUE)))
  expect_false(any(grepl("GDC", names, fixed = TRUE)))
})

test_that("show_contained_ccls returns correct object", {
  ccls_all = show_contained_ccls()
  
  expect_message(show_contained_ccls(), "CELLMINER amount CCLs: 60")
  expect_is(ccls_all, "data.frame")
  expect_equal(dim(ccls_all)[1], 60)
  expect_equal(dim(ccls_all)[2], 6)
  expect_equal(unique(ccls_all$Library), "CELLMINER")
  expect_true(ccls_all[1,1] == "786_0")
  expect_true(ccls_all[60,1] == "UO_31")
})

test_that("show_contained_variants_in_library returns correct object", {
  g_ranges = show_contained_variants_in_library(library_name = "CELLMINER")
  
  expect_is(g_ranges, "GRanges")
  expect_equal(length(g_ranges), 57262)
  expect_equal(sum(grepl("CELLMINER", g_ranges$Member_CCLs)), 57262)
  expect_error(show_contained_variants_in_library(library_name = "COSMIC"), "Could not find library!")
})

test_that("show_contained_variants_for_ccl returns the right object", {
  variants = show_contained_variants_for_ccl(name_ccl = "786_0", library_name = "CELLMINER")
  
  expect_is(variants, "GRanges")
  expect_equal(length(variants), 584)
  expect_error(show_contained_variants_for_ccl(
      name_ccl = "786_O",
      library_name = "COSMIC"),
      "Could not find library!")
})

test_that("show_which_ccls_contain_variant returns the right object", {
  ccls = show_which_ccls_contain_variant(start = 92030762, end = 92030762, chromosome = 8, library_name = "CELLMINER")
  names = unlist(strsplit(ccls$Member_CCLs, ",", fixed = TRUE))
  
  expect_is(ccls, "GRanges")
  expect_equal(length(names), 8)
  expect_true(any(grepl("786_0_CELLMINER", names, fixed = TRUE)))
  expect_true(any(grepl("PC_3_CELLMINER", names, fixed = TRUE)))
  expect_error(show_which_ccls_contain_variant(
    start = 92030762,
    end = 92030762,
    chromosome = 8,
    library_name = "COSMIC"),
    "Could not find library!")
})