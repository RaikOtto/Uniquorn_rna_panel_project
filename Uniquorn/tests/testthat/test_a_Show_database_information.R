context("Correct access and structure of stored database information")
names = read_library_names(ref_gen = "GRCH37")

test_that("show_contained_ccls returns correct object", {
    names = read_library_names(ref_gen = "GRCH37")
    ccls_all = show_contained_ccls()
  
    expect_is(ccls_all, "data.frame")
    expect_equal(dim(ccls_all)[2], 6)
    expect_equal(sort(unique(ccls_all$Library)), sort(names))
    expect_is(ccls_all[,1], "character")
    expect_is(ccls_all[,2], "character")

})

test_that("show_contained_variants_in_library returns correct object", {
    names = read_library_names(ref_gen = "GRCH37")
    if("CELLMINER" %in% names){
        g_ranges = show_contained_variants_in_library(
            library_name = "CELLMINER"
        )
  
        expect_is(g_ranges, "GRanges")
        expect_equal(length(g_ranges), 57262)
        expect_equal(sum(grepl("CELLMINER", g_ranges$Member_CCLs)), 57262)
    }
})

test_that("show_contained_variants_for_ccl returns the right object", {
    names = read_library_names(ref_gen = "GRCH37")
    if("CELLMINER" %in% names){
        variants = show_contained_variants_for_ccl(
            name_ccl = "786_0", library_name = "CELLMINER"
        )
  
        expect_is(variants, "GRanges")
        expect_equal(length(variants), 584)
  }
})

test_that("show_which_ccls_contain_variant returns the right object", {
    names = read_library_names(ref_gen = "GRCH37")
    if("CELLMINER" %in% names){
        ccls = show_which_ccls_contain_variant(start = 92030762,
            end = 92030762, chromosome = 8, library_name = "CELLMINER")
        names = unlist(strsplit(ccls$Member_CCLs, ",", fixed = TRUE))
  
        expect_is(ccls, "GRanges")
        expect_equal(length(names), 8)
        expect_true(any(grepl("786_0_CELLMINER", names, fixed = TRUE)))
        expect_true(any(grepl("PC_3_CELLMINER", names, fixed = TRUE)))
    }
})