library("Uniquorn")
library("testthat")

HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")

ident_result = identify_vcf_file( HT29_vcf_file, ref_gen = "GRCH37", verbose = TRUE )

expect_that( ident_result, is_a("data.frame") )

expect_that( dim(ident_result)[1], equals( 60  ) )
expect_that( dim(ident_result)[2], equals( 10  ) )

expect_that( as.logical( ident_result$Conf_score_sig[1] ),  equals( rep(TRUE,1)  ) )
expect_that( as.logical( ident_result$Conf_score_sig[2:60] ), equals( rep(FALSE,60-1)  ) )
