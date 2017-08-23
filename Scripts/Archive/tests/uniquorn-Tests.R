library("Uniquorn")
library("testthat")

HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")

ident_result = identify_vcf_file( HT29_vcf_file, ref_gen = "GRCH37", verbose = TRUE )

expect_that( ident_result, is_a("data.frame") )

expect_that( dim(ident_result)[1], is_more_than( 59 ) )
expect_that( dim(ident_result)[1], is_less_than( 2100 ) )
expect_that( dim(ident_result)[2], is_more_than( 9  ) )
expect_that( dim(ident_result)[2], is_less_than( 15  ) )

expect_that( 
    length( 
        ident_result$Conf_score_sig[
            ident_result$Conf_score_sig == TRUE 
        ]
    ),
    is_more_than( 0  ) 
)
expect_that( 
    length( 
        ident_result$Conf_score_sig[
            ident_result$Conf_score_sig == TRUE 
            ]
    ),
    is_less_than( 3  ) 
)