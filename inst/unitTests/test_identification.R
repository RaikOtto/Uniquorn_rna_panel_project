test_identification = function(){
  
  set.seed(1)

  HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")
  
    ident_result = identify_vcf_file(
        HT29_vcf_file,
        ref_gen = "GRCH37",
        verbose = TRUE
    )
  
  checkTrue( class(ident_result) == "data.frame" )
  
  checkTrue( dim(ident_result)[1] >= 60 )
  checkTrue( dim(ident_result)[2] == 10  )

    checkTrue(
      length( 
          ident_result$Conf_score_sig[
              ident_result$Conf_score_sig == TRUE 
              ]
      ) > 0
    )
    checkTrue( 
        length( 
            ident_result$Conf_score_sig[
                as.logical( ident_result$Conf_score_sig ) == FALSE 
                ]
        ) <= 2100
    )
}