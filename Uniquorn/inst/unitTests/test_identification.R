test_identification = function(){
  
  set.seed(1)

  HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")
  
    ident_result = identify_vcf_file(
        HT29_vcf_file,
        ref_gen = "GRCH37",
        verbose = TRUE,
        robust_mode = TRUE
    )
  
  RUnit::checkTrue( class(ident_result) == "data.frame" )
  
  RUnit::checkTrue( dim(ident_result)[1] >= 60 )
  RUnit::checkTrue( dim(ident_result)[2] == 8  )

  RUnit::checkTrue(
      length( 
          ident_result$Identification_sig[
              ident_result$Identification_sig == TRUE 
              ]
      ) == 3
    )
  RUnit::checkTrue( 
        length( 
            ident_result$Identification_sig[
                as.logical( ident_result$Identification_sig ) == FALSE 
                ]
        ) <= 2100
    )
  
  RUnit::checkTrue(number_ccls[names(number_ccls) == "CELLMINER"] == 60)
  
  add_custom_vcf_to_database( 
      vcf_input_files = HT29_vcf_file,
      library_name = "CELLMINER"
  )
  
  contained_ccls = show_contained_cls()
  number_ccls = table( contained_ccls$Library )
  
  RUnit::checkTrue(number_ccls[names(number_ccls) == "CELLMINER"] == 61)
  
  remove_ccls_from_database(
      ccl_names = "HT29_GZ",
      library_name = "CELLMINER"
  )
  
  contained_ccls = show_contained_cls()
  number_ccls = table( contained_ccls$Library )
  
  RUnit::checkTrue(number_ccls[names(number_ccls) == "CELLMINER"] == 60)
}