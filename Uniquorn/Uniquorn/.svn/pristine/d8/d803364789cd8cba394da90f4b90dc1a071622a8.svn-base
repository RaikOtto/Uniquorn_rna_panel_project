test_identification = function(){
  
  set.seed(1)

  HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")
  
  ident_result = identify_vcf_file( HT29_vcf_file, ref_gen = "GRCH37", verbose = TRUE )
  
  checkTrue( class(ident_result) == "data.frame" )
  
  checkTrue( dim(ident_result)[1] >= 60 )
  checkTrue( dim(ident_result)[2] == 10  )

  checkTrue( ( as.logical( ident_result$Conf_score_sig[1] )  == rep( TRUE, 1 ) ) || FALSE  )
  checkTrue( ( as.logical( ident_result$Conf_score_sig[2:60] ) ==( rep( FALSE, dim(ident_result)[1]-1) ) ) || FALSE )
}