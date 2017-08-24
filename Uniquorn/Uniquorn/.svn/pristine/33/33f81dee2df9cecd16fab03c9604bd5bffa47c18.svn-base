test_identification = function(){
  
  set.seed(1)

  HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")
  
  ident_result = identify_vcf_file( HT29_vcf_file, ref_gen = "GRCH37" )
  
  checkTrue( class(ident_result) == "data.frame" )
  
  checkTrue( dim(ident_result)[1] == 60 )
  checkTrue( dim(ident_result)[2] == 8  )

}