library(Uniquorn)

for (i_file in list.files("~/Uniquorn_data/vcfHG19/", pattern = ".vcf$", full.names = T)){
  
  ident_file = paste0(i_file,"_uniquorn_ident.tab")
    
  if (! file.exists(ident_file)){
      print(i_file)
      identify_vcf_file(i_file)
  }
}
