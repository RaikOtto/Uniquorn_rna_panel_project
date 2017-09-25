library(stringr)

gdc_raw_samples = list.files("~/Uniquorn_data/GDC_results//vcf_hg19//", full.names = F)
gdc_raw_samples_path = list.files("~/Uniquorn_data/GDC_results//vcf_hg19//", full.names = T)

ega_raw_samples = list.files("~/Uniquorn_data/EGA_results///vcf_hg19//", full.names = F)
ega_raw_samples_path = list.files("~/Uniquorn_data/EGA_results///vcf_hg19//", full.names = T)

for (i in 1:length(ega_raw_samples)){
    
    i_file = ega_raw_samples[i]
    #if ( ! grepl(i_file,pattern = ".vcf$" ,ignore.case = F) )
    #    continue
    
    #i_file = str_replace(i_file,pattern = "\\.","/")
    #i_file = basename(i_file)
    i_file = str_replace(i_file,pattern = ".hg19","")
    i_file = str_to_upper(i_file)
    i_file_path = dirname(ega_raw_samples_path[i])
    new_name = paste(i_file_path, i_file, sep ="/")
    new_name = str_replace( new_name, pattern = ".VCF$", ".vcf" )
    file.rename(ega_raw_samples_path[i],new_name)
}
