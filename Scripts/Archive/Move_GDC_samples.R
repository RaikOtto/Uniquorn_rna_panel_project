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

## cellminer

cellminer_raw_samples = list.files("~/Uniquorn_data/benchmark_vcf_files/raw_files/",pattern = ".CELLMINER.vcf", full.names = T)

for( i_file  in cellminer_raw_samples){
    
    print(i_file)
    vcf_t = read.table(i_file, comment.char = "", fill = T, header = F,stringsAsFactors = F, col.names = 1:10)
    
    header = vcf_t[ grep(vcf_t[,1], pattern = "^##"),]
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    
    vcf_geno = vcf_t[ grep(vcf_t[,1], pattern = "^##",invert = T),]
    
    header_line = vcf_geno[1,]
    vcf_geno = vcf_geno[-1,]

    chroms = as.character(vcf_geno[,1])
    chroms[ grep(chroms, pattern = "^chr", invert = T)] = paste("chr",chroms[ grep(chroms, pattern = "^chr", invert = T)], sep ="")
    
    vcf_geno[,1] = chroms
    header_line[1] = paste( paste( c(header[,1]), collapse= "\r\n", sep = ""), header_line[1], sep ="\r\n" ) 
    
    vcf_geno = rbind( header_line, vcf_geno )
    
    write.table(vcf_geno, i_file, sep ="\t",quote = F, row.names = F, col.names = F)
}


#i_files = list.files("~/Uniquorn_data/benchmark_vcf_files/raw_files/",pattern = ".vcf",full.names = T)
i_files = list.files("~/Uniquorn_data/benchmark_vcf_files/panel_ident_files_regularized/0.5/",pattern = ".ident.tsv",full.names = T,recursive = F)

for(i_file in i_files){
  
    if(str_detect(i_file,pattern = "GDC")) 
        library_name = "GDC"
    if(str_detect(i_file,pattern = "EGA")) 
        library_name = "EGA"
    if(str_detect(i_file,pattern = "COSMIC")) 
        library_name = "COSMIC"
    if(str_detect(i_file,pattern = "CCLE")) 
        library_name = "CCLE"
    if(str_detect(i_file,pattern = "CELLMINER")) 
        library_name = "CELLMINER"
    
    #i_file_mod = str_replace( i_file, pattern = paste(c(".",library_name,".vcf"),collapse= "",sep =""),"")
    i_file_mod = str_replace( i_file, pattern = paste(c(".",library_name,".ident.tsv"),collapse= "",sep =""),"")
    i_file_mod = str_replace_all(i_file_mod, pattern = "\\.", "_")
    #i_file_mod = paste(i_file_mod, paste(c(".",library_name,".vcf"),collapse= "",sep =""),sep = "")
    i_file_mod = paste(i_file_mod, paste(c(".",library_name,".ident.tsv"),collapse= "",sep =""),sep = "")
    i_file_mod = str_replace(i_file_mod,pattern = "/1_0/","/1.0/")
    i_file_mod = str_replace(i_file_mod,pattern = "/0_5/","/0.5/")
    i_file_mod = str_replace(i_file_mod,pattern = "/0_25/","/0.25/")
    i_file_mod = str_replace(i_file_mod,pattern = "/0_0/","/0.0/")
    file.rename(i_file,i_file_mod)
}

