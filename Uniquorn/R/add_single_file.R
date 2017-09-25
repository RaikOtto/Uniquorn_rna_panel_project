# pre processing

add_single_file = function( 
    vcf_file_path,
    library,
    sim_list,
    sim_list_stats,
    n_threads
){

    # reading vcf
    
    vcf_name = str_to_upper( basename( vcf_file_path) )
    vcf_name = str_replace_all(vcf_name, pattern = "\\.","_")
    vcf_name = str_replace_all(vcf_name, pattern = "_VCF","")
    vcf_full_name = paste(vcf_name, library, sep ="_")
    
    name_present = base::grepl( vcf_full_name, sim_list_stats$CL )
    if( sum( name_present ) > 0 ){ 
       print( paste0( c("Fingerprint with name ",vcf_full_name, 
                      " already present in database. Please change the 
                        name or remove the old cancer cell line."), 
                    collapse = "" )  )
        return(as.data.frame(matrix(as.character(), ncol=2,nrow=0)))
    } else {
    
        print(paste(c("Adding CCL ",vcf_name," to library ",library),sep="",collapse= ""))
        
        vcf_fingerprint = as.character( parse_vcf_file(vcf_file_path, n_threads ) )
        vcf_fingerprint_length = length(vcf_fingerprint)
        
        vcf_fingerprint = data.frame( 
            "Fingerprint" = vcf_fingerprint,
            "CL"  = rep(
                vcf_full_name,
                vcf_fingerprint_length
            )
        )
        
        print(paste0("Finished parsing file: ", vcf_file_path))

        return(vcf_fingerprint)
    }
}