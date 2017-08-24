# pre processing

add_single_file = function( 
    vcf_file_path,
    library,
    sim_list,
    sim_list_stats,
    n_threads
){

    # reading vcf
    
    name_present = base::grepl( library, sim_list_stats$CL )
    
    if( sum( name_present ) > 0 ){ 
        stop( paste0( c("Fingerprint with name ",library, 
                        " already present in database. Please change the 
                        name or remove the old cancer cell line."), 
                      collapse = "" )  )
    }

    vcf_fingerprint = as.character( parse_vcf_file(vcf_file_path, n_threads ) )
    
    vcf_fingerprint = data.frame( 
        "Fingerprint" = vcf_fingerprint,
        "CL"  = rep( library, length(as.character(vcf_fingerprint)) )
    )
    
    print(paste0("Finished parsing file: ", vcf_file_path))

    return(vcf_fingerprint)
}