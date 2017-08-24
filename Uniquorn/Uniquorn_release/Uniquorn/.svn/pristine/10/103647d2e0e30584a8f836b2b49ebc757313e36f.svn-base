#' parse_cosmic_genotype_data
#' 
#' Parses cosmic genotype data
#' 
#' @param sim_list Variable containing mutations & cell line
#' @param cosmic_file Path to cosmic clp file in hard disk
#' @return The R Table sim_list which contains the CoSMIC CLP fingerprints 
parse_cosmic_genotype_data = function( cosmic_file, sim_list ){

    split_coords = function(vec){
        
        split_vec = stringr::str_split( vec, ":" )
        chrom = split_vec[1]
    }
    
    exclude_cols_cosmic = c(rep("NULL",4),"character",rep("NULL",18),"character",rep("NULL",14))
    

    if ( ! grepl("CosmicCLP_MutantExport", cosmic_file)){ # MutantExport

        stop("Warning. This is not the recommended COSMIC genotype file! The recommended file is the 'CosmicCLP_MutantExport.tsv.gz' file.")
        exclude_cols_cosmic = c(rep("NULL",4),"character",rep("NULL",13),"character",rep("NULL",13))
    }
    
    cosmic_genotype_tab = utils::read.csv2( cosmic_file, sep ="\t", colClasses = exclude_cols_cosmic)
    
    coords = as.character( sapply( cosmic_genotype_tab[,2], FUN = stringr::str_replace_all, ":|-", "_" ) )
    cls    = stringr::str_replace_all( stringr::str_to_upper(cosmic_genotype_tab[,1]), "/|(|])| ", "" )
    cls[ cls == "KM-H2" ] = "KMH2"
    cls[ cls == "KMH-2" ] = "KMH2ALTERNATIVE"
    cls    = sapply( cls, FUN = function( cl_name ){ return( paste(cl_name, "COSMIC", sep = "_") ) } )
    
    unify = function( CL, coords, cls ){
        
        Fingerprint = unique( coords[ cls == CL] )
        return_val = cbind(Fingerprint, CL)
        
        return( return_val )
    }
    
    new_coords = sapply( unique(cls), coords, cls, FUN = unify )
    
    tmp_list = lapply( new_coords, "names<-", value = c("V1", "V2"))
    new_sim_list = do.call("rbind", lapply(tmp_list, data.frame, stringsAsFactors = FALSE))
    sim_list = rbind( sim_list, new_sim_list )
 
    return(sim_list)
}

#' parse_ccle_genotype_data
#' 
#' Parses ccle genotype data
#' 
#' @param sim_list Variable containing mutations and cell line
#' @param ccle_file Path to CCLE file on hard disk
#' @return The R Table sim_list which contains the CCLE fingerprints
parse_ccle_genotype_data = function( ccle_file, sim_list ){
    
    exclude_cols_ccle = c(rep("NULL",4),rep("character", 3),rep("NULL",8),"character",rep("NULL",35))
    
    ccle_genotype_tab = utils::read.csv2( ccle_file, sep ="\t", colClasses = exclude_cols_ccle)
    
    coords = as.character( apply(
        ccle_genotype_tab,
        FUN = function( vec ){ return( paste0( c(
            as.character( vec[1] ),
            as.character( vec[2] ),
            as.character( vec[3] )
            ), collapse = "_" ) ) },
        MARGIN = 1
    ) )
    
    cls    = sapply( ccle_genotype_tab[,4], FUN = stringr::str_split, "_" )
    cls    = as.character( sapply( cls, FUN = function(vec){ return(vec[1]) }) )
    cls    = as.character( sapply( cls, FUN = function( cl_name ){ return( paste(cl_name, "CCLE", sep = "_") ) } ) )
    
    new_sim_list = data.frame( coords, cls )
    colnames(new_sim_list) = colnames(sim_list)
    sim_list = rbind( sim_list, new_sim_list )
    
    return(sim_list)
}

#' split_add
#' 
#' @param vcf_matrix_row row of the vcf file
#' @return Transformed entry of vcf file, reduced to start and length
split_add = function( vcf_matrix_row ){
    
    reference  = as.character( unlist( vcf_matrix_row[4] ) )
    length_ref = length( unlist( stringr::str_split( reference  ,"") ))
    variations = as.character( unlist( stringr::str_split( unlist( vcf_matrix_row[5] ), "," ) ) )
    length_variations = length(unlist( stringr::str_split( unlist( vcf_matrix_row[5] ), "," ) ) )
    
    chromosome = rep( as.character( vcf_matrix_row[1] ), length(variations)  )
    start      = as.integer( rep( vcf_matrix_row[2], length(variations)  ) )
    
    fingerprint = as.character()
    for ( i in seq( length_variations ) ){ 
        
        start_var  = as.character( start[i] )
        variation  = variations[i]
        length_var = length( unlist( stringr::str_split( variation,"") ))
        
        length_alt = max( length_var, length_ref )
        end_var    = as.character( as.integer( as.integer(start_var) + length_alt - 1 ) )
        
        chrom      = stringr::str_replace( stringr::str_to_upper( stringr::str_trim( as.character( unlist( chromosome[i] ) ) ) ), "CHR", "" )
        
        fingerprint = c( 
            fingerprint, 
            paste0( c( 
                as.character(chrom),
                as.character(start_var),
                as.character(end_var)
            ),
            collapse = "_" )
        )
    }
    
    return( fingerprint )
}

#' parse_vcf_file
#' 
#' Parses the vcf file and filters all information except for the start and length of variations/ mutations.
#' 
#' @param vcf_file_path Path to the vcf file on the operating system
#' @usage 
#' parse_vcf_file( vcf_file_path )
#' @return Loci-based DNA-mutational fingerprint of the cancer cell line as found in the input VCF file 
parse_vcf_file = function( vcf_file_path  ){
    
    if ( base::file.exists( vcf_file_path ) ){
        
        vcf_handle = utils::read.table( vcf_file_path, sep ="\t", header = FALSE, comment.char = "#", fill = TRUE )
        
        fingerprint = apply( vcf_handle, FUN = split_add, MARGIN = 1  )
        fingerprint = unique( as.character( unlist(fingerprint ) ) )
        
        return( fingerprint )
        
    } else {
        
        stop( paste0( "Did not find VCF file: ", vcf_file_path  ) )
    }
}