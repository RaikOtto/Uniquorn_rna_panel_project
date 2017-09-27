#' Build Fingerprint From VCF File
#' 
#' Intern utility function. Wrapper function to filter VCF file using 
#' \code{parse_vcf_file} and build fingerprint with identifier.
#' 
#' @param vcf_file_path character string giving the path to the vcf file
#'  on the operating system.
#' @param sim_list_stats table containing all cancer cell lines currently
#' contained in the database.
#' @param cl_id character string giving the cancer cell line identifier.
#' @param n_threads integer specifying the number of threads to be used.
#' @import BiocParallel
#' @usage 
#' add_single_file(vcf_file_path, sim_list_stats, cl_id, n_threads)
#' @return Loci-based DNA-mutational fingerprint of the cancer cell line
#'  associated with the respective identifier.

add_single_file = function( 
    vcf_file_path,
    sim_list_stats,
    cl_id
){
    if(any(sim_list_stats[, CL == cl_id])){ 
        stop("Fingerprint with name ", cl_id, 
             " already present in database. Please change the", 
             " name or remove the old cancer cell line.")
    }
    vcf_fingerprint = parse_vcf_file(vcf_file_path)
    vcf_fingerprint = data.table( 
        "Fingerprint" = vcf_fingerprint,
        "CL"  = rep(cl_id, length(as.character(vcf_fingerprint)) )
    )
    return(vcf_fingerprint)
}