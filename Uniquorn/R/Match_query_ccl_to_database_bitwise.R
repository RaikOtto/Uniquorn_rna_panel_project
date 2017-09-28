match_query_ccl_to_database_bitwise = function(
    vcf_fingerprint,
    ref_gen = "GRCH37",
    library,
    chrom
){
  
    package_path    = system.file("", package = "Uniquorn")
    message(paste(c("Sample: ",cl_id,", Chr: ",chrom),collapse = "", sep =""))
    
    try( expr = "mut_mat = readRDS(rdata_path)")
    
    if (! exists("mut_mat")){
      
      mut_mat = as.logical(mut_mat)
      #mut_vec = mut_vec[empty_vec_chrom]
      mut_mat = mut_vec
      
    } else {
      
      #empty_vec_chrom = empty_vec_chrom | mut_vec
      mut_mat = cbind(mut_vec)
    }
    mut_mat = as.bit(mut_mat)
    
    vars_chr = vcf_fingerprint[ grep(vcf_fingerprint, pattern = paste(c( "^",chrom,"_"),collapse = "",sep ="")) ]
    starts = as.integer(sapply(vars_chr, FUN = function(vec){return(as.character(unlist(str_split(vec,pattern = "_")))[2])}))
    ends = as.integer(sapply(vars_chr, FUN = function(vec){return(as.character(unlist(str_split(vec,pattern = "_")))[3])}))
    starts_tmp = starts
    starts = starts[ (!is.na(starts_tmp)) | (!is.na(ends)) ]
    ends = ends[(!is.na(starts_tmp)) | (!is.na(ends)) ]
    
    if (max(starts) > chrom_length)
      message(paste(c("Warning, max var position greater than chromosome length: 2**28")))
    
    mut_vec = rep(FALSE, chrom_length)
    mut_vec[ c(starts,ends) ] = TRUE
    
}