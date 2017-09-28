match_query_ccl_to_database_bitwise = function(
    vcf_fingerprint,
    ref_gen = "GRCH37",
    chr
){
  
    length_frame = data.frame(
        "chr1"=247197891,"chr2"=242713278,"chr3"=199439629,"chr4"=191246650,"chr5"=180727832,
        "chr6"=170735623,"chr7"=158630410,"chr8"=146252219,"chr9"=140191642,"chr10"=135347681,
        "chr11"=134361903,"chr12"=132289533,"chr13"=114110907,"chr14"=106354309,"chr15"=100334282,
        "chr16"=88771793,"chr17"=78646005,"chr18"=76106388,"chr19"=63802660,"chr20"=62429769,
        "chr21"=46935585, "chr22"=49396972,"chrX"=154908521,"chrY"=57767721)
    
    package_path    = system.file("", package = "Uniquorn")
  
    message(paste(c("Chr: ",chr),collapse = "", sep =""))
    
    database_path =  paste( c( package_path,"/", "Uniquorn_variant_bit_vectors.",ref_gen,".",paste("chr",chr, sep =""),".RData"),collapse = "", sep ="/" )
    
    vars_chr = vcf_fingerprint[ grep(vcf_fingerprint, pattern = paste( c( "^",chr,"_"),collapse = "",sep ="")) ]
    vars_chr = vars_chr[!is.na(vars_chr)]
    starts = as.integer(sapply(vars_chr, FUN = function(vec){
        return(as.character(unlist(stringr::str_split(as.character(vec),pattern = "_")))[2])
    }))
    ends = as.integer(sapply(vars_chr, FUN = function(vec){
        return(as.character(unlist(stringr::str_split(as.character(vec),pattern = "_")))[3])
    }))
    starts_tmp = starts
    starts = starts[ (!is.na(starts_tmp)) | (!is.na(ends)) ]
    ends = ends[(!is.na(starts_tmp)) | (!is.na(ends)) ]
    
    mut_vec = as.bit(rep(FALSE, length_frame[[chr]]))
    mut_vec[ c(starts,ends) ] = TRUE
    mut_vec = as.bit(mut_vec)
    
    # reading database
    
    mut_vecs = readRDS(database_path)
    cl_id = names(mut_vecs)[2]
    sum( ( as.bit(mut_vec) & as.bit(mut_vecs[[cl_id]]) ))
    sapply( names(mut_vecs), FUN = function(mut_vec_ccl){
        #return(sum( ( as.bit(mut_vec) & as.bit(mut_vecs[[mut_vec_ccl]]) )))
        print(sum(as.bit(mut_vecs[[mut_vec_ccl]])))
    } )
    
}