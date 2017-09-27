#' add_ccl_variants_to_monet_db
#
#' This function adds the variants of parsed custom CCLs to a monet DB instance
#'  
#' @param vcf_input_files a character vector containing the input vcf files.
#'  This may be one or many vcf files.
#' @param ref_gen a character string specifying the reference genome version.
#'  All training sets are associated with a reference genome version.
#'  Default is \code{"GRCH37"}. 
#' @param library a character string giving the name of the library to add the
#'  cancer cell lines to. Default is \code{"CUSTOM"}. 
#'  Library name will be automatically added as a suffix to the identifier.
#' @param n_threads an integer specifying the number of threads to be used.
#' @param test_mode Is this a test? Just for internal use
#' @return Message wheather the adding was successful
#' @import BiocParallel doParallel bit bitOpts
#' @usage 
#' add_custom_vcf_to_database(vcf_input_files, ref_gen = "GRCH37", library = "CUSTOM",
#'                            n_threads = 1, test_mode = FALSE)
#' @examples 
#' HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package = "Uniquorn");
#' add_custom_vcf_to_database(vcf_input_files = HT29_vcf_file,
#'                            library = "CUSTOM",
#'                            ref_gen = "GRCH37",
#'                            n_threads = 1,
#'                            test_mode = TRUE)
add_ccl_variants_to_db_bitwise = function(
    vcf_input_files,
    ref_gen = "GRCH37",
    library = "CUSTOM",
    n_threads = 1,
    test_mode = FALSE)
{
    
    vcf_fingerprint = parse_vcf_file(vcf_input_file)
    
    for ( chrom in c(1:22,"X","Y")){
        if (n_threads > 1){
            doParallel::registerDoParallel(n_threads)
            foreach::foreach(
                vcf_input_file = vcf_input_files
            ) %dopar% {
                write_vcf_file_bitwise(vcf_fingerprint, ref_gen = ref_gen, test_mode = test_mode, library = library, chrom = chrom)
            }
            doParallel::stopImplicitCluster()
        } else {
          
            for (vcf_input_file in vcf_input_files){
                write_vcf_file_bitwise(vcf_fingerprint, ref_gen = ref_gen, test_mode = test_mode, library = library, chrom = chr)
            }
        }
    }
}

#' write_vcf_file_to_monet_db
#
#' This function adds the variants of parsed custom CCLs to a monet DB instance
#'  
#' @param vcf_input_files a character vector containing the input vcf files.
#'  This may be one or many vcf files.
#' @param ref_gen a character string specifying the reference genome version.
#'  All training sets are associated with a reference genome version.
#'  Default is \code{"GRCH37"}. 
#' @param library a character string giving the name of the library to add the
#'  cancer cell lines to. Default is \code{"CUSTOM"}. 
#'  Library name will be automatically added as a suffix to the identifier.
#' @param n_threads an integer specifying the number of threads to be used.
#' @param test_mode Is this a test? Just for internal use
#' @return Message wheather the adding was successful
#' @import DBI dplyr
#' @usage 
#' write_vcf_file_to_monet_db(vcf_input_files, ref_gen = "GRCH37", library = "CUSTOM",
#'  cl_id, n_threads = 1, test_mode = FALSE)
#' @examples 
#' HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package = "Uniquorn");
#' write_vcf_file_to_monet_db(vcf_input_files = HT29_vcf_file,
#'                            library = "CUSTOM",
#'                            ref_gen = "GRCH37",
#'                            cl_id = "Test_id"
#'                            n_threads = 1,
#'                            test_mode = TRUE)
write_vcf_file_bitwise = function(
  vcf_fingerprint,
  ref_gen = "GRCH37",
  library,
  chrom,
  test_mode = FALSE)
{
  
    cl_id = gsub("^.*/", "", vcf_input_file)
    cl_id = gsub(".vcf", "", cl_id, fixed = TRUE)
    cl_id = gsub(".hg19", "", cl_id, fixed = TRUE)
    cl_id = paste0(cl_id, "_", library)
    cl_id = toupper(cl_id)
    cl_id = stringr::str_replace_all(cl_id, pattern = "\\.", "_")
  
    length_frame = data.frame(
        "chr1"=247197891,"chr2"=242713278,"chr3"=199439629,"chr4"=191246650,"chr5"=180727832,
        "chr6"=170735623,"chr7"=158630410,"chr8"=146252219,"chr9"=140191642,"chr10"=135347681,
        "chr11"=134361903,"chr12"=132289533,"chr13"=114110907,"chr14"=106354309,"chr15"=100334282,
        "chr16"=88771793,"chr17"=78646005,"chr18"=76106388,"chr19"=63802660,"chr20"=62429769,
        "chr21"=46935585, "chr22"=49396972,"chrX"=154908521,"chrY"=57767721)
    
    library(dplyr)
    package_path = system.file("", package = "Uniquorn")
    database_path =  paste( c( package_path,"/","Uniquorn.RSQLite"), sep ="", collapse= "")
    con = DBI::dbConnect(RSQLite::SQLite(), path = database_path)
    
    message(paste(c("Sample: ",cl_id,", Chr: ",chr),collapse = "", sep =""))
    table_name = paste( library, chrom,sep = "_" )
    
    vars_chr = vcf_fingerprint[ grep(vcf_fingerprint, pattern = paste(c( "^",chrom,"_"),collapse = "",sep ="")) ]
    starts = as.integer(sapply(vars_chr, FUN = function(vec){return(as.character(unlist(str_split(vec,pattern = "_")))[2])}))
    ends = as.integer(sapply(vars_chr, FUN = function(vec){return(as.character(unlist(str_split(vec,pattern = "_")))[3])}))
    starts_tmp = starts
    starts = starts[ (!is.na(starts_tmp)) | (!is.na(ends)) ]
    ends = ends[(!is.na(starts_tmp)) | (!is.na(ends)) ]

    chrom_length = length_frame[[paste("chr",chrom, sep ="")]]
    mut_vec = as.raw( rep(0, chrom_length) )
    
    mut_vec[ c(starts,ends) ] = as.raw(1)

    try( 
      expr = paste( c("mut_data = tbl(con, ",table_name,")"), sep = "", collapse = ""),
      silent=T
    )
    
    if (! exists("mut_data")){ 
      mut_mat = data.frame( matrix(as.raw(), ncol= ))
    } else {
      
    }
    mut_mat$cl_id = mut_vec

    copy_to(con, mut_data, name = table_name)
    
    file.remove(database_path_wait)
}
