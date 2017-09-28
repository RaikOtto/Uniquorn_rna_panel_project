#' add_ccl_variants_to_db
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
add_ccl_variants_to_db = function(
    vcf_input_files,
    ref_gen = "GRCH37",
    library = "CUSTOM",
    n_threads = 1,
    test_mode = FALSE)
{
    
    if (n_threads > 1){
      
        doParallel::registerDoParallel(n_threads)
      
        foreach::foreach(
            vcf_input_file = vcf_input_files
        ) %dopar% {
            g_query = parse_vcf_file(vcf_input_file)
            parse_vcf_query_into_db(
                vcf_fingerprint,
                ref_gen = ref_gen,
                test_mode = test_mode,
                library = library
            )
        }
        
        doParallel::stopImplicitCluster()
        
    } else {
      
        for (vcf_input_file in vcf_input_files){
            
            g_query = parse_vcf_file(vcf_input_file)
            parse_vcf_query_into_db(
                vcf_fingerprint,
                ref_gen = ref_gen,
                test_mode = test_mode,
                library = library
            )
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
parse_vcf_query_into_db = function(
    g_query,
    ref_gen = "GRCH37",
    library,
    test_mode = FALSE)
{
    
    cl_id = gsub("^.*/", "", vcf_input_file)
    cl_id = gsub(".vcf", "", cl_id, fixed = TRUE)
    cl_id = gsub(".hg19", "", cl_id, fixed = TRUE)
    cl_id = paste0(cl_id, "_", library)
    cl_id = toupper(cl_id)
    cl_id = stringr::str_replace_all(cl_id, pattern = "\\.", "_")
    mcols( g_query )$Member_CCLs = rep(cl_id,length(g_query@seqnames ))
    message(paste(c("Sample: ",cl_id,", Library: ",library),collapse = "", sep =""))
    
    package_path = system.file("", package = "Uniquorn")
    rdata_path = paste( c( package_path,"/",library,"_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    library_path =  paste( c( package_path,"/Libraries_Ref_gen_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    
    # save new vcf
    
    try( expr = "g_mat = readRDS(rdata_path)")
    
    if (! exists("g_mat")){
      
        g_mat = g_query
        colnames(mcols(g_mat)) = "Member_CCLs"
      
    } else {
      
        g_mat_new = union(g_mat,g_query)
        mcols(g_mat_new)$Member_CCLs = rep("",nrow(mcols(g_mat_new)))
        
        fo_query = findOverlaps(
            query = g_query,
            subject = g_mat_new,
            select = "arbitrary"
        )
        fo_g_mat = findOverlaps(
            query = g_mat,
            subject = g_mat_new,
            select = "arbitrary"
        )
        mcols( g_mat_new )$Member_CCLs[fo_query] = elementMetadata(g_query)$Member_CCLs
        mcols( g_mat_new )$Member_CCLs[fo_g_mat]   =
            paste( mcols( g_mat_new )$Member_CCLs[fo_g_mat], elementMetadata(g_mat)$Member_CCLs, sep = ",")
        mcols( g_mat_new )$Member_CCLs = str_replace( mcols( g_mat_new )$Member_CCLs, pattern = "^,","" )
        g_mat = g_mat_new
        
    }
    
    if (! test_mode)
        saveRDS(g_mat,rdata_path)
    
    # synchronize libraries
    
    try( expr = "libraries = readRDS(library_path)")
    if (!exists("libraries")){
        libraries = c(library)
    } else {
        libraries = unique(c(libraries, library))
    }
    saveRDS(libraries,library_path)
}