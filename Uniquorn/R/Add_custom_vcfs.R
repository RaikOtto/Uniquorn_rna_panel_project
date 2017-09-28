#' add_new_ccls_to_ccl_library
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
#' @import BiocParallel
#' @usage 
#' add_new_ccls_to_ccl_library(vcf_input_files, ref_gen = "GRCH37", library = "CUSTOM",
#'                            n_threads = 1, test_mode = FALSE)
#' @examples 
#' HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package = "Uniquorn");
#' add_new_ccls_to_ccl_library(vcf_input_files = HT29_vcf_file,
#'                            library = "CUSTOM",
#'                            ref_gen = "GRCH37",
#'                            n_threads = 1,
#'                            test_mode = TRUE)
#' @export
add_new_ccls_to_ccl_library = function(
    vcf_input_files,
    ref_gen = "GRCH37",
    library_name = "CUSTOM",
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
                g_query,
                ref_gen = ref_gen,
                test_mode = test_mode,
                library_name = library_name
            )
        }
        
        doParallel::stopImplicitCluster()
        
    } else {
      
        for (vcf_input_file in vcf_input_files){
            
            g_query = parse_vcf_file(vcf_input_file)
            parse_vcf_query_into_db(
                g_query,
                ref_gen = ref_gen,
                test_mode = test_mode,
                library_name = library_name
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
#' @param library_name a character string giving the name of the library to add the
#'  cancer cell lines to. Default is \code{"CUSTOM"}. 
#'  Library name will be automatically added as a suffix to the identifier.
#' @param n_threads an integer specifying the number of threads to be used.
#' @param test_mode Is this a test? Just for internal use
#' @return Message wheather the adding was successful
#' @import GenomicRanges stringr
#' @usage 
#' write_vcf_file_to_monet_db(vcf_input_files, ref_gen = "GRCH37", library_name = "CUSTOM",
#'  cl_id, n_threads = 1, test_mode = FALSE)
#' @examples 
#' HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package = "Uniquorn");
#' write_vcf_file_to_monet_db(vcf_input_files = HT29_vcf_file,
#'                            library_name = "CUSTOM",
#'                            ref_gen = "GRCH37",
#'                            cl_id = "Test_id"
#'                            n_threads = 1,
#'                            test_mode = TRUE)
parse_vcf_query_into_db = function(
    g_query,
    ref_gen = "GRCH37",
    library_name,
    test_mode = FALSE)
{
    
    cl_id = mcols( g_query )$Member_CCLs[1]
    message(paste(c("Sample: ",cl_id,", Library: ",library_name),collapse = "", sep =""))
    
    package_path = system.file("", package = "Uniquorn")
    rdata_path = paste( c( package_path,"/",library_name,"_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    library_path =  paste( c( package_path,"/Libraries_Ref_gen_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    
    # save new vcf
    
    if (!file.exists(rdata_path)){
        g_mat = GenomicRanges::GRanges(
          seqnames = NULL,
          IRanges(
            start = NULL,
            end = NULL
          )
        )
    } else {
        g_mat = readRDS(rdata_path)
    }
    
    g_mat_new = unique(c( GenomicRanges::GRanges(g_mat), GenomicRanges::GRanges(g_query)))
    mcols(g_mat_new)$Member_CCLs = rep("",nrow(mcols(g_mat_new)))
        
    fo_query = findOverlaps(
        query = g_query,
        subject = g_mat_new,
        select = "arbitrary",
        type = "equal"
    )
    
    fo_g_mat = findOverlaps(
        query = g_mat,
        subject = g_mat_new,
        select = "arbitrary",
        type = "equal"
    )
    
    mcols( g_mat_new )$Member_CCLs[fo_query] = elementMetadata(g_query)$Member_CCLs
    mcols( g_mat_new )$Member_CCLs[fo_g_mat]   =
        paste( 
            mcols( g_mat_new )$Member_CCLs[fo_g_mat],
            elementMetadata(g_mat)$Member_CCLs, sep = ","
        )
    mcols( g_mat_new )$Member_CCLs = str_replace( mcols( g_mat_new )$Member_CCLs, pattern = "^,","" )
    g_mat = g_mat_new
        
    if (! test_mode)
        saveRDS(g_mat,rdata_path)
    
    # synchronize libraries
    library_path =  paste( c( package_path,"/Libraries_Ref_gen_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    try( expr = "library_names = readRDS(library_path)")
    if (!exists("library_names")){
        library_names = c(library_name)
    } else {
        library_names = unique(c(library_names, library_name))
    }
    saveRDS(library_names,library_path)
}