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
#' @import BiocParallel GenomicRanges IRanges
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
add_custom_vcf_to_database = function(
    vcf_input_files,
    ref_gen = "GRCH37",
    library_name = "CUSTOM",
    n_threads = 1,
    test_mode = FALSE)
{
    
    #ccl_list = read_ccl_list(ref_gen = ref_gen, library_name = library_name)
  
    if (n_threads > 1){
      
        doParallel::registerDoParallel(n_threads)
      
        foreach::foreach(
            vcf_input_file = vcf_input_files
        ) %dopar% {
            g_query = parse_vcf_file(
                vcf_input_file,
                ref_gen = ref_gen
            )
            
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
            
            g_query = parse_vcf_file(
                vcf_input_file,
                ref_gen = ref_gen
            )
            
            parse_vcf_query_into_db(
                g_query,
                ref_gen = ref_gen,
                test_mode = test_mode,
                library_name = library_name
            )
        }
    }
    print("Finished")
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
#' @import GenomicRanges stringr IRanges
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
    test_mode = FALSE
){
    
    cl_id = mcols( g_query )$Member_CCLs[1]
    cl_id = str_replace_all( cl_id, pattern = paste("_",library_name,sep = ""),"" )
    
    cl_data =  show_contained_cls( verbose = FALSE)
    if (cl_id %in% cl_data$CCL[ cl_data$Library == library_name  ] ){
        print(
            paste(c("CCL ", cl_id, " already contained in DB ",library_name,
                    ". Remove first or change name"),
            sep="",collapse = "")
        )
        return()
    }
    
    message(paste(c("Sample: ",cl_id,", Library: ",library_name),collapse = "", sep =""))
    
    g_mat = read_mutation_grange_objects(
        library_name = library_name,
        ref_gen = ref_gen,
        mutational_weight_inclusion_threshold = 0
    )
    
    if ( "Member_CCLs" %in% names(mcols(g_mat))  ){
      
        g_mat_new = unique(c( GenomicRanges::GRanges(g_mat), GenomicRanges::GRanges(g_query)))
        fo_g_mat = findOverlaps(
            query = g_mat,
            subject = g_mat_new,
            select = "arbitrary",
            type = "equal"
        )
        
    } else {
        
        g_mat_new = unique(GenomicRanges::GRanges(g_query))
        fo_g_mat = GenomicRanges::GenomicRangesList()
    }
    mcols(g_mat_new)$Member_CCLs = rep("",nrow(mcols(g_mat_new)))
    
    fo_query = findOverlaps(
        query = g_query,
        subject = g_mat_new,
        select = "arbitrary",
        type = "equal"
    )
    
    old_members = as.character( 
        sapply( as.character(elementMetadata(g_query)$Member_CCLs), function(vec){
                ccl_ids = as.character(unlist(str_split(vec,pattern = ",")))
                ccl_ids = unique(ccl_ids)
                ccl_ids = paste(ccl_ids, collapse = ",", sep ="")
            return(ccl_ids)
    }))
      
    mcols( g_mat_new )$Member_CCLs[fo_query] = old_members
    
    if ( "Member_CCLs" %in% names(mcols(g_mat)) )
    
        mcols( g_mat_new )$Member_CCLs[fo_g_mat]   =
          paste( 
            mcols( g_mat_new )$Member_CCLs[fo_g_mat],
            elementMetadata(g_mat)$Member_CCLs, sep = ","
          )
    
    
    mcols( g_mat_new )$Member_CCLs = stringr::str_replace( 
      mcols( g_mat_new )$Member_CCLs,
      pattern = "^,",""
    )
    g_mat = g_mat_new
    
    write_w0_and_split_w0_into_lower_weights(
        g_mat = g_mat,
        ref_gen = ref_gen,
        library_name = library_name
    )
    
    print(paste(c("Finished parsing ",cl_id,", library: ",library_name),sep="",collapse= ""))
  }
    