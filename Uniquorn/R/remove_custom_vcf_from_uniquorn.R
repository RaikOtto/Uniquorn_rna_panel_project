#' Remove Cancer Cell Line
#'
#' This function removes a cancer cell line training fingerprint (VCF file) from 
#' the database. The names of all training sets can 
#' be seen by using the function \code{show_contained_cls}.
#' 
#' @param ccl_names A character vector giving the names of the cancer cell line
#'  identifiers to be removed. Can be one or many
#' @param ref_gen A character vector specifying the reference genome version.
#'  All training sets are associated with a reference genome version. 
#'  Default is \code{"GRCH37"}.
#' @param library_name Name of the library from which the ccls are to be removed
#' @param test_mode Signifies if this is a test run
#' @import GenomicRanges stringr IRanges
#' @usage 
#' remove_ccls_from_database(ccl_names, ref_gen = "GRCH37", test_mode = FALSE)
#' @examples 
#' remove_ccls_from_database(ccl_names = c("HT29_CELLMINER"),
#'                                 ref_gen = "GRCH37",
#'                                 test_mode = TRUE)
#' @return Message that indicates whether the removal was succesful.
#' @export
remove_ccls_from_database = function( 
    ccl_names,
    ref_gen = "GRCH37",
    library_name,
    test_mode = FALSE
){
    message("Reference genome: ", ref_gen)
  
    if (length(library_name)> 1)
        stop("Cannot process more than one library per run, please 
            provide single libraries.")
    
    g_mat = read_mutation_grange_objects(
        library_name = library_name,
        ref_gen = ref_gen,
        mutational_weight_inclusion_threshold = 0
    )
    
    member_ccls = as.character(GenomicRanges::mcols(g_mat)$Member_CCLs)
    
    for (ccl_id in ccl_names){
      
        cl_lib = paste(ccl_id, library_name, sep = "_")
        message(paste("Deleting ", cl_lib))
        search_term = paste(c(
            paste( "(",cl_lib,",",")", sep = "" ),
            paste( "(",cl_lib,")", sep = "" )
        ),collapse= "|")
        
        member_ccls = sapply(
            member_ccls, 
            FUN = str_replace_all, 
            pattern = search_term,
            replacement = ""
        )
    }
    
    non_empty_vec = member_ccls != ""
    g_mat = g_mat[non_empty_vec]
    member_ccls = member_ccls[ non_empty_vec ]
    GenomicRanges::mcols(g_mat)$member_ccls = member_ccls
    
    message(paste(
      c( "Excluded ",
         as.character(sum(non_empty_vec == FALSE)),
         " variant entries after removal."),
        collapse = "", sep ="" )
    )
    
    message(paste0("Finished removing all ccls. Recalculating DB"))
    
    write_w0_and_split_w0_into_lower_weights(
        g_mat = g_mat,
        ref_gen = ref_gen,
        library_name = library_name
    )
    
    message(paste0("Finished removing all cancer cell lines"))
}


#' Remove entire Library from Database
#'
#' This function removes a entire library from the database by removing all 
#' associated cancer cell line fingerprints from the database.
#' 
#' @param library a character vector giving the names of the library to be
#'  removed.
#' @param ref_gen a character vector specifying the reference genome version.
#'  All training sets are associated with a reference genome version. 
#'  Default is \code{"GRCH37"}.
#' @param test_mode is this a test? Just for internal use.
#' @import DBI
#' @usage 
#' remove_library_from_database(library, ref_gen = "GRCH37", test_mode = FALSE)
#' @examples 
#' remove_custom_vcf_from_database(library = "CELLMINER",
#'                                 ref_gen = "GRCH37",
#'                                 test_mode = FALSE)
#' @return Message that indicates whether the removal was succesful.
#' @export
remove_library_from_database = function( 
  library,
  ref_gen = "GRCH37",
  test_mode = FALSE
){
  message("Reference genome: ", ref_gen)
  
  sim_list_stats = initiate_db_and_load_data(request_table = "sim_list_stats",
                                             subset = "*", ref_gen = ref_gen)
  sim_list       = initiate_db_and_load_data(request_table = "sim_list",
                                             subset = c("FINGERPRINT", "CL"),
                                             ref_gen = ref_gen)
  
  if(any(sim_list_stats[, CL %like% library])){
    message(paste0("Found ", sum(sim_list_stats[, CL %like% library])), 
                   " training sets for library ", library, ". Removing all.")
  } else {
    stop("No training set for library ", library, " found in database.")
  }
  sim_list = sim_list[!(CL %like% library)]

  message("Removed all samples. Re-calculating the Cancer cell line data.")
  res_vec = re_calculate_cl_weights(sim_list = sim_list, ref_gen = ref_gen)
  message("Finished aggregating, saving to database.")
  
  write_data_to_db(content_table = as.data.frame(res_vec[1]),
                   "sim_list",
                   ref_gen = "GRCH37", 
                   overwrite = TRUE,
                   test_mode = test_mode)
  write_data_to_db(content_table = as.data.frame(res_vec[2]),
                   "sim_list_stats",
                   ref_gen = "GRCH37",
                   overwrite = TRUE, 
                   test_mode = test_mode)
  
  message("Removing library ", library, " and all associated",
          " cancer cell lines done.")
}