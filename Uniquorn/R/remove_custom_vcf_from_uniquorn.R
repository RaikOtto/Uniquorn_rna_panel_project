#' Remove Cancer Cell Line
#'
#' This function removes a cancer cell line training fingerprint (VCF file) from the database. The names of all training sets can 
#' be seen by using the function \code{show_contained_cls}.
#' 
#' @param name_cl a character vector giving the names of the cancer cell line
#'  identifiers to be removed. Can be one or many
#' @param ref_gen a character vector specifying the reference genome version.
#'  All training sets are associated with a reference genome version. 
#'  Default is \code{"GRCH37"}.
#' @param test_mode is this a test? Just for internal use
#' @import DBI
#' @usage 
#' remove_custom_vcf_from_database(name_cl, ref_gen = "GRCH37", test_mode = FALSE)
#' @examples 
#' remove_custom_vcf_from_database(name_cl = "HT29_CELLMINER",
#'                                 ref_gen = "GRCH37",
#'                                 test_mode = TRUE)
#' @return Message that indicates whether the removal was succesful.
#' @export
remove_custom_vcf_from_database = function( 
    name_cl,
    ref_gen = "GRCH37",
    test_mode = FALSE
){
    message("Reference genome: ", ref_gen)
    
    sim_list_stats = initiate_db_and_load_data(request_table = "sim_list_stats",
                                               subset = "*", ref_gen = ref_gen)
    sim_list       = initiate_db_and_load_data(request_table = "sim_list",
                                               subset = c("FINGERPRINT", "CL"),
                                               ref_gen = ref_gen)
    for(cl_id in name_cl){
      
      cl_id = toupper(cl_id)
      if(any(sim_list_stats[, CL == cl_id])){
        message(paste0("Found CL ", cl_id,
                       ", removing from database."))
      } else {
        stop(paste0("No training set for cancer cell line found for the name: ",
                    cl_id))
      }
    sim_list = sim_list[CL != cl_id]
    }
    message("Removed all samples. Re-calculating the Cancer cell line data")
    
    res_vec = re_calculate_cl_weights(sim_list = sim_list, ref_gen = ref_gen)
    
    message("Finished aggregating, saving to database")
    
    write_data_to_db(content_table = as.data.frame(res_vec[1]),
                     "sim_list",
                     ref_gen = ref_gen, 
                     overwrite = TRUE,
                     test_mode = test_mode)
    write_data_to_db(content_table = as.data.frame(res_vec[2]),
                     "sim_list_stats",
                     ref_gen = ref_gen,
                     overwrite = TRUE, 
                     test_mode = test_mode)
    
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