#' Add Custom VCF File
#
#' This function adds a custom vcf file to the database. The new file will
#' either be added to an existing library or a new entry will be created if
#' such a library is not already contained in the database.
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
#' @import DBI BiocParallel doParallel
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
#' @export
add_custom_vcf_to_database_Bitwise = function(
    vcf_input_files,
    ref_gen = "GRCH37",
    library = "CUSTOM",
    n_threads = 1,
    test_mode = FALSE)
    {
      library = toupper(library)
      if ((library == "CUSTOM") | (library == "")){
          message("Proceeding with 'CUSTOM' as library name")
      }
      
      #Get CLs for ref. genome and all unique panels currently contained
      base::message("Reference genome: ", ref_gen)
      
      # Check for existance of supplied vcf files
      check_file = function(vcf_input_file){
        if (!file.exists(vcf_input_file)){
          message(paste0("Skipping file, because could not find: ",
                         vcf_input_file))
          return(FALSE)
        } else{
          message(paste0("Found following file, parsing: ", vcf_input_file))
          return(TRUE)
        }
      }
      index_vcf = lapply(vcf_input_files, function(x) check_file(x))
      vcf_input_files = vcf_input_files[unlist(index_vcf)]
      
      add_ccl_variants_to_db_bitwise(
          vcf_input_files = vcf_input_files,
          ref_gen = ref_gen,
          library = library,
          n_threads = n_threads,
          test_mode = test_mode
      )
      
      #Recalculate weights for parsed CLs
      message("Aggregating over parsed Cancer Cell Line data")
      res_vec = re_calculate_cl_weights( 
          sim_list = sim_list, 
          ref_gen = ref_gen
      )
      message("Finished aggregating, saving to database")
      
      if (!test_mode){
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
      }
      message("Finished adding cancer cell lines")
}
