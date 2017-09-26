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
add_custom_vcf_to_database = function(
  vcf_input_files,
  ref_gen = "GRCH37",
  library = "CUSTOM",
  n_threads = 1,
  test_mode = FALSE)
  {
    library = toupper(library)
    if (library == "CUSTOM"){
        message("No library name entered, proceeding with 'CUSTOM' as library")
    }
    
    #Get CLs for ref. genome and all unique panels currently contained
    base::message("Reference genome: ", ref_gen)
    
    sim_list_stats = initiate_db_and_load_data(request_table = "sim_list_stats",
                                               subset = "*", ref_gen = ref_gen)
    sim_list       = initiate_db_and_load_data(request_table = "sim_list",
                                               subset = c("FINGERPRINT", "CL"),
                                               ref_gen = ref_gen)
    panels = sim_list_stats[, unique(gsub("^.*_", "" , unique(CL)))]
    
    if(! library %in% panels){
      message("Specified library not present in database, creating new entry.")
    }
    
    if (nrow(sim_list_stats) == 0)
        warning("Identification might be spurious due to low ",
                "amount of training sample. No cancer cell line data stored ",
                 "at this point for reference genome: ", ref_gen)
    
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
    
    # Add vcf files in parallel
    if (n_threads > 1){

      # Register parallel backend and compute fingerprints in parallel
      doParallel::registerDoParallel(n_threads)
      all_fingerprints = foreach::foreach(
          vcf_input_file = vcf_input_files,
          .combine = rbind
      ) %dopar% {
               
        #Create CL identifier from input file name and library
        cl_id = gsub("^.*/", "", vcf_input_file)
        cl_id = gsub(".vcf", "", cl_id, fixed = TRUE)
        cl_id = gsub(".hg19", "", cl_id, fixed = TRUE)
        cl_id = paste0(cl_id, "_", library)
        cl_id = toupper(cl_id)
                          
        #Create new CL mutational fingerprint from input vcf
        vcf_fingerprint = add_single_file(
          vcf_file_path = vcf_input_file,
          sim_list_stats = sim_list_stats,
          cl_id = cl_id
        )
      }
      doParallel::stopImplicitCluster()
      #Add new fingerprints to existing CL list
      sim_list = rbind(sim_list, all_fingerprints)
    
    } else {
    
      # Add vcf files sequentially
      for (vcf_input_file in vcf_input_files){
      
        #Create CL identifier from input file name and library
        cl_id = gsub("^.*/", "", vcf_input_file)
        cl_id = gsub(".vcf", "", cl_id, fixed = TRUE)
        cl_id = gsub(".hg19", "", cl_id, fixed = TRUE)
        cl_id = paste0(cl_id, "_", library)
        cl_id = toupper(cl_id)
      
        #Create new CL mutational fingerprint from input vcf
        vcf_fingerprint = add_single_file(
            vcf_file_path = vcf_input_file,
            sim_list_stats = sim_list_stats,
            cl_id = cl_id,
            n_threads = n_threads
        )
       
        #Add new fingerprint to existing CL list
        sim_list = rbind(sim_list, vcf_fingerprint)
      }
    }
    
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
