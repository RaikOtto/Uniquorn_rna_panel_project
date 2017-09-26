#' Adds a custom vcf files to the three existing cancer cell line panels
#' @param vcf_input_files Input vcf file.s This may be one or many vcf files
#' @param ref_gen Reference genome version. All training sets are associated 
#' with a reference genome version. Default: GRCH37
#' @param library The name of the library to add the CCLs to. Standard is '_CUSTOM' 
#' will automatically be added as suffix.
#' @param test_mode Is this a test? Just for internal use
#' @param n_threads Specifies number of threads to be used
#' @return Message if the adding has succeeded
#' @import DBI BiocParallel
#' @usage 
#' add_custom_vcfs_to_database( 
#'  vcf_input_files, 
#' ref_gen = "GRCH37", 
#' library = "", 
#' test_mode = FALSE,
#' n_threads = 1)
#' @examples 
#' HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn");
#' add_custom_vcfs_to_database( 
#' vcf_input_files = HT29_vcf_file,
#' library = "",
#' ref_gen = "GRCH37",
#' test_mode = TRUE,
#' n _threads = 1)
#' @export
add_custom_vcfs_to_database = function( 
    vcf_input_files,
    ref_gen = "GRCH37",
    library = "",
    test_mode = FALSE,
    n_threads = 1
    ){
    
    if (n_threads >1)
        BiocParallel::register(BiocParallel::MulticoreParam( n_threads ))
    
    library = stringr::str_to_upper(library)
    
    if ( library == "" ){
        print("No library name entered, please enter a library name")
        library = "CUSTOM"
    }
    
    base::print( base::paste0( "Reference genome: ",ref_gen ) )
    
    sim_list_stats = initiate_db_and_load_data( 
        ref_gen = ref_gen, 
        request_table = "sim_list_stats"
    )
    
    sim_list = initiate_db_and_load_data( 
        ref_gen = ref_gen, 
        request_table = "sim_list"
    )
    sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,]
    sim_list = sim_list[, base::which( 
        colnames(sim_list) != "Ref_Gen"
    ) ]
    sim_list = sim_list[, base::which( 
        colnames(sim_list) != "Weight"  )]
    
    sim_list_stats = sim_list_stats[ sim_list_stats$Ref_Gen == ref_gen,]
    
    if ( dim(sim_list_stats)[1] == 0 )
        warning( paste( "Warning: Identification might be spurious due to 
                        low amount of training sample. No cancer cell line data stored 
                        at this point for reference genome:", ref_gen, sep =" " ) )
    
    
    for (vcf_input_file in vcf_input_files){
        
        if ( ! file.exists(vcf_input_file)){
            warning(paste0("Skipping file, because could not find: ", vcf_input_file))
            next
        } else {
            print(paste0("Found following file, parsing: ",vcf_input_file))
        }
        
        vcf_fingerprint = add_single_file(
            vcf_file_path = vcf_input_file,
            library = library,
            sim_list = sim_list,
            sim_list_stats = sim_list_stats,
            n_threads = n_threads
        )
        
        sim_list = rbind( sim_list, vcf_fingerprint)

    }
    
    print("Aggregating over parsed Cancer Cell Line data")
    
    res_vec = re_calculate_cl_weights( 
        sim_list = sim_list, 
        ref_gen = ref_gen
    )
    
    print("Finished aggregating, saving to database")
    
    if (! test_mode){
    
            write_data_to_db( 
                content_table = res_vec[1], 
                "sim_list",       
                ref_gen = ref_gen, 
                overwrite = TRUE,
                test_mode = test_mode 
            )
            
            write_data_to_db( 
                content_table = res_vec[2], 
                "sim_list_stats", 
                ref_gen = ref_gen,
                overwrite = TRUE, 
                test_mode = test_mode 
            )
    }
    
    print("Finished adding CCLs")
}