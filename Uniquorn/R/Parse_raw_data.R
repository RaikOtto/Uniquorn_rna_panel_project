#' initiate_canonical_databases
#' 
#' Parses data into r list variable
#' 
#' @param cosmic_file The path to the cosmic DNA genotype data file. 
#' Ensure that the right reference genome is used
#' @param ccle_file The path to the ccle DNA genotype data file. 
#' Ensure that the right reference genome is used
#' @param ref_gen Reference genome version
#' @return Returns message if parsing process has succeeded
#' @import DBI R.utils RSQLite stringr
#' @usage 
#' initiate_canonical_databases(
#' cosmic_file = "CosmicCLP_MutantExport.tsv",
#' ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
#' ref_gen = "GRCH37")
#' @examples 
#' initiate_canonical_databases(
#' cosmic_file = "CosmicCLP_MutantExport.tsv",
#' ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
#' ref_gen = "GRCH37")
#' @export
initiate_canonical_databases = function(
    cosmic_file = "CosmicCLP_MutantExport.tsv",
    ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
    ref_gen = "GRCH37"
    ){

    message("Reference genome: ", ref_gen)

    # pre-processing
    sim_list = initiate_db_and_load_data(request_table = "sim_list",
                                         subset = c("FINGERPRINT", "CL"),
                                         ref_gen = ref_gen, load_default_db = TRUE)

    # Parse CoSMIC file
    if (file.exists(cosmic_file)){
        message("Found CoSMIC: ", file.exists(cosmic_file))
      
        if (grepl(".gz$", cosmic_file, ignore.case = TRUE)){
            gunzip(cosmic_file, overwrite = TRUE)
            cosmic_file = gsub(".gz$", "", cosmic_file, ignore.case = TRUE)
        }
        sim_list = parse_cosmic_genotype_data(cosmic_file, sim_list)
    }
    t_t_index = grep(sim_list$CL, pattern = "T-T_COSMIC")
    sim_list$CL[t_t_index] = "TT2_COSMIC"
  
    # Parse CCLE file
    if (file.exists(ccle_file)){
      message("Found CCLE: ", file.exists(ccle_file))
      sim_list = parse_ccle_genotype_data(ccle_file, sim_list)
    }
    
    if ((!file.exists(cosmic_file)) & (!file.exists(ccle_file))){ 
        warning("Did neither find CCLE & CoSMIC CLP file! Aborting.")
    } else {
        message("Finished parsing, aggregating over parsed Cancer Cell Line data")
        res_vec = re_calculate_cl_weights( 
            sim_list = sim_list, 
            ref_gen = ref_gen 
        )
      
        message("Finished aggregating, saving to database")
        write_data_to_db( 
            content_table = as.data.frame(res_vec[1]), 
            "sim_list",       
            ref_gen = ref_gen,
            overwrite = TRUE 
        )
        write_data_to_db( 
            content_table = as.data.frame(res_vec[2]), 
            "sim_list_stats", 
            ref_gen = ref_gen, 
            overwrite = TRUE 
        )
        message("Initialization of Uniquorn DB finished")
    }
}