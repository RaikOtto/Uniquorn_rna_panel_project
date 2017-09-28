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

    # Parse CoSMIC file
    if (file.exists(cosmic_file)){
        message("Found CoSMIC: ", file.exists(cosmic_file))
      
        if (grepl(".gz$", cosmic_file, ignore.case = TRUE)){
            gunzip(cosmic_file, overwrite = TRUE)
            cosmic_file = gsub(".gz$", "", cosmic_file, ignore.case = TRUE)
        }
      
        parse_cosmic_genotype_data( cosmic_file, ref_gen = ref_gen )
    }

    # Parse CCLE file
    if (file.exists(ccle_file)){
      message("Found CCLE: ", file.exists(ccle_file))
      parse_ccle_genotype_data(ccle_file, ref_gen = ref_gen)
    }
    
    if ((!file.exists(cosmic_file)) & (!file.exists(ccle_file))){ 
        warning("Did neither find CCLE & CoSMIC CLP file! Aborting.")
    } else {
        message("Finished parsing, aggregating over parsed Cancer Cell Line data")
        
        message("Finished aggregating, saving to database")
        
        message("Initialization of Uniquorn DB finished")
    }
}