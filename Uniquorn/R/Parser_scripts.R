#' parse_cosmic_genotype_data
#' 
#' Parses cosmic genotype data
#' 
#' @param sim_list Variable containing mutations & cell line
#' @param cosmic_file Path to cosmic clp file in hard disk
#' @return The R Table sim_list which contains the CoSMIC CLP fingerprints 
parse_cosmic_genotype_data = function(cosmic_file, sim_list){

    # Only read in columns specified with subset
    subset = c(5, 24)
    
    if (!grepl("CosmicCLP_MutantExport", cosmic_file)){ # MutantExport
        stop("Warning. This is not the recommended COSMIC genotype file!",
             " The recommended file is the 'CosmicCLP_MutantExport.tsv.gz' file.")
        subset = c(5, 19)
    }
  
    cosmic_genotype_tab = data.table::fread(cosmic_file, select = subset, 
                                            sep = "\t", showProgress = FALSE)
    colnames(cosmic_genotype_tab) = c("sample", "position")
    
    # Extract and process coordinates and CL IDs
    coords = cosmic_genotype_tab[, gsub(":|-", "_", position)]
    cls = cosmic_genotype_tab[, gsub("/|(|])| ", "", sample, ignore.case = TRUE)]
    cls[cls == "KM-H2"] = "KMH2"
    cls[cls == "KMH-2"] = "KMH2ALTERNATIVE"
    cls = paste(cls, "COSMIC", sep = "_")
    
    # Build and arrange new sim list with fingerprint and CL ID
    new_sim_list = unique(data.table(Fingerprint = coords, CL = cls, stringsAsFactors = FALSE))
    new_sim_list = new_sim_list[, .(Fingerprint), by = CL]
    new_sim_list = setcolorder(new_sim_list, c("Fingerprint", "CL"))

    sim_list = rbind(sim_list, new_sim_list)
    return(sim_list)
}

#' parse_ccle_genotype_data
#' 
#' Parses ccle genotype data
#' 
#' @param sim_list Variable containing mutations and cell line
#' @param ccle_file Path to CCLE file on hard disk
#' @return The R Table sim_list which contains the CCLE fingerprints
parse_ccle_genotype_data = function(ccle_file, sim_list){
  
    # Only read in columns specified with subset
    subset = c(5, 6, 7, 16)
    ccle_genotype_tab = data.table::fread(ccle_file, select = subset,
                                          sep = "\t", showProgress = FALSE)
    
    # Extract and process coordinates and CL IDs
    coords = paste(ccle_genotype_tab[,Chromosome], 
                   ccle_genotype_tab[,Start_position],
                   ccle_genotype_tab[,End_position], sep = "_")
    
    cls = ccle_genotype_tab[, gsub("\\_.*", "", Tumor_Sample_Barcode)]
    cls = paste(cls, "CCLE", sep = "_")
    
    # Build and arrange new sim list with fingerprint and CL ID
    new_sim_list = unique(data.table(Fingerprint = coords, CL = cls, stringsAsFactors = FALSE))
    new_sim_list = new_sim_list[, .(Fingerprint), by = CL]
    new_sim_list = setcolorder(new_sim_list, c("Fingerprint", "CL"))
    
    sim_list = rbind(sim_list, new_sim_list)
    return(sim_list)
}

#' Filter Parsed VCF Files
#' 
#' Intern utility function. Filters the parsed VCF file for all informations
#' except for the start and length of variations/mutations.
#' 
#' @param vcf_file_path character string giving the path to the vcf file
#'  on the operating system.
#' @usage 
#' parse_vcf_file(vcf_file_path)
#' @return Loci-based DNA-mutational fingerprint of the cancer cell line
#'  as found in the input VCF file.
parse_vcf_file = function(vcf_file_path){
        
      vcf_handle = data.table::fread(paste0("grep -v ^# ", vcf_file_path),
                                     sep = "\t", select = c(1,2,4,5),
                                     fill = TRUE)
      
      # split variants with many alt seq in separate rows
      vcf_handle = vcf_handle[, list(V5 = unlist(strsplit(V5, ","))),
                              by = .(V1, V2, V4)]
      
      # process variants
      chr = vcf_handle[, gsub("chr", "", V1, ignore.case = TRUE)]
      start_var = vcf_handle[, V2]
      length_ref = vcf_handle[, nchar(V4)]
      length_var = vcf_handle[, nchar(V5)]
      length_max = pmax(length_ref, length_var)
      end_var = start_var + (length_max - 1)
      
      # build fingerprint and return
      fingerprint = paste(chr, start_var, end_var, sep = "_")
      return(unique(fingerprint))
}