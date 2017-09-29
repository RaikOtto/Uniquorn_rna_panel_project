#' parse_cosmic_genotype_data
#' 
#' Parses cosmic genotype data
#' 
#' @param cosmic_file Path to cosmic clp file in hard disk
#' @return The R Table sim_list which contains the CoSMIC CLP fingerprints 
parse_cosmic_genotype_data = function(cosmic_file, ref_gen = "GRCH37"){

    # Only read in columns specified with subset
    library_name = "COSMIC"
    package_path = system.file("", package = "Uniquorn")
    rdata_path = paste( c( package_path,"/",library_name,"_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    library_path =  paste( c( package_path,"/Libraries_Ref_gen_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    
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
    print("Parsing Cosmic Coordinates, that might take some time")
    coords = cosmic_genotype_tab[, gsub(":|-", "_", position)]
    seq_name = as.character(sapply( coords, FUN = function(vec){return(as.character(unlist(str_split(vec,"_")))[1])}))
    starts = as.character(sapply( coords, FUN = function(vec){return(as.character(unlist(str_split(vec,"_")))[2])}))
    ends = as.character(sapply( coords, FUN = function(vec){return(as.character(unlist(str_split(vec,"_")))[3])}))
    
    cls = cosmic_genotype_tab[, gsub("/|(|])| ", "", sample, ignore.case = TRUE)]
    cls[cls == "KM-H2"] = "KMH2"
    cls[cls == "KMH-2"] = "KMH2ALTERNATIVE"
    cls = paste(cls, "COSMIC", sep = "_")
    
    c_matches = match(coords, unique(coords), nomatch = 0)
    
    print("Aggregating Cosmic CCL names")
    new_cls = aggregate( cls, by = list(c_matches), FUN= paste, sep = "", collapse= ",")
    
    # Extract and process coordinates and CL IDs
    g_query = GenomicRanges::GRanges(
      seqnames = seq_name,
      IRanges(
        start = as.integer( starts ),
        end = as.integer( ends )
      )
    )
    g_query = unique(g_query)
    mcols(g_query)$Member_CCLs = new_cls$x
    
    if (file.exists(library_path)){
      library_names = readRDS(library_path)
      library_names = unique(c(library_names, "COSMIC"))
      saveRDS(library_names,library_path)
    } else {
      saveRDS("COSMIC",library_path)
    }
    saveRDS(g_query, rdata_path)
    print("Finished parsing Cosmic")
}

#' parse_ccle_genotype_data
#' 
#' Parses ccle genotype data
#' 
#' @param ccle_file Path to CCLE file on hard disk
#' @return The R Table sim_list which contains the CCLE fingerprints
parse_ccle_genotype_data = function(ccle_file, ref_gen = "GRCH37"){
  
    library_name = "CCLE"
    rdata_path = paste( c( package_path,"/",library_name,"_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    
    # Only read in columns specified with subset
    subset = c(5, 6, 7, 16)
    ccle_genotype_tab = data.table::fread(
        ccle_file,
        select = subset,
        sep = "\t",
        showProgress = FALSE
    )
    
    cls = ccle_genotype_tab[, gsub("\\_.*", "", Tumor_Sample_Barcode)]
    cls = paste(cls, "CCLE", sep = "_")
    
    coords    = paste(ccle_genotype_tab$Chromosome,ccle_genotype_tab$Start_position,ccle_genotype_tab$End_position,sep="_")
    c_matches = match(coords,unique(coords), nomatch = 0)
    new_cls <<- c()
    
    print("Aggregating CCLE CCL names")
    new_cls = aggregate( cls, by = list(c_matches), FUN= paste, sep = "", collapse= ",")
    
    # Extract and process coordinates and CL IDs
    g_query = GenomicRanges::GRanges(
        seqnames = paste( "chr", ccle_genotype_tab$Chromosome, sep ="" ),
        IRanges(
            start = ccle_genotype_tab$Start_position,
            end = ccle_genotype_tab$End_position
        )
    )
    g_query = unique(g_query)
    mcols(g_query)$Member_CCLs = new_cls$x
    
    if (file.exists(library_path)){
        library_names = readRDS(library_path)
        library_names = unique(c(library_names, "CCLE"))
        saveRDS(library_names,library_path)
    } else {
        saveRDS("CCLE",library_path)
    }
    saveRDS(g_query, rdata_path)
    print("Finished parsing CCLE")
}

#' subset_dbs_into_mutational_weights
#' 
#' Splits a new db into the mutational weights and creates new subfolders
#' 
#' @param db_g_range GenomicRanges file
#' @param library_name Name of the new library
#' 
#' @return The R Table sim_list which contains the CoSMIC CLP fingerprints 
subset_dbs_into_mutational_weights = function(
    db_g_range,
    library_name,
    ref_gen,
    mutational_weight_inclusion_threshold = c(1,.5,.25)
){
  
  
  
}

#' Filter Parsed VCF Files
#' 
#' Intern utility function. Filters the parsed VCF file for all informations
#' except for the start and length of variations/mutations.
#' 
#' @param vcf_file character string giving the path to the vcf file
#'  on the operating system.
#' @usage 
#' parse_vcf_file(vcf_file)
#' @import GenomicRanges IRanges
#' @return Loci-based DNA-mutational fingerprint of the cancer cell line
#'  as found in the input VCF file.identify_vcf_files
parse_vcf_file = function(vcf_file){
        
      vcf_handle = data.table::fread(paste0("grep -v ^# ", vcf_file),
                                     sep = "\t", select = c(1,2,4,5),
                                     fill = TRUE)
      
      # split variants with many alt seq in separate rows
      vcf_handle = vcf_handle[, list(V5 = unlist(strsplit(V5, ","))),
                              by = .(V1, V2, V4)]
      
      # process variants
      chroms = vcf_handle[, gsub("chr", "", V1, ignore.case = TRUE)]
      start_var = vcf_handle[, V2]
      length_ref = vcf_handle[, nchar(V4)]
      length_var = vcf_handle[, nchar(V5)]
      length_max = pmax(length_ref, length_var)
      end_var = start_var + (length_max - 1)
      
      chroms_pure = grep(chroms, pattern = "_", invert = T)
      chroms      = chroms[ chroms_pure ]
      start_var   = start_var[chroms_pure]
      end_var     = end_var[chroms_pure]
      chroms      = paste("chr",chroms, sep ="")
      
      # build fingerprint and return
      g_query = GenomicRanges::GRanges(
          seqnames = chroms,
          IRanges(
              start = start_var,
              end = end_var
          )
      )
      
      cl_id = gsub("^.*/", "", vcf_file)
      cl_id = gsub(".vcf", "", cl_id, fixed = TRUE)
      cl_id = gsub(".hg19", "", cl_id, fixed = TRUE)
      cl_id = toupper(cl_id)
      cl_id = stringr::str_replace_all(cl_id, pattern = "\\.", "_")
      mcols( g_query )$Member_CCLs = rep(cl_id,length(g_query@seqnames ))
      
      return(g_query)
}