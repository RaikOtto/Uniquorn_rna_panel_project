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
    seq_name = as.character(sapply( coords, FUN = function(vec){
        return(as.character(unlist(str_split(vec,"_")))[1])}))
    starts = as.character(sapply( coords, FUN = function(vec){
        return(as.character(unlist(str_split(vec,"_")))[2])}))
    ends = as.character(sapply( coords, FUN = function(vec){
        return(as.character(unlist(str_split(vec,"_")))[3])}))
    
    cls = cosmic_genotype_tab[, gsub("/|(|])| ", "", sample, ignore.case = TRUE)]
    cls[cls == "KM-H2"] = "KMH2"
    cls[cls == "KMH-2"] = "KMH2ALTERNATIVE"
    
    c_matches = match(coords, unique(coords), nomatch = 0)
    
    print("Aggregating Cosmic CCL names")
    new_cls = aggregate( cls, by = list(c_matches), FUN= paste, sep = "", collapse= ",")
    
    # Extract and process coordinates and CL IDs
    g_mat = GenomicRanges::GRanges(
        seqnames = seq_name,
        IRanges(
            start = as.integer( starts ),
            end = as.integer( ends )
        )
    )
    g_mat = unique(g_mat)
    mcols(g_mat)$Member_CCLs = new_cls$x
    mcols(g_mat)$Member_CCLs = str_replace_all(mcols(g_mat)$Member_CCLs, 
        pattern = paste( "_",library_name,sep =""), "")
    
    print("Normalizing CCL names")
    #n = sapply( mcols(g_mat)$Member_CCLs, FUN = function(vec){
    #    return(
    #        paste(
    #            unique(
    #                as.character(unlist(str_split(mcols(g_mat)$Member_CCLs,
    #                pattern  =",")))
    #            ), collapse = ","
    #        )
    #    )
    #} )
    
    write_mutation_grange_objects(
        g_mat = g_mat,
        library_name = library_name,
        ref_gen = ref_gen, 
        mutational_weight_inclusion_threshold = 0
    )
    
    ccl_stats <<- data.frame(
      "CCL" = names(table(as.character(
        unlist(str_split(new_cls$x, pattern = ","))))),
      "Mut_count_0" = as.character(table(as.character(unlist(
        str_split(new_cls$x, pattern = ",")))))
    )
    
    for(mwit in c(.25,.5,1.0)){
      
      hit_index = which( str_count(
        mcols(g_mat)$Member_CCLs, pattern = ","
      ) <= as.integer( round(1/mwit) ) )
      
      mwit_g_mat = g_mat[hit_index]
      
      write_mutation_grange_objects(
        mutational_weight_inclusion_threshold = mwit,
        g_mat = mwit_g_mat,
        library_name = library_name,
        ref_gen = ref_gen
      )
      
      mut_t = table(as.character(unlist(str_split( 
        mcols(mwit_g_mat)$Member_CCLs,pattern = "," ))))
      
      ccl_count = rep("0", nrow(ccl_stats))
      ccl_match = match( ccl_stats$CCL, names(mut_t), nomatch = 0 )
      ccl_count[ccl_match] = mut_t[ccl_match != 0]
      ccl_stats = cbind(ccl_stats, ccl_count)
    }
    colnames(ccl_stats) = c("CCL","W0","W25","W05","W1")
    
    write_ccl_list(ccl_list = ccl_stats,ref_gen = ref_gen,library_name = library_name)
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
    
    # Only read in columns specified with subset
    subset = c(5, 6, 7, 16)
    ccle_genotype_tab = data.table::fread(
        ccle_file,
        select = subset,
        sep = "\t",
        showProgress = FALSE
    )
    
    cls = ccle_genotype_tab[, gsub("\\_.*", "", Tumor_Sample_Barcode)]
    cls = str_replace(cls, paste( "_",library_name, sep = "" ), "")
    
    coords    = paste(
        str_replace( ccle_genotype_tab$Chromosome, pattern = "^chr", "" ),
        ccle_genotype_tab$Start_position,
        ccle_genotype_tab$End_position,
    sep="_")
    c_matches = match(coords,unique(coords), nomatch = 0)
    new_cls <<- c()
    
    print("Aggregating CCLE CCL names")
    new_cls = aggregate( cls, by = list(c_matches), FUN= paste, sep = "", collapse= ",")
    
    # Extract and process coordinates and CL IDs
    g_mat = GenomicRanges::GRanges(
        seqnames = str_replace( ccle_genotype_tab$Chromosome, pattern = "^chr", "" ),
        IRanges(
            start = ccle_genotype_tab$Start_position,
            end = ccle_genotype_tab$End_position
        )
    )
    g_mat = unique(g_mat)
    mcols(g_mat)$Member_CCLs = new_cls$x
    
    write_mutation_grange_objects(
        g_mat = g_mat,
        library_name = library_name,
        ref_gen = ref_gen, 
        mutational_weight_inclusion_threshold = 0
    )
    
    ccl_stats <<- data.frame(
        "CCL" = names(table(as.character(
            unlist(str_split(new_cls$x, pattern = ","))))),
        "Mut_count_0" = as.character(table(as.character(unlist(
            str_split(new_cls$x, pattern = ",")))))
    )
    
    for(mwit in c(.25,.5,1.0)){
      
        hit_index = which( str_count(
            mcols(g_mat)$Member_CCLs, pattern = ","
        ) <= as.integer( round(1/mwit) ) )
        
        mwit_g_mat = g_mat[hit_index]
        
        write_mutation_grange_objects(
            mutational_weight_inclusion_threshold = mwit,
            g_mat = mwit_g_mat,
            library_name = library_name,
            ref_gen = ref_gen
        )
        
        mut_t = table(as.character(unlist(str_split( 
            mcols(mwit_g_mat)$Member_CCLs,pattern = "," ))))
        
        ccl_count = rep("0", nrow(ccl_stats))
        ccl_match = match( ccl_stats$CCL, names(mut_t), nomatch = 0 )
        ccl_count[ccl_match] = mut_t[ccl_match != 0]
        ccl_stats = cbind(ccl_stats, ccl_count)
    }
    colnames(ccl_stats) = c("CCL","W0","W25","W05","W1")
 
    write_ccl_list(ccl_list = ccl_stats,ref_gen = ref_gen,library_name = library_name)
    
    print("Finished parsing CCLE")
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
      chroms      = str_replace(chroms,pattern = "^chr","")
      start_var   = start_var[chroms_pure]
      end_var     = end_var[chroms_pure]
      
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
