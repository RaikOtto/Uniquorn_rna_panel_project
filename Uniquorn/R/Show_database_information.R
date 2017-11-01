#' Display Contained Cancer Cell Lines
#' 
#' This function displays the names, amount of mutations/variations and the overall weight of the mutations of all contained cancer cell line fingerprints 
#' for a chosen reference genome and optional library.
#' 
#' @param ref_gen a character vector specifying the reference genome version. All training sets are associated with a reference genome version. Default is \code{"GRCH37"}.
#' @param verbose Should DB informations be printed
#' @return R table which contains the identifiers of all cancer cell line samples
#'  which match the specified parameters (reference genome and library).
#' @usage 
#' show_contained_cls(ref_gen, library = NULL)
#' @examples
#' ##Show all contained cancer cell lines for reference GRCH37
#' show_contained_cls(ref_gen = "GRCH37", verbose = TRUE)
#' 
#' ##Show just cancer cell lines contained in library CELLMINER
#' show_contained_cls(ref_gen = "GRCH37", library = "CELLMINER")
#' @import GenomicRanges stringr
#' @export
show_contained_cls = function(
    ref_gen = "GRCH37",
    verbose = TRUE
){
    
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste( c( package_path,"/Libraries/",ref_gen),
        sep ="", collapse= "")
  
    if ( ! dir.exists( library_path ))
        stop("No libraries found!")
    
    libraries = list.dirs(library_path,full.names = F)
    libraries = libraries[ libraries != ""]
    
    ccls_all <<- data.frame(
        "CCL" = as.character(),
        "Library" = as.character(),
        "W0" = as.character(),
        "W25" = as.character(),
        "W05" = as.character(),
        "W1" = as.character()
    )
    
    for (library_name in libraries){
       
        stats_path = paste( c( library_path,"/",library_name,
            "/CCL_List_Uniquorn_DB.RData"), sep ="", collapse= "")
        g_library = readRDS(stats_path)
        cl_id_no_library = str_replace_all( g_library$CCL, pattern = 
            paste("_",library_name,sep = ""),"" )
        g_library$CCL = cl_id_no_library
        
        if (verbose)
            print(paste(
                c(library_name," amount CCLs: ", as.character(nrow(g_library))),
                sep = "", collapse = "")
        )
        g_library$Library = rep(library_name, nrow(g_library))
        g_library = g_library[c("CCL","Library","W0","W25","W05","W1")]
        ccls_all <<- rbind(ccls_all, g_library )
    }
    return(ccls_all)
}

#' All variants contained in reference library
#' 
#' This function shows all variants contained in a reference library
#' for a given inclusion weight. Default inclusion weight is 0 
#' (all variants)
#' 
#' @param ref_gen a character vector specifying the reference genome 
#' version. All training sets are associated with a reference genome version. 
#' Default is \code{"GRCH37"}.
#' @param library_name Name of the reference library
#' @param mutational_weight_inclusion_threshold Include only mutations 
#' with a weight of at least x. Range: 0.0 to 1.0. 1= unique to CL. 
#' ~0 = found in many CL samples. 
#' @usage 
#' show_contained_variants_in_library(ref_gen)
#' @examples
#' show_contained_variants_in_library(
#' ref_gen = "GRCH37",
#' library_name = "CELLMINER",
#' mutational_weight_inclusion_threshold = 0)
#' @return Returns a GenomicRanges object that contains the variants
#' @export
show_contained_variants_in_library = function(
    ref_gen = "GRCH37",
    library_name,
    mutational_weight_inclusion_threshold = 0
){
  
    message(paste( 
      c("Entered reference genome: ", ref_gen, 
        " entered library name: ", library_name),
        collapse= "",
        sep = "")
    )
    
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste( c( package_path,"/Libraries/",ref_gen), sep ="", collapse= "")
  
    if ( ! dir.exists( library_path ))
        stop("No libraries found!")
  
    libraries = list.dirs(library_path,full.names = F)
    libraries = libraries[ libraries != ""]
    
    if (! (library_name %in% libraries ))
        stop("Could not find library")
    
    g_mat = read_mutation_grange_objects(
        library_name = library_name,
        ref_gen = ref_gen,
        mutational_weight_inclusion_threshold = 
          mutational_weight_inclusion_threshold
    )
  
    return(g_mat)
}

#' Variants In Cancer Cell Line
#' 
#' This function shows all mutations present in the database 
#' for a selected cancer cell line and reference genome.
#' 
#' @param name_cl a character vector giving the identifier 
#' of the cancer cell line for which mutations will be shown.
#' @param ref_gen a character vector specifying the reference 
#' genome version. All training sets are associated with a 
#' reference genome version. Default is \code{"GRCH37"}.
#' @param library_name Name of the reference library 
#' @param mutational_weight_inclusion_threshold Include only mutations 
#' with a weight of at least x. Range: 0.0 to 1.0. 1= unique to CL. 
#' ~0 = found in many CL samples. 
#' @import stringr
#' @usage 
#' show_contained_variants_for_cl(name_cl, 
#' ref_gen,
#' library_name,
#' mutational_weight_inclusion_threshold)
#' @examples 
#' show_contained_variants_for_cl(
#' name_cl = "SK_OV_3",
#' ref_gen = "GRCH37",
#' library_name = "CELLMINER",
#' mutational_weight_inclusion_threshold = 0)
#' @return R table which contains all mutations associated with the 
#' defined cancer cell line and reference genome.
#' @export
show_contained_variants_for_cl = function(
    name_cl,
    ref_gen = "GRCH37",
    library_name,
    mutational_weight_inclusion_threshold = 0
){

    message(paste( 
      c("Entered reference genome: ", ref_gen, 
        " entered library name: ", library_name,
        " entered name of ccl: ", name_cl
      ),
      collapse= "",
      sep = "")
    )
    
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste( c( package_path,"/Libraries/",ref_gen), sep ="", collapse= "")
    
    if ( ! dir.exists( library_path ))
      stop("No libraries found!")
    
    libraries = list.dirs(library_path,full.names = F)
    libraries = libraries[ libraries != ""]
    
    if (! (library_name %in% libraries ))
      stop("Could not find library")
    
    g_mat = read_mutation_grange_objects(
      library_name = library_name,
      ref_gen = ref_gen,
      mutational_weight_inclusion_threshold = 
        mutational_weight_inclusion_threshold
    )
    
    g_query_index = stringr::str_detect(g_mat$Member_CCLs, pattern = name_cl)
    g_query = g_mat[g_query_index,]
    
    return(g_query)
}

#' Cancer Cell Lines With Specific Mutation
#' 
#' This function displays all cancer cell lines in the database which 
#' contain a specified mutation. Closed interval coordinates. 
#' Format mutation: CHR_START_STOP, e.g. 1_123_123.
#' 
#' @param mutation_name a character vector giving the name of the mutation 
#' in the format CHROMOSOME_START_STOP, e.g. \code{"11_244501_244510"}.
#' @param ref_gen a character vector specifying the reference genome 
#' version. All training sets are associated with a reference genome 
#' version. Default is \code{"GRCH37"}.
#' @usage 
#' show_which_cls_contain_mutation(mutation_name, ref_gen)
#' @examples 
#' show_which_cls_contain_mutation(mutation_name = "10_103354427_103354427",
#'                                 ref_gen = "GRCH37")
#' @import DBI 
#' @return R table which contains all cancer cell line samples which 
#' contain the specified mutation with respect to the specified 
#' reference genome version.
#' @export
show_which_cls_contain_mutation = function(mutation_name, ref_gen = "GRCH37"){
  
    message("Reference genome: ", ref_gen)
    
    sim_list = initiate_db_and_load_data(request_table = "sim_list",
                                         subset = c("FINGERPRINT", "CL"),
                                         ref_gen = ref_gen)
    sim_list = sim_list[Fingerprint %like% mutation_name]

    if (nrow(sim_list) == 0){
        message("Could not find any cancer cell line for the mutation ",
                mutation_name, " in the database.")
    } else {
        message("Found ", nrow(sim_list),
                " many cancer cell lines for mutation ", mutation_name,
                " for reference genome ", ref_gen, ".")
    }
    return(sim_list)  
}


read_library_names = function(
    ref_gen
){
  
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste( c( package_path,"/Libraries/",ref_gen,"/"), 
          sep ="", collapse= "")
    
    library_names = list.dirs(library_path, full.names = F)
    library_names = library_names[library_names!= ""]
    
    return(library_names)
}