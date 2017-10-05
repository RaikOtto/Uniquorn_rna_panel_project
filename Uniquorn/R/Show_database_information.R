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
    library_path =  paste( c( package_path,"/Libraries/",ref_gen), sep ="", collapse= "")
  
    if ( ! dir.exists( library_path ))
        stop("No libraries found!")
    
    libraries = list.dirs(library_path,full.names = F)
    libraries = libraries[ libraries != ""]
    
    ccls_all <<- c()
    
    for (library_name in libraries){
        
        stats_path = paste( c( library_path,"/",library_name,"/CCL_List_Uniquorn_DB.RData"), sep ="", collapse= "")
        g_library = readRDS(stats_path)
        
        if (verbose)
            print(paste(
                c(library_name," amount CCLs: ", as.character(nrow(g_library))),
                sep = "", collapse = "")
            )
        ccls_all <<- c(ccls_all,g_library )
    }
    return(ccls_all)
}

#' All Mutations For Reference Genome
#' 
#' This function shows all training-set mutations for a selected reference genome, i.e. the mutations that are being used
#' for identification of a query cancer cell line.
#' 
#' @param ref_gen a character vector specifying the reference genome version. All training sets are associated with a reference genome version. Default is \code{"GRCH37"}.
#' @usage 
#' show_contained_mutations(ref_gen)
#' @examples 
#' show_contained_mutations(ref_gen = "GRCH37")
#' @return R Table which contains all mutations associated with a specified reference genome.
#' @export
show_contained_mutations = function(ref_gen = "GRCH37"){
  
    message("Reference genome: ", ref_gen)
    
    sim_list = initiate_db_and_load_data(request_table = "sim_list",
                                             subset = c("FINGERPRINT", "CL"),
                                             ref_gen = ref_gen)
    
    message("Found ", nrow(sim_list), " many cancer cell lines associated",
            " mutations for reference genome ", ref_gen, ".")
  
    return(sim_list)  
}

#' Mutations In Cancer Cell Line
#' 
#' This function shows all mutations present in the database for a selected cancer cell line and reference genome.
#' 
#' @param name_cl a character vector giving the identifier of the cancer cell line for which mutations will be shown.
#' @param ref_gen a character vector specifying the reference genome version. All training sets are associated with a reference genome version. Default is \code{"GRCH37"}.
#' @import DBI
#' @usage 
#' show_contained_mutations_for_cl(name_cl, 
#' ref_gen)
#' @examples 
#' show_contained_mutations_for_cl(name_cl = "SK_OV_3_CELLMINER_mutations",
#'                                 ref_gen = "GRCH37")
#' @return R table which contains all mutations associated with the defined cancer cell line and reference genome.
#' @export
show_contained_mutations_for_cl = function(name_cl, ref_gen = "GRCH37"){

    message("Reference genome: ", ref_gen)
    
    sim_list = initiate_db_and_load_data(request_table = "sim_list",
                                       subset = c("FINGERPRINT", "CL"),
                                       ref_gen = ref_gen)
    sim_list = sim_list[CL %like% name_cl]

    if (nrow(sim_list) == 0){
        message("Could not find the cancer cell line ", name_cl,
                " in the database.")
    } else {
        message("Found ", nrow(sim_list), " many mutations for cancer cell line",
                name_cl, " for reference genome ", ref_gen )
    }
    return(sim_list)
}

#' Cancer Cell Lines With Specific Mutation
#' 
#' This function displays all cancer cell lines in the database which contain a specified mutation. Closed interval coordinates. Format mutation: CHR_START_STOP, e.g. 1_123_123.
#' 
#' @param mutation_name a character vector giving the name of the mutation in the format CHROMOSOME_START_STOP, e.g. \code{"11_244501_244510"}.
#' @param ref_gen a character vector specifying the reference genome version. All training sets are associated with a reference genome version. Default is \code{"GRCH37"}.
#' @usage 
#' show_which_cls_contain_mutation(mutation_name, ref_gen)
#' @examples 
#' show_which_cls_contain_mutation(mutation_name = "10_103354427_103354427",
#'                                 ref_gen = "GRCH37")
#' @import DBI 
#' @return R table which contains all cancer cell line samples which contain the specified mutation with respect to the specified reference genome version.
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
