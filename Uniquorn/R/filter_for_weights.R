#' filter_for_weights
#' 
#' Filter the reference set
#' 
#' \code{filter_for_weights} parses vcf file and output 
#' basic information
#' @param ref_gen Reference genome version. All training sets are 
#' associated with a reference genome version. Default: GRCH37
#' @param verbose Print additional information
#' @param sim_list Contains the mutations
#' @param sim_list_stats Contains the overal mutation statistics
#' @param mutational_weight_inclusion_threshold Lower bound for mutational weight to be 
#' included
#' @usage 
#' filter_for_weights( 
#'     mutational_weight_inclusion_threshold,
#'     ref_gen,
#'     verbose,
#'     sim_list,
#'     sim_list_stats)
#' @import data.table
#' @return Filtered reference sets
filter_for_weights = function( 
    mutational_weight_inclusion_threshold,
    ref_gen,
    verbose,
    sim_list,
    sim_list_stats
    ){
    
    # filter for weights
    if (mutational_weight_inclusion_threshold != 0.0){
        if (verbose){    
            message("Adjusted mutational inclusion weight, only using mutations",
                    " that have a weight higher than ", 
                    mutational_weight_inclusion_threshold, ".")
        }
        # Exlcude mutations and recalculate weights
        sim_list = sim_list[Weight >= mutational_weight_inclusion_threshold]
        sim_list_stats = sim_list[, .(Count=.N), by = CL]
        all_weights = sim_list[, .(All_weights=sum(Weight)), by = CL]
        data.table::setkey(sim_list_stats, CL)
        data.table::setkey(all_weights, CL)
        sim_list_stats = sim_list_stats[all_weights]
        sim_list_stats = sim_list_stats[, .(CL, Count, round(All_weights, 1))]
        
        if (verbose){
          message("Found ", nrow(sim_list), " many mutations with mutational ",
                  "weight of at least ", mutational_weight_inclusion_threshold, ".")
        }
        res_list = list("sim_list" = sim_list, "sim_list_stats" = sim_list_stats)
        return(res_list)
    }
}