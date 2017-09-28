#' calculate_similarity_results
#' 
#' @param sim_list Contains reference mutation data
#' @param sim_list_stats Contains global reference mutation stats
#' @param found_mut_mapping Mapping to mutations from query to reference mutation set
#' @param minimum_matching_mutations Minimal amount of required matching mutations
#' @param p_value Required maximal p-value
#' @param q_value Required maximal q-value
#' @param confidence_score Threshold above which a positive prediction occurs
#' default 3.0
#' @param vcf_fingerprint The start and end positions of variants in the query
#' @param panels The reference libraries
#' @import data.table
#' @return Results table
calculate_similarity_results = function(
    sim_list,
    sim_list_stats,
    found_mut_mapping,
    minimum_matching_mutations,
    p_value,
    q_value,
    confidence_score,
    vcf_fingerprint,
    panels
){
    
    if (length(found_mut_mapping) != 0){
        # Extract CLs with found mutations and compute abs and rel mutation hit
        candidate_hits_abs = sim_list[, .(CL = CL[found_mut_mapping])][, .(Count=.N), by = CL]
        candidate_hits_abs_all = merge(sim_list_stats[,1], candidate_hits_abs,
                                       by = "CL", all = TRUE)$Count
        candidate_hits_abs_all[is.na(candidate_hits_abs_all)] = 0
        candidate_hits_rel = round(candidate_hits_abs_all / 
                                   sim_list_stats$Count * 100, 1)
        cl_absolute_mutation_hits = sim_list_stats$Count
    }
    
    # p and q-value calculation 
    p_values = calculate_p_and_q_values(
        candidate_hits_abs_all,
        cl_absolute_mutation_hits,
        sim_list,
        sim_list_stats,
        minimum_matching_mutations,
        p_value,
        q_value,
        vcf_fingerprint,
        panels = panels
    )
    q_values       = stats::p.adjust(p_values, "BH")
    conf_score_vec = round(-log(q_values), 2)
    conf_score_vec[is.infinite(conf_score_vec)] = 100
    
    # threshold non weighted
    passed_threshold_vec_p_value = rep(FALSE, length(p_values))
    passed_threshold_vec_p_value[p_values <= p_value] = TRUE
    passed_threshold_vec_q_value = rep(FALSE, length(q_values))
    passed_threshold_vec_q_value[q_values <= q_value] = TRUE
    
    output_cl_names = gsub(paste(panels, collapse  = "|"), "", sim_list_stats$CL)
    panel_vec = sim_list_stats[, gsub("^.*_", "" , unique(CL))]
    
    passed_threshold_vec_p_value[ 
        as.integer(candidate_hits_abs_all)  < 
        as.integer(minimum_matching_mutations) 
    ] = FALSE
    
    passed_threshold_vec_q_value[ 
        as.integer(candidate_hits_abs_all)  < 
        as.integer(minimum_matching_mutations) 
    ] = FALSE
    
    passed_threshold_vec_con_score =
        as.double(conf_score_vec) >=
        as.double(confidence_score)
    
    passed_threshold_vec_con_score[
        as.integer(candidate_hits_abs_all)  < 
        as.integer(minimum_matching_mutations) 
    ] = FALSE
    
    res_table = data.frame(
        "CL"                        = output_cl_names,
        "CL_source"                 = panel_vec,
        "Found_muts"                = as.character(
            candidate_hits_abs_all
        ),
        "Count_mutations"           = as.character(
            cl_absolute_mutation_hits
        ),
        "P_values"                  = as.character(
            p_values
        ),
        "Q_values"                  = as.character(
            q_values
        ),
        "Conf_score"                = as.character(
            conf_score_vec
        ),
        "P_value_sig"               = as.character(
            passed_threshold_vec_p_value
        ),
        "Q_value_sig"               = as.character(
            passed_threshold_vec_q_value
        ),
        "Conf_score_sig"            = as.character( 
            passed_threshold_vec_con_score
        ),
        stringsAsFactors = FALSE
    )
    
    res_table = res_table[ 
        order( 
            as.integer(
                as.character(
                    res_table$Found_muts
                )
            ),
            decreasing = TRUE
        ),
    ]
    
    return(res_table)
}