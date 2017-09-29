match_query_ccl_to_database = function(
    g_query,
    ref_gen = "GRCH37",
    library_name,
    mutational_weight_inclusion_threshold
){
  
    g_mat = read_mutation_grange_objects(
        ref_gen = ref_gen,
        library_name = library_name,
        mutational_weight_inclusion_threshold = mutational_weight_inclusion_threshold
    )
    # IMPLEMENT IMPORT OF CLLs FOR LIBRARY HERE
    
    fo_query = findOverlaps(
        query = g_query,
        subject = g_mat,
        select = "arbitrary",
        type = "equal"
    )
    match = subsetByOverlaps(g_mat, g_query)
    hit_ccls = as.character(unlist(str_split(mcols(match)$Member_CCLs, pattern = "," )))
    
    #sort(table(hit_ccls),decreasing = T)
    return(sort(table(hit_ccls), decreasing = T))
}