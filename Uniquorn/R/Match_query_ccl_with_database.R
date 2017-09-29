match_query_ccl_to_database = function(
    g_query,
    ref_gen = "GRCH37",
    library_name,
    mutational_weight_inclusion_threshold
){
  
    package_path = system.file("", package = "Uniquorn")
    rdata_path = paste( c( package_path,"/",library_name,"_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    
    # IMPLEMENT IMPORT OF CLLs FOR LIBRARY HERE
    
    if (! file.exists(rdata_path))
        stop(paste("DB not found: ",rdata_path))
    
    g_mat = readRDS(rdata_path)

    fo_query = findOverlaps(
        query = g_query,
        subject = g_mat,
        select = "arbitrary",
        type = "equal"
    )
    match = subsetByOverlaps(g_mat, g_query)
    hit_ccls = as.character(unlist(str_split(mcols(match)$Member_CCLs, pattern = "," )))
    
    table(hit_ccls)
    return(hit_ccls)
}