match_query_ccl_to_database = function(
    g_query,
    ref_gen = "GRCH37",
    library
){
  
    package_path = system.file("", package = "Uniquorn")
    rdata_path = paste( c( package_path,"/",library,"_",ref_gen,"_Uniquorn_DB.RData"), sep ="", collapse= "")
    
    message(paste(c("Sample: ",cl_id,", Library: ",library),collapse = "", sep =""))
    
    try( expr = "g_mat = readRDS(rdata_path)")
    
    if (! exists("gut_mat"))
      stop(paste("Cannot find database:", rdata_path, sep =" "))
      
    fo_query = findOverlaps(
        query = g_query,
        subject = g_mat,
        select = "arbitrary"
    )
    mcols(g_mat)$Member_CCLs[fo_query]
    
}