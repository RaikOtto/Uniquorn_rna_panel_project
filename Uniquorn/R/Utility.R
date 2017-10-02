read_ccl_list = function(
    library_name,
    ref_gen
){
  
  package_path = system.file("", package = "Uniquorn")
  
  CCL_List_data_path =  paste( c(
    package_path,"/Libraries/",
    ref_gen,"/",library_name,
    "/CCL_List_Uniquorn_DB.RData"),
    sep ="", collapse= ""
  )
  
  return(readRDS(CCL_List_data_path))
}

write_ccl_list = function(
  ccl_list,
  library_name,
  ref_gen
){
  
    package_path = system.file("", package = "Uniquorn")
    
    CCL_List_data_path =  paste( c(
      package_path,"/Libraries/",
      ref_gen,"/",library_name,
      "/CCL_List_Uniquorn_DB.RData"),
      sep ="", collapse= ""
    )
    saveRDS(ccl_list,CCL_List_data_path)
}

write_mutation_grange_objects = function(
    g_mat,
    library_name,
    mutational_weight_inclusion_threshold,
    ref_gen
){
  
    package_path = system.file("", package = "Uniquorn")
  
    rdata_path =  paste( c(
        package_path,"/Libraries/",
        ref_gen,"/",library_name,
        "/W",as.character(mutational_weight_inclusion_threshold),
        "_Uniquorn_DB.RData"),
        sep ="", collapse= ""
    )
  
    dir.create(
        paste( 
            c(package_path,"/Libraries/",ref_gen,"/"),
            collapse = "",
            sep =""
        ),
        showWarnings = F)
    dir.create(
        paste( 
            c(package_path,"/Libraries/",ref_gen,"/",library_name),
            collapse = "",
            sep =""
        ),
        showWarnings = F
    )
    saveRDS(g_mat,rdata_path)
}

read_mutation_grange_objects = function(
    library_name,
    mutational_weight_inclusion_threshold,
    ref_gen
){
  
    package_path = system.file("", package = "Uniquorn")
    
    rdata_path =  paste( c(
        package_path,"/Libraries/",
        ref_gen,"/",library_name,
        "/W",as.character(mutational_weight_inclusion_threshold),
        "_Uniquorn_DB.RData"),
        sep ="",
        collapse= ""
    )
    
    if (! file.exists( rdata_path )){
        message(paste(c("Database not found: ", rdata_path,", creating new database."),collapse = "",sep =""))
        g_mat = GenomicRanges::GRanges(
            seqnames = NULL,
            IRanges(
                start = NULL,
                end = NULL
            )
        )
        dir.create(
          paste( 
            c(package_path,"/Libraries/"),
            collapse = "",
            sep =""
          ),
          showWarnings = F
        )
        dir.create(
            paste( 
                c(package_path,"/Libraries/",ref_gen,"/"),
                collapse = "",
                sep =""
            ),
        showWarnings = F)
        dir.create(
            paste( 
                c(package_path,"/Libraries/",ref_gen,"/",library_name),
                collapse = "",
                sep =""
            ),
        showWarnings = F)
    } else {
        g_mat = readRDS(rdata_path)
    }
    
    return(g_mat)
}

read_library_names = function(
    ref_gen
){
  
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste( c( package_path,"/Libraries/",ref_gen,"/"), sep ="", collapse= "")
  
    library_names = list.dirs(library_path, full.names = F)
    library_names = library_names[library_names!= ""]
    
    return(library_names)
}
