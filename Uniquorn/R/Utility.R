write_mutation_grange_objects = function(
    g_mat,
    library_name,
    mutational_weight_inclusion_threshold,
    ref_gen,
    type = ""
){
  
    package_path = system.file("", package = "Uniquorn")
  
    rdata_path =  paste( c(
        package_path,"/Libraries/",
        ref_gen,"/",library_name,
        "/W",as.character(mutational_weight_inclusion_threshold),
        "_Uniquorn_DB.RData"),
        sep ="", collapse= ""
    )
    
    if ( type != ""){
        rdata_path = str_replace(rdata_path, pattern = "_Uniquorn_DB.RData",
            paste( c( ".",type,"_Uniquorn_DB.RData"), sep = "",collapse = "" )
        )
    }
  
    dir.create(
        paste( 
            c(package_path,"/Libraries/",ref_gen,"/"),
            collapse = "",
            sep =""
        ),
        showWarnings = FALSE)
    dir.create(
        paste( 
            c(package_path,"/Libraries/",ref_gen,"/",library_name),
            collapse = "",
            sep =""
        ),
        showWarnings = FALSE
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
        message(paste(c("Database not found: ", rdata_path,
                        ", creating new database."),collapse = "",sep =""))
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
          showWarnings = FALSE
        )
        dir.create(
            paste( 
                c(package_path,"/Libraries/",ref_gen,"/"),
                collapse = "",
                sep =""
            ),
        showWarnings = FALSE)
        dir.create(
            paste( 
                c(package_path,"/Libraries/",ref_gen,"/",library_name),
                collapse = "",
                sep =""
            ),
        showWarnings = FALSE)
    } else {
        g_mat = readRDS(rdata_path)
    }
    
    return(g_mat)
}
