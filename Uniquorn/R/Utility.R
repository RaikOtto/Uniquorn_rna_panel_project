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
        sep ="", collapse= "")
    
    if (! file.exists( rdata_path )){
        message(paste(c("Database not found: ", rdata_path),collapse = "",sep =""))
        g_mat = GenomicRanges::GRanges(
            seqnames = NULL,
            IRanges(
                start = NULL,
                end = NULL
            )
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


#' initiate_db_and_load_data
#' 
#' Intern utility function, load database and return the sim_list and sim_list_stats tables.
#' 
#' @inheritParams identify_vcf_file
#' @param request_table character string giving the name of the table to be
#'  extracted from the database.
#' @param load_default_db logical indicating whether the default database
#'  should be used as source for the data.
#' @return Returns the sim_list and sim_list_stats tables.
#' @usage 
#' initiate_db_and_load_data(ref_gen, request_table, load_default_db)
#' @import DBI RSQLite
initiate_db_and_load_data = function(ref_gen, request_table, subset, load_default_db = FALSE){
    
    package_path = system.file("", package = "Uniquorn")
    default_database_path =  paste( package_path, "uniquorn_db_default.sqlite", sep ="/")
    database_path = paste(package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/")

    drv = RSQLite::SQLite()
    
    if(file.exists(database_path) & (!load_default_db) & (file.size(database_path) != 0))
        con = DBI::dbConnect(drv, dbname = database_path)
    else
        con = DBI::dbConnect(drv, dbname = default_database_path)

    #res = DBI::dbGetQuery(con, paste0("SELECT ", paste(subset, collapse = ", "),
    #                             " FROM ", request_table, " WHERE Ref_Gen LIKE ",
    #                             paste0("'%", ref_gen, "%'")))
    res = DBI::dbGetQuery(con, paste0("SELECT ", paste(subset, collapse = ", "),
                                 " FROM ", request_table))
    res = data.table::setDT(res)
    
    DBI::dbDisconnect(con)
    
    return(res)
}

#' write_data_to_db
#' 
#' Intern utility function, writes to database the sim_list and sim_list_stats variables
#' 
#' @inheritParams identify_vcf_file
#' @param content_table Tables to be written in db
#' @param table_name Name of the table to be written into the DB
#' @param overwrite Overwrite the potentially existing table
#' @param test_mode Is this a test? Just for internal use 
#' @return the sim_list and sim_list_stats variable
#' @usage 
#' write_data_to_db( 
#'    content_table, 
#'    table_name,
#'    ref_gen,
#'    overwrite,
#'    test_mode )
#' @import DBI RSQLite
write_data_to_db = function( 
    content_table, 
    table_name, 
    ref_gen = "GRCH37",
    overwrite = TRUE, 
    test_mode = FALSE 
    ){
    
    package_path    = system.file("", package = "Uniquorn")
    
    database_path =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )

    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    if( test_mode )
        con = DBI::dbConnect(drv, dbname = "")
    
    DBI::dbWriteTable( con, name = table_name, value = as.data.frame( content_table ), overwrite = overwrite )
    
    DBI::dbDisconnect(con)

}

#' check_files
#' 
#' Intern utility function, checks for existance of user specified vcf files.
#' 
#' @param vcf_input_file a character string giving the path to the VCF file.
#' @return logical; indicating whether the specified file could be located on
#' the system.
check_files = function(vcf_input_file){
    if (!file.exists(vcf_input_file)){
        return(FALSE)
    } else{
        return(TRUE)
    }
}

#' Re-calculate sim_list_weights
#' 
#' This function re-calculates the weights of mutation after a change of the training set
#' @inheritParams write_data_to_db
#' @return A list containing both the sim_list at pos 1 and sim_list_stats at pos 2 data frames.
#' @param sim_list R Table which contains a mapping 
#' from mutations/ variations to their containing CLs
re_calculate_cl_weights = function(sim_list, ref_gen){
    
    panels = sim_list[, unique(gsub("^.*_", "" , unique(CL)))]
    message(paste("Distinguishing between panels:",
        paste0(c(panels), collapse = ", "), sep = " "))
    
    if(!exists("sim_list_global"))
        sim_list_global = data.table()
    if(!exists("sim_list_stats_global"))
      sim_list_stats_global = data.table()
    
    for (panel in panels) {
        
      message(paste0("Current panel: ", panel))
      
      sim_list_panel = subset(sim_list, CL %like% panel)
      sim_list_stats_panel = sim_list_panel[, .(Count=.N), by = CL]
        
      message("Aggregating over mutational frequency to obtain mutational weight")
      weights_panel = sim_list_panel[, .(Weight=.N), by = Fingerprint]
      weights_panel$Weight = weights_panel[, 1.0 / Weight]
      data.table::setkey(weights_panel, Fingerprint)
      data.table::setkey(sim_list_panel, Fingerprint)
      sim_list_panel = sim_list_panel[weights_panel]
        
      # calculate weights
      aggregation_all_panel = sim_list_panel[, .(All_weights=sum(Weight)), by = CL]
      data.table::setkey(aggregation_all_panel, CL)
      data.table::setkey(sim_list_stats_panel, CL)
      sim_list_stats_panel = sim_list_stats_panel[aggregation_all_panel]
    
      # prepare final result tables
      Ref_Gen = rep(ref_gen, nrow(sim_list_panel))
      sim_list_panel = cbind(sim_list_panel, Ref_Gen)
      Ref_Gen = rep(ref_gen, nrow(sim_list_stats_panel))
      sim_list_stats_panel = cbind(sim_list_stats_panel, Ref_Gen)
        
      # store panel results
      sim_list_global = rbind(sim_list_global, sim_list_panel)
      sim_list_stats_global = rbind(sim_list_stats_global, sim_list_stats_panel)
      
    }
    return(list(sim_list_global, sim_list_stats_global))
}
