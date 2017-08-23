library("stringr")

inclusion_weight = 1.0
benchmark_type   = "non_regularized"
#benchmark_type   = "regularized"

run_benchmark = function(
  inclusion_weight,
  benchmark_type
){
  
    source("~/Dropbox/PhD/Uniquorn_project/Scripts/utility.R")
    gold_t = read.table( file = "~/Dropbox/PhD/Uniquorn_project/Pub/Goldstandard.tab",sep="\t", header = TRUE)
  
    if ( benchmark_type == "regularized"){
      
        ident_result_files_path = paste( 
            "~/Uniquorn_data//benchmark_vcf_files/ident_files_regularized/",
            paste( as.character( inclusion_weight ), "", sep ="/" ), sep = "/"
        )
        auc_file_path = paste0(
            c("~/Dropbox/PhD/Uniquorn_project/AUC_regularized/auc_",as.character(inclusion_weight),"_relaxed.tab" ),
            collapse = ""
        )
        
    } else if ( benchmark_type == "non_regularized") {
    
        ident_result_files_path = paste( 
            "~/Uniquorn_data//benchmark_vcf_files/ident_files_non_regularized//",
            paste( as.character( inclusion_weight ), "", sep ="/" ), sep = "/"
        )
        auc_file_path = paste0(
            c("~/Dropbox/PhD/Uniquorn_project/AUC_non_regularized/auc_",as.character(inclusion_weight),"_relaxed.tab" ),
            collapse = ""
        )
    }
  
  # pre process
  
  b_files =  list.files(
      ident_result_files_path,
      pattern = ".vcf_uniquorn_ident.tab",
      full.names = T
    )
  
  ## benchmark results positive predictions
  
  build_tables()
 
  parse_identification_data = function( b_file ){
    
    print(c(b_file, as.character(which(b_files %in% b_file))))

    b_table = read.table(b_file, sep ="\t", header = T, stringsAsFactors = FALSE)
    
    if ( grepl( "COSMIC", b_file )  )
      source_file = "COSMIC"
    
    if ( grepl( "CCLE", b_file )  )
      source_file = "CCLE"
    
    if ( grepl( "CELLMINER", b_file )  )
      source_file = "CELLMINER"
    
    if ( grepl( "CUSTOM", b_file )  )
        source_file = "CUSTOM"

    b_name        = tail(unlist(str_split(b_file, pattern = "/")  ),1)
    b_name_source = str_replace( b_name, ".vcf_uniquorn_ident.tab", "" )
    b_name = parse_string( b_name, source_file, query_name = T)
    
    identical_cls = as.character( gold_t$Name_identical[ gold_t[,1] == b_name_source ] )
    identical_cls = as.character( unlist( str_split( identical_cls, "," ) ) )
    
    related_cls   = as.character( gold_t$Related       [ gold_t[,1] == b_name_source ] )
    related_cls   = as.character( unlist( str_split( related_cls, "," ) ) )
    
    match_name_query = rep( b_name, dim(b_table)[1] )
      
    match_name_training = parse_string( as.character( b_table$CL ), "" )

    same_identity       = paste( 
        as.character( b_table$CL ), 
        b_table$CL_source, 
        sep = "_" 
    ) %in% identical_cls
    
    related_identity       = paste( 
      as.character( b_table$CL ), 
      b_table$CL_source, 
      sep = "_" 
    ) %in% related_cls
    
    to_be_found = as.character( same_identity | related_identity )
    
    auc_table <<- rbind(
      auc_table, 
      cbind(
        as.character( b_table$Conf_score ),
        as.character( b_table$Found_muts),
        as.character( to_be_found)
      ) 
    )
    colnames( auc_table) = c("Conf_score","Mutations_found","Should_be_found")
    
    write.table( auc_table, auc_file_path, col.names = T, row.names = F, quote = F)
  }
  sapply( b_files, FUN = parse_identification_data)
  
  colnames( auc_table) = c("Conf_score","Mutations_found","Should_be_found")
  
  write.table( auc_table, auc_file_path, col.names = T, row.names = F, quote = F, sep = "\t")
  
}

run_benchmark(
  inclusion_weight = inclusion_weight,
  benchmark_type = benchmark_type
)
