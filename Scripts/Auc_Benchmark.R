library("stringr")
library("devtools")
library("argparse")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

parser = ArgumentParser()
parser$add_argument('-iw', "--inclusion_weight", type="double")
parser$add_argument('-p', "--panel_mode", action="store_true", default = "FALSE")
parser$add_argument('-r', "--robust_mode", action="store_true", default = "FALSE")
args = parser$parse_args()

inclusion_weight     = args$inclusion_weight
panel_mode           = args$panel_mode
robust_mode          = args$robust_mode

regex_term = "\\(|\\*|\\-|\\-|\\_|\\+|\\-|\\)|\\-|\\:|\\[|\\]|\\."

parse_identification_data = function( b_file, gold_t,auc_file_path,ident_result_files_path ){
  
    print(b_file)
    
    b_table = read.table(b_file, sep ="\t", header = T, stringsAsFactors = FALSE)
    library_names = Uniquorn::read_library_names(ref_gen = "GRCH37")
    
    source_file = as.character(unlist(
      str_split(b_file, pattern = "\\.")
    ))[2]
    
    b_name = tail(unlist(str_split(b_file, pattern = "/")  ),1)
    b_name = str_to_upper(b_name)
    b_name_plane = str_replace( b_name, pattern = paste( c(".",source_file,".IDENT.TSV"), sep ="", collapse = ""),"")
    b_name_plane = str_replace_all( b_name_plane, pattern = regex_term, "")
    b_name_full = str_replace( str_replace( b_name, pattern = ".IDENT.TSV",""), "\\.", "_" )
    query_vec = rep(b_name_full, nrow(b_table))
    
    b_name_ccls = str_to_upper( as.character(b_table$CCL) )
    b_name_ccls_plane = str_to_upper( as.character(sapply( b_name_ccls, FUN = str_replace_all, pattern = regex_term,"")) )
    b_name_ccls_full = paste( as.character(b_name_ccls_plane), b_table$Library, sep ="_" )
    
    ### new
    
    gold_cls = str_to_upper(as.character(gold_t$CL_plane))
    gold_cls = str_replace_all( gold_cls, pattern = regex_term, "")
    
    ident_cls = str_to_upper(as.character(b_name_plane))
    ident_cls = str_replace_all(ident_cls, pattern = regex_term,"")
    identifier_index = which( gold_cls == ident_cls )
    
    identical_cls = str_to_upper( as.character( gold_t$Name_identical[ identifier_index ] ) )
    identical_cls = as.character( unlist( str_split( identical_cls, "," ) ) )
    identical_cls = identical_cls[ identical_cls != ""]
    
    split_ident = lapply( identical_cls, FUN = function(vec){return( unlist(str_split(vec, pattern = "_")  ) )} )
    ident_len = sapply(split_ident, FUN = length)
    ident_lib = as.character( sapply( split_ident, FUN =  tail, 1 ) )
    ident_ident <<- c()
    for(i in 1:length(identical_cls)){
      ident_ident <<- c( ident_ident, paste( c( head( as.character( unlist(split_ident[i])), ident_len[i] - 1)),collapse= "", sep ="" ) )
    }
    ident_ident = str_to_upper( as.character(sapply( ident_ident, FUN = str_replace_all, pattern = regex_term,"")) )
    ident_ident = paste(ident_ident, ident_lib, sep = "_")
    identical_cls_vec = b_name_ccls_full %in% ident_ident
    
    related_cls   = as.character( gold_t$Related[ identifier_index ] )
    related_cls   = as.character( unlist( str_split( related_cls, "," ) ) )
    related_cls   = related_cls[ related_cls != ""]
    
    split_related = lapply( related_cls, FUN = function(vec){return( unlist(str_split(vec, pattern = "_")  ) )} )
    related_len = sapply(split_related, FUN = length)
    related_lib = as.character( sapply( split_related, FUN =  tail, 1 ) )
    related_ident <<- c()
    for(i in 1:length(related_cls)){
      if (length(related_cls) == 0)
        next()
      related_ident <<- c( related_ident, head( as.character( unlist(split_related[i])), related_len[i] - 1) )
    }
    related_ident = str_to_upper( as.character(sapply( related_ident, FUN = str_replace_all, pattern = regex_term,"")) )
    related_ident = paste(related_ident, related_lib, sep = "_")
    
    related_cls_vec = b_name_ccls_full %in% related_ident
    
    identical_related_vcf = identical_cls_vec | related_cls_vec
    
    auc_table_new = data.frame(
      "P_values" = b_table$P_values,
      "Should_be_found" = identical_related_vcf
    )
    
    auc_table <<- rbind(
      auc_table, 
      auc_table_new
    )
    
    write.table( auc_table, auc_file_path, col.names = T, row.names = F, quote = F, sep ="\t")
}

run_benchmark = function(
    inclusion_weight,
    panel_mode,
    robust_mode
){
  
    source("~/Uniquorn_rna_panel_project//Scripts/utility.R")
    gold_t = read.table( file = "~/Uniquorn_rna_panel_project/Misc//Goldstandard.tsv",sep="\t", header = TRUE)
  
    inclusion_weight_path = as.character( inclusion_weight )
    inclusion_weight_path = str_replace( as.character( inclusion_weight_path ), pattern ="\\.", "_" )
    ident_result_files_path = paste( 
        "~/Uniquorn_data//benchmark_vcf_files/ident_files_regularized/",
        paste( as.character( inclusion_weight_path ), "", sep ="/" ), sep = "/"
    )
    auc_file_path = paste0(
        c("~/Uniquorn_data/benchmark_vcf_files/AUC/auc_",as.character(inclusion_weight),".tsv" ),
        collapse = ""
    )
    
    seen_obj_path = paste0(
        c("~/Uniquorn_data/benchmark_vcf_files/AUC/auc_",as.character(inclusion_weight),".RDS" ),
        collapse = ""
    )
        
    # pre process
    
    b_files =  list.files(
        ident_result_files_path,
        pattern = ".ident.tsv",
        full.names = T
    )
    b_files[sapply(b_files, file.size) > 300000]
    
    if (file.exists(seen_obj_path)) {
        seen_obj <<- readRDS(seen_obj_path)
    } else {
        seen_obj <<- c()
    }
        
    
    ## benchmark results positive predictions
  
    build_tables()
 
    for (b_file in b_files){
        
        identifier = str_replace( tail( as.character(unlist(str_split(b_file, pattern = "/"))), 1 ), pattern =".ident.tsv", "" )
        
        if (identifier %in% seen_obj){
            next()
        } else {
            seen_obj <<- c(seen_obj, identifier)
            saveRDS(seen_obj, file = seen_obj_path)
            parse_identification_data(
                b_file,
                gold_t,
                auc_file_path,
                ident_result_files_path
            )
        }
    }

    write.table( auc_table, auc_file_path, col.names = T, row.names = F, quote = F, sep = "\t")
}

run_benchmark(
    inclusion_weight = inclusion_weight,
    panel_mode = panel_mode,
    robust_mode = robust_mode
)
