library("stringr")
library("BiocParallel")
library("doParallel")
library("foreach")
library("devtools")
<<<<<<< HEAD
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

inclusion_weight     = 1.0
only_first           = FALSE
exclude_self         = FALSE
p_value = .05
minimum_matching_mutations = 0
run_identification   = F
panel_mode = T
auc_mode             = FALSE
ref_gen              = "GRCH37"
#type_benchmark       = "non_regularized"
panel_mode           = FALSE
type_benchmark       = "regularized"
n_threads = 4

run_benchmark = function(
    inclusion_weight = inclusion_weight, 
    only_first = only_first,
    exclude_self = exclude_self,
    run_identification = run_identification, 
    #confidence_score = confidence_score,
    minimum_matching_mutations = minimum_matching_mutations,
    type_benchmark,
    ref_gen
=======
library("argparse")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

parser = ArgumentParser()
parser$add_argument('-iw', "--inclusion_weight", type="double")
parser$add_argument('-p', "--panel_mode", action="store_true", default = "FALSE")
parser$add_argument('-nt', "--number_threads", type="integer", default = "1")
args = parser$parse_args()

inclusion_weight     = args$inclusion_weight
panel_mode           = args$panel_mode

inclusion_weight      = 0.0
panel_mode            = FALSE

p_value               = .05
minimum_matching_mutations = 0
run_identification   = F
auc_mode             = FALSE

run_benchmark = function(
    inclusion_weight = inclusion_weight, 
    minimum_matching_mutations = minimum_matching_mutations
>>>>>>> 52e010a4c1d3331df1913c4c690c119fc97a08c9
){
  
    source("~/Uniquorn_rna_panel_project/Scripts/utility.R")
    
    build_path_variables(
<<<<<<< HEAD
        inclusion_weight = inclusion_weight,
        only_first = only_first,
        exclude_self = exclude_self,
        run_identification = run_identification,
        cellminer = F,
        type_benchmark = type_benchmark
    )
    out_path_ident = str_replace(
        benchmark_ident_file_path,
        pattern = "/Benchmark_results_regularized//0.5_Benchmark_identification_result.tab",
        replacement = paste(c("ident_files_regularized",inclusion_weight,""),sep="",collapse= "/")
    )
=======
      inclusion_weight = inclusion_weight,
      only_first = FALSE,
      exclude_self = FALSE,
      run_identification = FALSE,
      cellminer = FALSE,
      type_benchmark = "regularized"
    )HFHQC-584H7-0I87I
    out_path_ident = str_replace(
      benchmark_ident_file_path,
      pattern = "/Benchmark_results_regularized/0.5_Benchmark_identification_result.tab",
      replacement = paste(c("ident_files_regularized",inclusion_weight,""),sep="",collapse= "/")
    )
>>>>>>> 52e010a4c1d3331df1913c4c690c119fc97a08c9
    gold_t <<- read.table( file = "~//Uniquorn_rna_panel_project/Misc/Goldstandard.tsv",sep="\t", header = TRUE)
    
    if (panel_mode){
        ident_result_files_path <<- str_replace(ident_result_files_path,pattern = "ident_files","panel_ident_files")
        benchmark_ident_file_path <<- str_replace(benchmark_ident_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
        benchmark_res_file_path <<- str_replace(benchmark_res_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
    }
<<<<<<< HEAD

=======
    
>>>>>>> 52e010a4c1d3331df1913c4c690c119fc97a08c9
    b_files <<- list.files(ident_result_files_path , pattern = ".ident.tsv", full.names = T )
    
    ## benchmark results positive predictions
    
    build_tables()
<<<<<<< HEAD
  
    for (b_file in b_files)
        parse_identification_data(b_file)

    #sapply( b_files, FUN = parse_identification_data)
  
    res_table = res_table[ order( as.double( as.character(res_table$Found_muts_abs) ), decreasing = T),  ]
    write.table( res_ident_table, benchmark_ident_file_path, sep ="\t", quote = F, row.names = F  )
    write.table( res_table, benchmark_res_file_path, sep ="\t", quote =F, row.names = F )
}

parse_identification_data = function( b_file ){
  
    print(c(b_file, as.character(which(b_files %in% b_file))))
    
    b_table = read.table(b_file, sep ="\t", header = T, stringsAsFactors = FALSE)
    library_names = Uniquorn::read_library_names(ref_gen = ref_gen)
    
    for (library_name in library_names){
        if ( str_detect(pattern = library_name, b_file) )
            source_file <<- library_name
    }
    
    b_name = tail(unlist(str_split(b_file, pattern = "/")  ),1)
    b_name_plane = str_replace( b_name, pattern = paste( c(".",source_file,".ident.tsv"), sep ="", collapse = ""),"")
    b_name_full = str_replace( str_replace( b_name, pattern = ".ident.tsv",""), "\\.", "_" )
    query_vec = rep(b_name_full, nrow(b_table))

    b_name_ccls = as.character(b_table$CCL)
    b_name_ccls_plane = str_to_upper( as.character(sapply( b_name_ccls, FUN = str_replace_all, pattern = "[\\(\\*\\-\\_\\+\\)]","")) )
    b_name_ccls_full = paste( as.character(b_name_ccls), b_table$Library, sep ="_" )

    ### new
    
    identifier_index = which( as.character(gold_t$CL_plane) == as.character(b_name_plane)  )
    
    identical_cls = as.character( gold_t$Name_identical[ identifier_index ] )
    identical_cls = as.character( unlist( str_split( identical_cls, "," ) ) )
    identical_cls = identical_cls[ identical_cls != ""]
    identical_cls_vec = b_name_ccls_full %in% identical_cls
    
    related_cls   = as.character( gold_t$Related[ identifier_index ] )
    related_cls   = as.character( unlist( str_split( related_cls, "," ) ) )
    related_cls   = related_cls[ related_cls != ""]
    related_cls_vec = b_name_ccls_full %in% related_cls
    
    identical_related_vcf = identical_cls_vec | related_cls_vec
    
    res_table <<- data.frame(
        "Query_CCL"        = c( as.character( res_table$Query_CCL ) , as.character( query_vec ) ),
        "Library_CCL"      = c( as.character( res_table$Library_CCL ) , as.character( b_name_ccls_full ) ),
        "CCL_to_be_found"      = c( as.character( res_table$CCL_to_be_found ) , as.character( identical_related_vcf ) ),
        "CCL_passed_threshold" = c( as.character( res_table$CCL_passed_threshold), as.character( b_table$Identification_sig ) ),
        "Found_muts_abs"   = c( as.character( res_table$Found_muts_abs), as.character( b_table$Matches ) ),
        "P_value"          = c( as.character( res_table$P_value), as.character( b_table$P_values ) )
    )
    res_table = res_table[ as.character( res_table$Found_muts_abs ) != "0",]
    
    # identification table
    
    if (exclude_self)
        b_table = b_table[ b_name_ccls_full != b_name_full, ]
    
    to_be_found = b_name_ccls_full[identical_related_vcf]
    found = b_name_ccls_full[ b_table$Identification_sig ]

    true_pos_ident  = to_be_found[ to_be_found %in% found ]
    false_neg_ident = to_be_found[!(to_be_found %in% found)]
    
    false_pos_ident = found[!(found %in% to_be_found)]
    
    to_be_found_cls = concat_me( to_be_found  )
    true_pos_ident_cls  = concat_me( true_pos_ident  )
    false_pos_ident_cls = concat_me( false_pos_ident )
    false_neg_ident_cls = concat_me( false_neg_ident )
    found_cls           = concat_me( found  )
    
    res_ident_table <<- data.frame( 
        "Query"                     = c( as.character( res_ident_table$Query ),                    as.character( b_name_full ) ),
        "Expected"                  = c( as.character( res_ident_table$Expected ),                 as.character( to_be_found_cls  ) ),
        "Found"                     = c( as.character( res_ident_table$Found ),                    as.character( found_cls ) ),
        "True_positive"             = c( as.character( res_ident_table$True_positive ),            as.character( true_pos_ident_cls  ) ),
        "False_negative"            = c( as.character( res_ident_table$False_negative ),           as.character( false_neg_ident_cls  ) ),
        "False_positive"            = c( as.character( res_ident_table$False_positive ),           as.character( false_pos_ident_cls   )),
        stringsAsFactors = F
=======
    
    for (b_file in b_files)
      parse_identification_data(b_file)
    
    #sapply( b_files, FUN = parse_identification_data)
    
    res_table = res_table[ order( as.double( as.character(res_table$Found_muts_abs) ), decreasing = T),  ]
    write.table( res_ident_table, benchmark_ident_file_path, sep ="\t", quote = F, row.names = F  )
    write.table( res_table, benchmark_res_file_path, sep ="\t", quote =F, row.names = F )
  }
  
parse_identification_data = function( b_file ){
    
    print(c(b_file, as.character(which(b_files %in% b_file))))
    
    b_table = read.table(b_file, sep ="\t", header = T, stringsAsFactors = FALSE)
    library_names = Uniquorn::read_library_names(ref_gen = "GRCH37")
    
    for (library_name in library_names){
        if ( str_detect(pattern = library_name, b_file) )
            source_file <<- library_name
    }
    
    b_file
    
    b_name = tail(unlist(str_split(b_file, pattern = "/")  ),1)
    b_name_plane = str_replace( b_name, pattern = paste( c(".",source_file,".ident.tsv"), sep ="", collapse = ""),"")
    b_name_plane = str_replace_all( b_name_plane, pattern = ,"[\\(\\*\\-\\-\\_\\+\\-\\)-]", "")
    b_name_full = str_replace( str_replace( b_name, pattern = ".ident.tsv",""), "\\.", "_" )
    query_vec = rep(b_name_full, nrow(b_table))
    
    b_name_ccls = str_to_upper( as.character(b_table$CCL) )
    b_name_ccls_plane = str_to_upper( as.character(sapply( b_name_ccls, FUN = str_replace_all, pattern = "[\\(\\*\\-\\_\\+\\)]","")) )
    b_name_ccls_full = paste( as.character(b_name_ccls_plane), b_table$Library, sep ="_" )
    
    ### new
    
    gold_cls = str_to_upper(as.character(gold_t$CL_plane))
    ident_cls = str_to_upper(as.character(b_name_plane))
    ident_cls = str_replace_all(ident_cls, pattern = "\\.","")
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
    ident_ident = str_to_upper( as.character(sapply( ident_ident, FUN = str_replace_all, pattern = "[\\(\\*\\-\\_\\+\\)]","")) )
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
    related_ident = str_to_upper( as.character(sapply( related_ident, FUN = str_replace_all, pattern = "[\\(\\*\\-\\_\\+\\)]","")) )
    related_ident = paste(related_ident, related_lib, sep = "_")
    
    related_cls_vec = b_name_ccls_full %in% related_ident
    
    identical_related_vcf = identical_cls_vec | related_cls_vec
    
    res_table <<- data.frame(
      "Query_CCL"        = c( as.character( res_table$Query_CCL ) , as.character( query_vec ) ),
      "Library_CCL"      = c( as.character( res_table$Library_CCL ) , as.character( b_name_ccls_full ) ),
      "CCL_to_be_found"      = c( as.character( res_table$CCL_to_be_found ) , as.character( identical_related_vcf ) ),
      "CCL_passed_threshold" = c( as.character( res_table$CCL_passed_threshold), as.character( b_table$Identification_sig ) ),
      "Found_muts_abs"   = c( as.character( res_table$Found_muts_abs), as.character( b_table$Matches ) ),
      "P_value"          = c( as.character( res_table$P_value), as.character( b_table$P_values ) )
    )
    res_table = res_table[ as.character( res_table$Found_muts_abs ) != "0",]
    
    # identification table
    
    to_be_found = b_name_ccls_full[identical_related_vcf]
    found = str_to_upper( as.character( b_name_ccls_full[ b_table$Identification_sig ] ) )
    split_found = lapply( found, FUN = function(vec){return( unlist(str_split(vec, pattern = "_")  ) )} )
    found_len = sapply(split_found, FUN = length)
    found_lib = as.character( sapply( split_found, FUN =  tail, 1 ) )
    found_ident <<- c()
    
    if (length(found) > 0)
    for(i in 1:length(found)){
        found_ident <<- c( found_ident, head( as.character( unlist(split_found[i])), found_len[i] - 1) )
    }
    found_ident = str_to_upper( as.character(sapply( found_ident, FUN = str_replace_all, pattern = "[\\(\\*\\-\\_\\+\\)]","")) )
    found_ident = paste(found_ident, found_lib, sep = "_")
    
    true_pos_ident  = to_be_found[ to_be_found %in% found ]
    false_neg_ident = to_be_found[!(to_be_found %in% found)]
    
    false_pos_ident = found[!(found %in% to_be_found)]
    
    to_be_found_cls = concat_me( to_be_found  )
    true_pos_ident_cls  = concat_me( true_pos_ident  )
    false_pos_ident_cls = concat_me( false_pos_ident )
    false_neg_ident_cls = concat_me( false_neg_ident )
    found_cls           = concat_me( found  )
    
    res_ident_table <<- data.frame( 
      "Query"                     = c( as.character( res_ident_table$Query ),                    as.character( b_name_full ) ),
      "Expected"                  = c( as.character( res_ident_table$Expected ),                 as.character( to_be_found_cls  ) ),
      "Found"                     = c( as.character( res_ident_table$Found ),                    as.character( found_cls ) ),
      "True_positive"             = c( as.character( res_ident_table$True_positive ),            as.character( true_pos_ident_cls  ) ),
      "False_negative"            = c( as.character( res_ident_table$False_negative ),           as.character( false_neg_ident_cls  ) ),
      "False_positive"            = c( as.character( res_ident_table$False_positive ),           as.character( false_pos_ident_cls   )),
      stringsAsFactors = F
>>>>>>> 52e010a4c1d3331df1913c4c690c119fc97a08c9
    )
    
    #auc_table$False_neg[ index ] = as.character( as.integer( auc_table$False_neg[ index ] ) + length( tmp_false_neg ) )
    write.table( res_ident_table, benchmark_ident_file_path, sep ="\t", quote = F, row.names = F  )
    write.table( res_table, benchmark_res_file_path, sep ="\t", quote =F, row.names = F )
<<<<<<< HEAD

=======
  
>>>>>>> 52e010a4c1d3331df1913c4c690c119fc97a08c9
}

#inclusion_weight = 0
run_benchmark(
    inclusion_weight,
<<<<<<< HEAD
    only_first,
    exclude_self,
    run_identification = run_identification,
    #confidence_score = confidence_score,
    minimum_matching_mutations = minimum_matching_mutations,
    type_benchmark = type_benchmark,
    ref_gen = ref_gen
=======
    minimum_matching_mutations = minimum_matching_mutations
>>>>>>> 52e010a4c1d3331df1913c4c690c119fc97a08c9
)
