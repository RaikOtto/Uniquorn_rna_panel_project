library("stringr")
library("Uniquorn")

inclusion_weight      =1.0
only_first            = FALSE
exclude_self          = FALSE
run_identification    = F
cellminer_v2          = F
distinct_mode         = TRUE
identical_mode        = FALSE
distuinguished_panels = T
panel_mode            = FALSE
#type_benchmark        = 'non_regularized'
type_benchmark        = 'regularized'

run_small_statistics = function( 
  input_path_comparison = input_path_comparison,
  rigid_similarity = rigid_similarity, 
  inclusion_weight = inclusion_weight, 
  only_first = only_first, 
  exclude_self = exclude_self,
  run_identification = run_identification, 
  panel_mode,
  distuinguished_panels = distuinguished_panels,
  type_benchmark = type_benchmark
){
  
    source("~//Uniquorn_rna_panel_project//Scripts/utility.R")
    source("~/Uniquorn_rna_panel_project//Scripts/Benchmark_Merger_Skript_Relaxed.R")
    build_path_variables( 
        inclusion_weight = inclusion_weight, 
        only_first = only_first, 
        exclude_self = exclude_self, 
        run_identification = run_identification, 
        cellminer = F,
        type_benchmark = type_benchmark
    )
    
    print( c("Inclusion weight", inclusion_weight ) )
  
    if (type_benchmark == "regularized"){
        
        input_path_comparison = paste(c(
            "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_regularized/",
            inclusion_weight
            ,"_Benchmark_comparisons_result.tab"), sep ="", collapse = "")
        input_path_identification = paste( c(
          "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_regularized/",
          inclusion_weight,
          "_Benchmark_identification_result.tab"), sep ="", collapse = ""
        )
    }
    
    if (type_benchmark == "non_regularized"){
      
        input_path_comparison = paste(c(
          "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_non_regularized/",
          inclusion_weight
          ,"_Benchmark_comparisons_result.tab"), sep ="", collapse = "")
        input_path_identification = paste( c(
          "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_non_regularized/",
          inclusion_weight,
          "_Benchmark_identification_result.tab"), sep ="", collapse = ""
        )
    }
    
    output_table_path = str_replace( 
      input_path_comparison, 
      pattern = "_Benchmark_comparisons_result.tab",
      "_Benchmark_identification_result_aggregated.tab"
    )
    
    if (panel_mode){
        ident_result_files_path = str_replace(ident_result_files_path,pattern = "ident_files","panel_ident_files")
        benchmark_ident_file_path = str_replace(benchmark_ident_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
        benchmark_res_file_path = str_replace(benchmark_res_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
        output_table_path = str_replace(output_table_path,pattern = "Benchmark_results","panel_Benchmark_results")
        input_path_identification = str_replace(input_path_identification,pattern = "Benchmark_results","panel_Benchmark_results")
        
    }
    
    #run_relaxed_merging( 
    #    input_path_identification,
    #    output_table_path 
    #)

    t_rel = read.table( input_path_identification, sep ="\t", header = T)
  
    expected    = paste( c(as.character( t_rel$Expected[ t_rel$Expected!= ""] )), collapse = ", ",sep = "" )
    nr_expected = length((as.character(unlist(str_split(expected,", ")))))
    print(paste(c("Expected:",nr_expected),collapse = "",sep = "") )
    
    true_positives = as.character(unlist(str_split(t_rel$True_positive,", ")))
    true_positives = true_positives[ true_positives != ""  ]
    true_positives = true_positives[ !is.na(true_positives) ]
    nr_true_positives = length(true_positives)
    print(paste(c("TP:",as.character(nr_true_positives)),collapse = "",sep = ""))
    
    false_negatives = as.character(unlist(str_split(t_rel$False_negative,", ")))
    false_negatives = false_negatives[ false_negatives!= ""]
    false_negatives = false_negatives[ !is.na(false_negatives) ]
    nr_false_negatives = length( false_negatives )
    print(paste(c("FN: ",as.character(nr_false_negatives)),collapse = "",sep = ""))
    
    false_positives = t_rel$False_positive[t_rel$False_positive!= ""]
    nr_false_positive = length((as.character(unlist(str_split(false_positives,",")))))
    print(paste(c("FP: ",nr_false_positive),collapse = "",sep = ""))
    
    nr_true_negatives  = nrow(t_rel)**2 - nrow(t_rel) - nr_false_negatives
    
    TPR = round(nr_true_positives / (nr_true_positives + nr_false_positive),3) * 100
    print(paste(c("TPR: ",TPR),collapse = "",sep = ""))
    
    TNR = round(nr_true_negatives  / (nr_true_negatives + nr_false_negatives),10) * 100
    print(paste(c("TNR: ",TNR),collapse = "",sep = ""))
    
    PPV = round( nr_true_positives / ( nr_true_positives + nr_false_negatives ),3) * 100
    print(paste(c("PPV: ",PPV),collapse = "",sep = ""))
    
    F1 = round( (2* nr_true_positives) / ((2*nr_true_positives) + nr_false_positive + nr_false_negatives) ,3) * 100
    print(paste(c("F1: ",F1),collapse = "",sep = ""))
}
run_small_statistics(
    inclusion_weight = inclusion_weight,
    only_first = only_first,
    exclude_self = exclude_self,
    run_identification = run_identification,
    panel_mode = panel_mode,
    type_benchmark = type_benchmark
)
