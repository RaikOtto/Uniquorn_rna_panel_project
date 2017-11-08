library("stringr")
library("devtools")
library("argparse")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

parser = ArgumentParser()
parser$add_argument('-iw', "--inclusion_weight", type="double")
parser$add_argument('-p', "--panel_mode", action="store_true", default = "FALSE")
parser$add_argument('-nt', "--number_threads", type="integer", default = "1")
args = parser$parse_args()

#args$inclusion_weight = .5
inclusion_weight = args$inclusion_weight
panel_mode = args$panel_mode
#inclusion_weight = 0.25
#panel_mode = FALSE

run_small_statistics = function( 
  input_path_comparison = input_path_comparison,
  inclusion_weight = inclusion_weight, 
  panel_mode
){
    inclusion_weight = str_replace(inclusion_weight, pattern = "\\.", "_")
    source("~//Uniquorn_rna_panel_project//Scripts/utility.R")
    
    build_path_variables( 
        inclusion_weight = inclusion_weight, 
        only_first = FALSE, 
        exclude_self = FALSE, 
        run_identification = FALSE, 
        cellminer = FALSE,
        type_benchmark = "regularized"
    )
    
    print( paste( c("Inclusion weight: ", inclusion_weight ), sep = "", collapse = "" ) )
    
    input_path_comparison = paste(c(
      "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_regularized/",
      inclusion_weight
      ,"_Benchmark_comparisons_result.tab"), sep ="", collapse = "")
    input_path_identification = paste( c(
      "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_regularized/",
      inclusion_weight,
      "_Benchmark_identification_result.tab"), sep ="", collapse = ""
    )
    
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
    
    t_rel = read.table( input_path_identification, sep ="\t", header = T)
    print(paste(c("Number cases: ", nrow(t_rel) ), sep ="", collapse=""))
    
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
    
    Sensitivity = round(nr_true_positives / (nr_true_positives + nr_false_negatives),3) * 100
    print(paste(c("Sensitivity: ",Sensitivity),collapse = "",sep = ""))
    
    PPV = round( nr_true_positives / ( nr_true_positives + nr_false_positive ),3) * 100
    print(paste(c("PPV: ",PPV),collapse = "",sep = ""))
    
    F1 = round( (2* nr_true_positives) / ((2*nr_true_positives) + nr_false_positive + nr_false_negatives) ,3) * 100
    print(paste(c("F1: ",F1),collapse = "",sep = ""))
}
#panel_mode = FALSE
#inclusion_weight = 1.0
run_small_statistics(
    inclusion_weight = inclusion_weight,
    panel_mode = panel_mode
)
