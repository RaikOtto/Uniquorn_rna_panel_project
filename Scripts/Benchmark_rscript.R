library("stringr")
library("BiocParallel")
library("doParallel")
library("foreach")
library("devtools")
library("argparse")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

parser = ArgumentParser()
parser$add_argument('-iw', "--inclusion_weight", type="double")
parser$add_argument('-p', "--panel", action="store_true", default = "FALSE")
parser$add_argument('-nt', "--number_threads", type="integer", default = "1")
args = parser$parse_args()

only_first           = FALSE
exclude_self         = FALSE
p_value = .05

minimum_matching_mutations = args$inclusion_weight
panel_mode                 = args$panel_mode

run_identification   = F
auc_mode             = FALSE
ref_gen              = "GRCH37"
#type_benchmark       = "non_regularized"
type_benchmark       = "regularized"

run_identification = function(
    inclusion_weight,
    minimum_matching_mutations,
    type_benchmark = type_benchmark,
    ref_gen,
    panel_mode,
    number_threads
){
  
    source("~/Uniquorn_rna_panel_project/Scripts/utility.R")
    
    build_path_variables( 
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
    
    ### !!! ###
    
    raw_files_path = "~/Uniquorn_data/benchmark_vcf_files/raw_files/"
    
    build_path_variables( 
        inclusion_weight = inclusion_weight,
        only_first = only_first,
        exclude_self = exclude_self,
        run_identification = run_identification,
        cellminer = F,
        type_benchmark = type_benchmark
    )
    
    ### PANEL ### !!
    
    if (panel_mode){
        
        raw_files_path = "~/Uniquorn_data/benchmark_vcf_files/panel_raw_files/"
        out_path_ident = paste(c("~/Uniquorn_data/benchmark_vcf_files/panel_ident_files_",type_benchmark,
                               "/",inclusion_weight,"/"),sep="",collapse= "")
    } else {
        out_path_ident = paste(c("~/Uniquorn_data/benchmark_vcf_files/ident_files_",type_benchmark,
                               "/",inclusion_weight,"/"),sep="",collapse= "")
    }
    i_files = list.files( raw_files_path, pattern = ".vcf$", full.names = T )
    
    
    
    if( number_threads > 1 ){
      
        doParallel::registerDoParallel(number_threads)  
        foreach::foreach(
            i_file = i_files
        ) %dopar% {
            print(c(match(i_file , i_files),i_file))
            file_name = tail(as.character(unlist(str_split(i_file,pattern = "/"))),1)
            file_name = str_replace(file_name,pattern =".vcf",".ident.tsv")
            out_path_ident_file = paste(out_path_ident, file_name,sep ="")
            
            if( !file.exists(out_path_ident_file) ){
              
                file.create(out_path_ident_file)
                identify_vcf_file(
                    i_file,
                    mutational_weight_inclusion_threshold = inclusion_weight,
                    output_file = out_path_ident_file
                )
            }
        }
        
        doParallel::stopImplicitCluster()
        
    } else {
    
        for (i_file in i_files){
          
            print(c(match(i_file , i_files),i_file))
            file_name = tail(as.character(unlist(str_split(i_file,pattern = "/"))),1)
            file_name = str_replace(file_name,pattern =".vcf",".ident.tsv")
            out_path_ident_file = paste(out_path_ident, file_name,sep ="")
            
            if( !file.exists(out_path_ident_file) ){
                
                file.create(out_path_ident_file)
                identify_vcf_file(
                    i_file,
                    mutational_weight_inclusion_threshold = inclusion_weight,
                    output_file = out_path_ident_file
                )
            }
        }
    }
}

run_identification(
    inclusion_weight = args$inclusion_weight,
    minimum_matching_mutations = minimum_matching_mutations,
    type_benchmark = type_benchmark,
    ref_gen = ref_gen,
    panel_mode = args$panel,
    number_threads = args$number_threads
)
