library("stringr")
library("devtools")
library("argparse")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

parser = ArgumentParser()
parser$add_argument('-iw', "--inclusion_weight", type="double")
parser$add_argument('-p', "--panel_mode", action="store_true", default = "FALSE")
parser$add_argument('-r', "--robust_mode", action="store_true", default = "FALSE")
parser$add_argument('-cs', "--confidence_score", type="integer", default = "40")
args = parser$parse_args()

only_first           = FALSE
exclude_self         = FALSE

minimum_matching_mutations = args$inclusion_weight
panel_mode                 = args$panel_mode
robust_mode                = args$robust_mode
confidence_score           = args$confidence_score

run_identification   = F
auc_mode             = FALSE
ref_gen              = "GRCH37"
type_benchmark       = "regularized"

run_identification = function(
    inclusion_weight,
    minimum_matching_mutations,
    type_benchmark = type_benchmark,
    ref_gen,
    panel_mode,
    robust_mode,
    p_value
){
  
    source("~/Uniquorn_rna_panel_project/Scripts/utility.R")
  
    inclusion_weight_path = as.character( inclusion_weight )
    inclusion_weight_path = str_replace( as.character( inclusion_weight_path ), pattern ="\\.", "_" )
    
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
    
    #seen_obj_path = str_replace(out_path_ident, pattern = ".Benchmark_identification_result.tab", ".seen_obj.RDS")
    
    ### !!! ###
    
    raw_files_path = "~/Uniquorn_data/benchmark_vcf_files/raw_files/"
    
    inclusion_weight = str_replace(inclusion_weight, pattern = "\\.", "_")
    
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
    i_files = list.files( raw_files_path, pattern = ".vcf$", full.names = T ,ignore.case = TRUE)
    
    #if (file.exists(seen_obj_path)) {seen_obj <<- readRDS(seen_obj_path)
    #} else {seen_obj <<- c()}
    
    
    for (i_file in i_files){
      
        print(c(match(i_file , i_files),i_file))
        file_name = tail(as.character(unlist(str_split(i_file,pattern = "/"))),1)
        file_name = str_replace(file_name,pattern =".(vcf)|(.VCF)",".ident.tsv")
        out_path_ident_file = paste(out_path_ident, file_name,sep ="")
        
        #if (file_name %in% seen_obj){next()
        #} else {
        #    seen_obj <<- c(seen_obj, file_name)
        #    saveRDS(seen_obj, file = seen_obj_path)
        #}
        
        if( !file.exists(out_path_ident_file) ){
            
            file.create(out_path_ident_file)
            identify_vcf_file(
                i_file,
                mutational_weight_inclusion_threshold = inclusion_weight,
                output_file = out_path_ident_file,
                robust_mode = robust_mode,
                p_value = p_value
            )
        }
    }
}

run_identification(
    inclusion_weight = args$inclusion_weight,
    minimum_matching_mutations = minimum_matching_mutations,
    type_benchmark = type_benchmark,
    ref_gen = ref_gen,
    panel_mode = panel_mode,
    robust_mode = robust_mode,
    confidence_score = confidence_score
)
