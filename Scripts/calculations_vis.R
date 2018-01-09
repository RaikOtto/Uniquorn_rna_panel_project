library("ggplot2")
library("stringr")
library("devtools")
library("argparse")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

library_names = Uniquorn::read_library_names(ref_gen = "GRCH37")
mat_mat = show_contained_cls()

variant_t <<- data.frame(
    "Library" = as.character(),
    "Weight" = as.character(),
    "Count" = as.integer()
)

for ( var_weight in c("W0","W25","W05","W1")){
    for (library_name in library_names){
      
        lib_mat = mat_mat[ mat_mat$Library == library_name,]
        data_lib = as.integer( lib_mat[[var_weight]] )
        
        new_mat = data.frame(
            "Library" = rep(library_name, length(data_lib)),
            "Weight" = rep(var_weight, length(data_lib)),
            "Count" = data_lib
        )
        variant_t = rbind(variant_t, new_mat)
    }
}


