library("stringr")
source("~/Uniquorn_rna_panel_project//Scripts/utility.R")
library("readODS")
library("devtools")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()


### load ods benchmark

avg_file   = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 1)
rna_file   = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 2)
TruSight   = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 3)
ClearSight = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 4)
Hotspot    = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 5)

seq_type = factor( 
  c(rep("AVG",16),rep("RNA",16),rep("TruSight",16),rep("ClearSight",16),rep("Hotspot_v2",16)),
  levels = c("AVG","RNA","ClearSight","TruSight","Hotspot_v2"))
Weights = factor(
  as.character( rep( c( rep(1.0,4),rep(0.5,4),rep(0.25,4),rep(0.0,4) ), 5 ) ),
  levels = c("1","0.5","0.25","0")
)
TPs = as.double( c( avg_file[3,-1], rna_file[3,-1], TruSight[3,-1],ClearSight[3,-1],Hotspot[3,-1] ) )

# Sensitivity

Sensitivity = as.double( c( avg_file[6,-1], rna_file[6,-1], TruSight[6,-1],ClearSight[6,-1],Hotspot[6,-1] ) )
seq_type = factor( 
  c(rep("AVG",16),rep("RNA",16),rep("TruSight",16),rep("ClearSight",16),rep("Hotspot_v2",16)),
  levels = c("AVG","RNA","ClearSight","TruSight","Hotspot_v2"))


benchmark_Sensitivity = data.frame(
  "Weight" = Weights,
  "Sensitivity" = Sensitivity,
  "Seq_Type" = seq_type
)

# F1

F1s = as.double( c( avg_file[8,-1], rna_file[8,-1], TruSight[8,-1],ClearSight[8,-1],Hotspot[8,-1] ) )

benchmark_F1 = data.frame(
  "Weight" = Weights,
  "F1" = F1s,
  "Seq_Type" = seq_type
)

# PPV

PPVs = as.double( c( avg_file[7,-1], rna_file[7,-1], TruSight[7,-1],ClearSight[7,-1],Hotspot[7,-1] ) )

benchmark_PPV = data.frame(
  "Weight" = Weights,
  "PPV" = PPVs,
  "Seq_Type" = seq_type
)

### abs_mat

library_names = Uniquorn::read_library_names(ref_gen = "GRCH37")
mat_mat = show_contained_ccls()

abs_mat <<- data.frame(
  "Library" = as.character(),
  "Weight" = as.character(),
  "Count" = as.integer(),
  stringsAsFactors = F
)

for ( var_weight in c("W0","W25","W05","W1")){
  for (library_name in library_names){
    
    lib_mat = mat_mat[ mat_mat$Library == library_name,]
    data_lib = as.integer( lib_mat[[var_weight]] )

    new_mat = data.frame(
      "Library" = rep(library_name, length(data_lib)),
      "Weight"  = rep(var_weight,   length(data_lib)),
      "Count"   = data_lib
    )
    abs_mat   = rbind(abs_mat, new_mat)
  }
}
abs_mat$Library = as.character(abs_mat$Library)
abs_mat$Library[abs_mat$Library=="COSMIC"] = "CPG"
abs_mat$Library[abs_mat$Library=="EGA"] = "Klijn"
abs_mat$Library = factor(abs_mat$Library, levels = c("CCLE","CELLMINER","CPG","Klijn","GDC"))

### mean mat

mean_mat <<- data.frame(
  "Library" = as.character(),
  "Weight" = as.character(),
  "Count" = as.double(),
  stringsAsFactors = F
)
libraries = as.character( unique(as.character(abs_mat$Library ) ))

for ( var_weight in c("W0","W25","W05","W1")){
  for (library_name in libraries){
    
    lib_vec = abs_mat[ abs_mat$Library == library_name,]
    lib_vec = lib_vec[ lib_vec$Weight == var_weight, ]
    
    mean_mat <<- data.frame(
      "Library" = c( as.character( mean_mat$Library ), library_name ),
      "Weight"  = c( as.character( mean_mat$Weight  ), var_weight   ),
      "Count"   = as.double( c( as.double(mean_mat$Count), mean( as.double( lib_vec$Count)) ) )
    )
  }
}

