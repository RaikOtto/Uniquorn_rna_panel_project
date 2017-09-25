library("stringr")
source("~/Dropbox/PhD/Uniquorn_project/Benchmark.R")
library("Uniquorn")

c = show_contained_cls()

cellminer_cls = parse_string( c$CL[grepl("CELLMINER", c$CL)], "CELLMINER" )
ccle_cls = parse_string( c$CL[grepl("CCLE", c$CL)], "CCLE" )
cosmic_cls = parse_string( c$CL[grepl("COSMIC", c$CL)], "COSMIC" )

length( c( cellminer_cls, ccle_cls, cosmic_cls ) ) # 1989 CL images
length( unique( c( cellminer_cls, ccle_cls, cosmic_cls ) ) ) # 1342 genuine CL images

intersect_cosmic_ccle = cosmic_cls[ unique(cosmic_cls) %in% unique(ccle_cls) ]
intersect_cosmic_cellminer = cosmic_cls[ unique(cosmic_cls) %in% unique(cellminer_cls) ]
intersect_ccle_cellminer = ccle_cls[ unique(ccle_cls) %in% unique(cellminer_cls) ]

# real world benchmark integration

benchmark_res_t = read.table("/Users/raikotto/Dropbox/PhD/Uniquorn_project/benchmark_vcf_files//Benchmark_result_1_not_only_first_exclude_self_distinguished_panels_results.tab",sep ="\t", header =T)
#benchmark_res_t = read.table("/Users/raikotto/Dropbox/PhD/Uniquorn_project/benchmark_vcf_files//Benchmark_result_1_not_only_first_not_exclude_self_distinguished_panels_results.tab",sep ="\t", header =T)

result_table_intersect = benchmark_res_t[ 0, ]

# cosmic ccle
for (cl in c( intersect_cosmic_ccle )){
  
  benchmark_res_t_sub  = benchmark_res_t[ benchmark_res_t$Source_query == "COSMIC", ]
  selection = benchmark_res_t_sub[ benchmark_res_t_sub$CL_name_query == cl  ,]
  selection = selection[ selection$CL_name_training == cl  ,  ]
  
  result_table_intersect = rbind( result_table_intersect, selection  )
  
}

# cosmic cellminer
for (cl in c( intersect_cosmic_cellminer )){
  
  benchmark_res_t_sub  = benchmark_res_t[ 
    (benchmark_res_t$Source_query %in% c("COSMIC","CELLMINER")) &
      (benchmark_res_t$Source_training %in% c("COSMIC","CELLMINER"))
    ,]
  selection = benchmark_res_t_sub[ benchmark_res_t_sub$CL_name_query == cl  ,]
  selection = selection[ selection$CL_name_training == cl  ,  ]
  
  result_table_intersect = rbind( result_table_intersect, selection  )
  
}

# ccle cellminer
for (cl in c( intersect_ccle_cellminer )){
  
  benchmark_res_t_sub  = benchmark_res_t[ 
    (benchmark_res_t$Source_query %in% c("CCLE","CELLMINER")) &
      (benchmark_res_t$Source_training %in% c("CCLE","CELLMINER"))
    ,]
  selection = benchmark_res_t_sub[ benchmark_res_t_sub$CL_name_query == cl  ,]
  selection = selection[ selection$CL_name_training == cl  ,  ]
  
  result_table_intersect = rbind( result_table_intersect, selection  )
  
}

#write.table(result_table_intersect, "~/Dropbox/PhD/Uniquorn_project/Pub/Results/Benchmark/intersect_name_identical.tab",sep ="\t", quote = F, row.names = F)
write.table(result_table_intersect, "~/Dropbox/PhD/Uniquorn_project/Pub/Results/Benchmark/intersect_name_identical_count_self.tab",sep ="\t", quote = F, row.names = F)

### part not identified CLs

library("stringr")
source("~/Dropbox/PhD/Uniquorn_project/Benchmark.R")
library("Uniquorn")

c = show_contained_cls()

cellminer_cls = parse_string( c$CL[grepl("CELLMINER", c$CL)], "CELLMINER" )
ccle_cls = parse_string( c$CL[grepl("CCLE", c$CL)], "CCLE" )
cosmic_cls = parse_string( c$CL[grepl("COSMIC", c$CL)], "COSMIC" )

all_contained = unique(c(cellminer_cls,ccle_cls,cosmic_cls))
length(all_contained) # 1342

intersect_cosmic_ccle = cosmic_cls[ unique(cosmic_cls) %in% unique(ccle_cls) ]
intersect_cosmic_cellminer = cosmic_cls[ unique(cosmic_cls) %in% unique(cellminer_cls) ]
intersect_ccle_cellminer = ccle_cls[ unique(ccle_cls) %in% unique(cellminer_cls) ]

cls_in_cosmic_to_ccle_cellminer = intersect_cosmic_ccle[c(intersect_cosmic_ccle %in% intersect_cosmic_cellminer)]
contained_three = cls_in_cosmic_to_ccle_cellminer[ cls_in_cosmic_to_ccle_cellminer  %in% intersect_ccle_cellminer ]

length( contained_three ) # 34 in all three panels => 34*6 comparisons = 204

all_intersects = unique( c( intersect_cosmic_ccle, intersect_cosmic_cellminer, intersect_ccle_cellminer ) )
length(all_intersects) # 613

all_intersects[ !(all_intersects %in% contained_three)]

nmbr_identifications = (length( intersect_cosmic_ccle ) + length( intersect_cosmic_cellminer ) + length( intersect_ccle_cellminer )) * 2

length( unique(c(cellminer_cls,cosmic_cls,ccle_cls)) )
length( unique(c(intersect_cosmic_ccle, intersect_cosmic_cellminer,intersect_ccle_cellminer)) )

cellminer_only = cellminer_cls[!( cellminer_cls %in% unique(c(intersect_cosmic_cellminer, intersect_ccle_cellminer)))]
cosmic_only = cosmic_cls[!( cosmic_cls %in% unique(c(intersect_cosmic_cellminer, intersect_cosmic_ccle)))]
ccle_only = ccle_cls[!( ccle_cls %in% unique(c(intersect_ccle_cellminer, intersect_cosmic_ccle)))]

length(cellminer_only)
length(cosmic_only)
length(ccle_only)

name_not_twice_contained = c( unique(cosmic_cls),unique(ccle_cls), unique(cellminer_cls)  )
name_not_twice_contained = name_not_twice_contained[ which( !( name_not_twice_contained %in%  intersect_cosmic_ccle)      ) ]
name_not_twice_contained = name_not_twice_contained[ which( !( name_not_twice_contained %in%  intersect_cosmic_cellminer) ) ]
name_not_twice_contained = name_not_twice_contained[ which( !( name_not_twice_contained %in%  intersect_ccle_cellminer  ) ) ]

length(name_not_twice_contained)

same_name_t = read.table("~/Dropbox/PhD/Uniquorn_project/Pub/Results/Benchmark/intersect_name_identical.tab",sep ="\t", header =T)

candidates = seq(dim(benchmark_res_t)[1])

for (i in seq(dim(benchmark_res_t)[1])){
  
  cl1 = as.character( benchmark_res_t[i,1] )
  cl2 = as.character( benchmark_res_t[i,2] )
  if(cl1 == cl2)
    candidates[i] = 0
}

