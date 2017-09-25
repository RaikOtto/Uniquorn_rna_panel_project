
### venn diagrams

library("Uniquorn")
library("stringr")
library("limma")
library("VennDiagram")

### gold

gold_t = read.table("~/Dropbox/PhD/Uniquorn_project/Pub/Goldstandard.tab",sep ="\t", header = T, stringsAsFactors = F)

# rigid
all_cls = read.table("~/Dropbox/PhD/Uniquorn_project/benchmark_vcf_files/Benchmark_identification_result_1_not_only_first_not_exclude_self_distinguished_panels__Rigid_results.tab",sep ="\t", header = T, stringsAsFactors = F)

cls = paste( all_cls$Name_query, all_cls$Source_query, sep ="_" )

# relaxed

all_cls = read.table("~/Dropbox/PhD/Uniquorn_project/benchmark_vcf_files/Benchmark_identification_result_1_not_only_first_not_exclude_self_distinguished_panels__Relaxed_results.tab",sep ="\t", header = T, stringsAsFactors = F)

cls = paste( all_cls$Name_query, all_cls$Source_query, sep ="_" )

# both

member_cosmic = grepl("_COSMIC",    cls)
member_ccle   = grepl("_CCLE",      cls)
member_cell   = grepl("_CELLMINER", cls)

expected_list = sapply( all_cls$Expected, FUN = str_split, ", " )

intersect_cosmic = grepl( "_COSMIC",   all_cls$Expected   )
intersect_ccle   = grepl( "_CCLE",      all_cls$Expected )
intersect_cell   = grepl( "_CELLMINER", all_cls$Expected )

intersect_cosmic_ccle    =     intersect_cosmic   &     intersect_ccle   & (! intersect_cell )
intersect_cosmic_cell    =     intersect_cosmic   & ( ! intersect_ccle ) &    intersect_cell
intersect_ccle_cell      = ( ! intersect_cosmic ) &     intersect_ccle   &    intersect_cell

cosmic_exclusive         =     intersect_cosmic   & ( ! intersect_ccle ) & (! intersect_cell )
ccle_exclusive           = ( ! intersect_cosmic ) &     intersect_ccle   & (! intersect_cell )
cell_exclusive           = ( ! intersect_cosmic ) &  (! intersect_ccle ) &    intersect_cell

intersect_all            =     intersect_cosmic   &     intersect_ccle   &    intersect_cell

nr_intersect_cosmic_ccle =     length( unique( str_split( all_cls$Expected[ intersect_cosmic_ccle ], ", " ) ) ) #* 2
nr_intersect_cosmic_cell =     length( unique( str_split( all_cls$Expected[ intersect_cosmic_cell ], ", " ) ) ) #* 2
nr_intersect_ccle_cell   =     length( unique( str_split( all_cls$Expected[ intersect_ccle_cell   ], ", " ) ) ) #* 2

nr_intersect_all         =     length( unique( str_split( all_cls$Expected[ intersect_all ], ", " ) ) ) #* 6

nr_cosmic_exclusive      = sum( cosmic_exclusive )
nr_ccle_exclusive        = sum( ccle_exclusive   )
nr_cell_exclusive        = sum( cell_exclusive   )

nr_intersect_cosmic_ccle
nr_intersect_cosmic_cell
nr_intersect_ccle_cell

nr_intersect_all

sum_nr_self_intersect_cosmic = sum( grepl( all_cls$Expected, pattern = ".*(_COSMIC)+.*(_COSMIC)+" ) )
sum_nr_self_intersect_ccle   = sum( grepl( all_cls$Expected, pattern = ".*(_CCLE)+.*(_CCLE)+" ) )
sum_nr_self_intersect_cell   = sum( grepl( all_cls$Expected, pattern = ".*(_CELLMINER)+.*(_CELLMINER)+" ) )

# venn

member_mat = matrix( integer(), nrow = 0, ncol = 3 )

member_mat = rbind( member_mat, cbind( rep( 1, nr_intersect_cosmic_ccle*2 ), rep( 1, nr_intersect_cosmic_ccle*2 ), rep( 0, nr_intersect_cosmic_ccle*2 ) ) )
member_mat = rbind( member_mat, cbind( rep( 1, nr_intersect_cosmic_cell*2 ), rep( 0, nr_intersect_cosmic_cell*2 ), rep( 1, nr_intersect_cosmic_cell*2 ) ) )

member_mat = rbind( member_mat, cbind( rep( 0, nr_intersect_ccle_cell*2   ), rep( 1, nr_intersect_ccle_cell*2   ), rep( 1, nr_intersect_ccle_cell*2   ) ) )

member_mat = rbind( member_mat, cbind( rep( 1, nr_intersect_all*6         ), rep( 1, nr_intersect_all*6         ), rep( 1, nr_intersect_all*6         ) ) )

cosmic_reps         = nr_cosmic_exclusive + ( sum( member_cosmic ) - nr_cosmic_exclusive ) + sum_nr_self_intersect_cosmic
ccle_reps           = nr_ccle_exclusive   + ( sum( member_ccle   ) - nr_ccle_exclusive   ) + sum_nr_self_intersect_ccle
cell_reps           = nr_cell_exclusive   + ( sum( member_cell   ) - nr_cell_exclusive   ) + sum_nr_self_intersect_cell

cosmic_substitution = cbind( rep( 1, cosmic_reps ), rep( 0, cosmic_reps ), rep( 0, cosmic_reps ) )
ccle_substitution   = cbind( rep( 0, ccle_reps   ), rep( 1, ccle_reps   ), rep( 0, ccle_reps   ) )
cell_substitution   = cbind( rep( 0, cell_reps   ), rep( 0, cell_reps   ), rep( 1, cell_reps   ) )

member_mat          = rbind( member_mat,cosmic_substitution,ccle_substitution,cell_substitution )
member_mat          = data.frame( member_mat )

colnames(member_mat) = c("COSMIC","CCLE","CELLMINER")

vennDiagram( member_mat )

venn.diagram(
  x = list(
    COSMIC = which( member_mat$COSMIC == 1 ),
    CCLE = which( member_mat$CCLE == 1 ),
    CELLMINER = which( member_mat$CELLMINER == 1 )
  ), 
  height=2000, width=2000, resolution=300,
  col = "transparent",
  margin = 0.0,
  fill = c("cornflowerblue", "green", "darkorchid1"), 
  alpha = 0.50,
  cex = 1.5,
  filename="~/Desktop/a.tiff",
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "darkorchid4"),
  #main = list("Possible true possitive identifications"),
  #cat.cex = 1.5, 
  #cat.pos = 0,
  #cat.dist = 0.07,
  cat.fontfamily = "serif",#rotation.degree = 270,
  label.col = "white"
)


# new try



#sapply( gold_t$Merged, FUN = function(   ) )