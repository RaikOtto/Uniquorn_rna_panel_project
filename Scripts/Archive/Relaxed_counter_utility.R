library("Uniquorn")
library("stringr")
library("limma")
library("VennDiagram")

g_t = read.table("~/Dropbox/PhD/Uniquorn_project/Pub/Goldstandard.tab",sep="\t",header=T, stringsAsFactors = F)
#g_t = read.table("~/Dropbox/PhD/Uniquorn_project/Pub/Results/Benchmark/Identification_relaxed.tab",sep="\t",header=T, stringsAsFactors = F)

h   = sapply( g_t$Merged, FUN = str_split, "," )
#h   = sapply( g_t$True_positives, FUN = str_split, "," )
j   = sapply( h, FUN = sort )
label_cls = sapply( j, FUN= paste0,  collapse = "," )

my_fun = function( input_vec ){
  
  input_vec_open = as.character(unlist( str_split( input_vec, ","  ) ))
  res_vec = matrix( numeric(), nrow = 0, ncol= 3 )

  res_vec = matrix( c( 
    length( grep( "COSMIC",    input_vec_open ) ),
    length( grep( "CCLE",      input_vec_open ) ),
    length( grep( "CELLMINER", input_vec_open ) )
  ), ncol = 3 )
  colnames( res_vec ) = c("COSMIC", "CCLE", "CELLMINER")
  rownames( res_vec ) = paste0( c( input_vec_open ), collapse = "," )
  
  return( res_vec )
}

d = t( matrix( sapply( g_t$Merged, FUN = my_fun ), nrow = 3 ) )
#d = t( matrix( sapply( g_t$True_positives, FUN = my_fun ), nrow = 3 ) )

tmp_res = aggregate( d, by = list( label_cls ), FUN = paste  )

splitter_merger_fun = function( input_vector ){
  
  input_vector = input_vector[-1]
  for( i in seq(3)){
    
    input_vector[i] = sum(as.integer( unlist( input_vector[i] ) ))
  }
  input_vector = matrix( as.integer(input_vector), ncol = 3, nrow = 1 )
  return( input_vector )
}

counts = t( apply( tmp_res, FUN = splitter_merger_fun, MARGIN = 1 ) )
rownames(counts) = tmp_res$Group.1
colnames(counts) = c("COSMIC","CCLE","CELLMINER")

counts = as.data.frame(counts)

vennDiagram(counts)

member_cosmic = counts$COSMIC > 0
member_ccle   = counts$CCLE > 0
member_cell   = counts$CELLMINER > 0

member_cosmic_ccle    =     member_cosmic   &     member_ccle   & (! member_cell )
member_cosmic_cell    =     member_cosmic   & ( ! member_ccle ) &    member_cell
member_ccle_cell      = ( ! member_cosmic ) &     member_ccle   &    member_cell

cosmic_exclusive      =     member_cosmic   & ( ! member_ccle ) & (! member_cell )
ccle_exclusive        = ( ! member_cosmic ) &     member_ccle   & (! member_cell )
cell_exclusive        = ( ! member_cosmic ) &  (! member_ccle ) &    member_cell

member_all            =     member_cosmic   &     member_ccle   &    member_cell

# numbers

nr_cosmic_exclusive      = sum( counts[ cosmic_exclusive,] )
nr_ccle_exclusive        = sum( counts[ ccle_exclusive,  ] )
#nr_cell_exclusive        = sum( counts[ cell_exclusive,  ] ); 
nr_cell_exclusive = 0

nr_idents_cosmic_ccle =     sum( counts[ member_cosmic_ccle, ])
nr_idents_cosmic_cell =     sum( counts[ member_cosmic_cell, ])
nr_idents_ccle_cell   =     sum( counts[ member_ccle_cell,   ])

nr_member_all         =     sum( counts[ member_all, ] )

all = nr_cosmic_exclusive + nr_ccle_exclusive + nr_cell_exclusive + nr_idents_cosmic_ccle + nr_idents_cosmic_cell + nr_idents_ccle_cell + nr_member_all 

# venn

member_mat = matrix( integer(), nrow = 0, ncol = 3 )

member_mat = rbind( member_mat, cbind( rep( 1, nr_idents_cosmic_ccle ), rep( 1, nr_idents_cosmic_ccle ), rep( 0, nr_idents_cosmic_ccle ) ) )
member_mat = rbind( member_mat, cbind( rep( 1, nr_idents_cosmic_cell ), rep( 0, nr_idents_cosmic_cell ), rep( 1, nr_idents_cosmic_cell ) ) )

member_mat = rbind( member_mat, cbind( rep( 0, nr_idents_ccle_cell   ), rep( 1, nr_idents_ccle_cell   ), rep( 1, nr_idents_ccle_cell   ) ) )

member_mat = rbind( member_mat, cbind( rep( 1, nr_member_all         ), rep( 1, nr_member_all         ), rep( 1, nr_member_all         ) ) )

cosmic_reps         = nr_cosmic_exclusive #+ ( sum( member_cosmic ) - nr_cosmic_exclusive )
ccle_reps           = nr_ccle_exclusive  # + ( sum( member_ccle   ) - nr_ccle_exclusive   )
cell_reps           = nr_cell_exclusive  # + ( sum( member_cell   ) - nr_cell_exclusive   )

cosmic_substitution = cbind( rep( 1, cosmic_reps ), rep( 0, cosmic_reps ), rep( 0, cosmic_reps ) )
ccle_substitution   = cbind( rep( 0, ccle_reps   ), rep( 1, ccle_reps   ), rep( 0, ccle_reps   ) )
cell_substitution   = cbind( rep( 0, cell_reps   ), rep( 0, cell_reps   ), rep( 1, cell_reps   ) )

#cosmic_self_finder  = cbind( rep( 1, sum_nr_self_intersect_cosmic),rep( 0, sum_nr_self_intersect_cosmic ),rep(0, sum_nr_self_intersect_cosmic )  )
#ccle_self_finder    = cbind( rep( 0, sum_nr_self_intersect_ccle ), rep( 1, sum_nr_self_intersect_ccle ),rep(0, sum_nr_self_intersect_ccle )  )
#cell_self_finder    = cbind( rep( 0, sum_nr_self_intersect_cell ), rep( 0, sum_nr_self_intersect_cell ),rep(1, sum_nr_self_intersect_cell )  )

#member_mat          = rbind( member_mat, cosmic_self_finder, ccle_self_finder, cell_self_finder )
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
