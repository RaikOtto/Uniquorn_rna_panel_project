
### venn diagrams

library("Uniquorn")
library("stringr")
library("limma")
library("VennDiagram")

### gold

gold_t = read.table("~/Dropbox/PhD/Uniquorn_project/Pub/Goldstandard.tab",sep ="\t", header = T, stringsAsFactors = F)

intersects     = as.character( gold_t$Merged )
naked_label    = as.character( gold_t$Naked_merged )
naked_label_nr = unique(as.character( unlist( str_split( naked_label, "," )  ) ) ) # 1341


# absolute numbers

seen_d <<- as.list( rep( F, length( gold_t$CL ) ) )
names( seen_d ) = as.character( gold_t$CL )

clean_intersects = c()

for (i in 1:dim(gold_t)[1] ) {
  
  print (i)
  entry = as.character( gold_t[i, c(1,4)] )
  entry = as.character( unlist( str_replace(entry, "\\[", "")))
  entry = as.character( unlist( str_replace(entry, "\\]", "")))
  
  res = as.character(
        unlist(
          str_replace(
            entry[2],
            pattern = entry[1],
            ""
          )
        )
    )
    clean_intersects  <<- c( clean_intersects, res )
}

counter_function = function( intersects ){
  
  expt_mat = as.vector( rep( 0, 7 ) )
  names( expt_mat ) = c("COSMIC","CCLE","CELL","COSMIC_CCLE","COSMIC_CELL","CCLE_CELL","ALL")
  
  for ( i in 1:length(intersects ) ) {

    entry = intersects[i]
    
    cosmic_in_query = grepl( "COSMIC",    gold_t$CL[i] )
    ccle_in_query   = grepl( "CCLE",      gold_t$CL[i] )
    cell_in_query   = grepl( "CELLMINER", gold_t$CL[i] )
    
    cosmic_in_entry = grepl( "COSMIC",    entry )
    ccle_in_entry   = grepl( "CCLE",      entry )
    cell_in_entry   = grepl( "CELLMINER", entry )
    
    cosmic_ccle     =     cosmic_in_entry   &     ccle_in_entry   & ( ! cell_in_entry )
    cosmic_cell     =     cosmic_in_entry   & ( ! ccle_in_entry ) &     cell_in_entry 
    ccle_cell       = ( ! cosmic_in_entry ) &     ccle_in_entry   &     cell_in_entry 
    all             =     cosmic_in_entry   &     ccle_in_entry   &     cell_in_entry

    if ( cosmic_in_query ){
      
        expt_mat["COSMIC"] = expt_mat["COSMIC"] + str_count( entry, pattern =  "COSMIC" )
      
        if ( cosmic_ccle ){
          
            expt_mat["COSMIC_CCLE"] = expt_mat["COSMIC_CCLE"] + str_count( entry, pattern =  "CCLE" )
            
        } else if (cosmic_cell) {
          
            expt_mat["COSMIC_CELL"] = expt_mat["COSMIC_CELL"] + str_count( entry, pattern =  "CELLMINER" )  
            
        } else if (all){
          
            expt_mat["ALL"] = expt_mat["ALL"] + str_count( entry, pattern =  "CCLE" )
            expt_mat["ALL"] = expt_mat["ALL"] + str_count( entry, pattern =  "CELLMINER" )  
        }
    } else if ( ccle_in_query ){
      
      expt_mat["CCLE"] = expt_mat["CCLE"] + str_count( entry, pattern =  "CCLE" )
      
      if ( ccle_cell ){
        
          expt_mat["CCLE_CELL"] = expt_mat["CCLE_CELL"] + str_count( entry, pattern =  "CELLMINER" )
          
      } else if (cosmic_ccle) {
        
        expt_mat["COSMIC_CCLE"] = expt_mat["COSMIC_CCLE"] + str_count( entry, pattern =  "COSMIC" )
        
      } else if (all){
        
        expt_mat["ALL"] = expt_mat["ALL"] + str_count( entry, pattern =  "COSMIC" )
        expt_mat["ALL"] = expt_mat["ALL"] + str_count( entry, pattern =  "CELLMINER" )
      }
    } else if ( cell_in_query ){
      
      expt_mat["CELL"] = expt_mat["CELL"] + str_count( entry, pattern =  "CELLMINER" )
      
      if ( ccle_cell ){
        
        expt_mat["CCLE_CELL"] = expt_mat["CCLE_CELL"] + str_count( entry, pattern =  "CCLE" )
        
      } else if (cosmic_cell) {
        
        expt_mat["COSMIC_CCLE"] = expt_mat["COSMIC_CCLE"] + str_count( entry, pattern =  "COSMIC" )
        
      } else if (all){
        
        expt_mat["ALL"] = expt_mat["ALL"] + str_count( entry, pattern =  "COSMIC" )
        expt_mat["ALL"] = expt_mat["ALL"] + str_count( entry, pattern =  "CCLE" )
      }
    } 
  }
    
  return(expt_mat)
}

c_vec = counter_function( intersects )

# intersect

member_mat = matrix( integer(), nrow = 0, ncol = 3 )

member_mat = rbind( member_mat, cbind( rep( 1, c_vec["COSMIC"]      ), rep( 0, c_vec["COSMIC"]      ), rep( 0, c_vec["COSMIC"] ) ) )
member_mat = rbind( member_mat, cbind( rep( 0, c_vec["CCLE"]        ), rep( 1, c_vec["CCLE"]        ), rep( 0, c_vec["CCLE"] ) ) )
member_mat = rbind( member_mat, cbind( rep( 0, c_vec["CELL"]        ), rep( 0, c_vec["CELL"]        ), rep( 1, c_vec["CELL"] ) ) )

member_mat = rbind( member_mat, cbind( rep( 1, c_vec["COSMIC_CCLE"] ), rep( 1, c_vec["COSMIC_CCLE"] ), rep( 0, c_vec["COSMIC_CCLE"] ) ) )
member_mat = rbind( member_mat, cbind( rep( 1, c_vec["COSMIC_CELL"] ), rep( 0,  c_vec["COSMIC_CELL"] ), rep( 1,  c_vec["COSMIC_CELL"] ) ) )

member_mat = rbind( member_mat, cbind( rep( 0,  c_vec["CCLE_CELL"]   ), rep( 1, c_vec["CCLE_CELL"]   ), rep( 1, c_vec["CCLE_CELL"] ) ) )

member_mat = rbind( member_mat, cbind( rep( 1, c_vec["ALL"]         ), rep( 1,  c_vec["ALL"]         ), rep( 1,  c_vec["ALL"]       ) ) )

member_mat          = data.frame( member_mat )

colnames(member_mat) = c("COSMIC","CCLE","CELLMINER")

vennDiagram( member_mat )

d = venn.diagram(
  x = list(
    COSMIC = which( member_mat$COSMIC == 1 ),
    CCLE = which( member_mat$CCLE == 1 ),
    CELLMINER = which( member_mat$CELLMINER == 1 )
  ), 
  height=2000, width=2000, resolution=300,
  col = "transparent",
  margin = 0.0,
  fill = c("blue","black", "red"), 
  alpha = .66,
  cex = 1.5,
  filename="~/Downloads//a.tiff",
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "darkorchid4"),
  lwd = 4,
  lty = 1,
  cat.fontfamily = "serif",#rotation.degree = 270,
  label.col = "white",
  euler.d = F,
  scaled = F,
  print.mode = 'percent'
)

