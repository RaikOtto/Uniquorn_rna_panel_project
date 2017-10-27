
### venn diagrams

library("Uniquorn")
library("stringr")
library("limma")
library("VennDiagram")

### gold

gold_t = read.table("~/Dropbox/PhD/Uniquorn_project/Pub/Goldstandard.tab",sep ="\t", header = T, stringsAsFactors = F)

intersects     = as.character( gold_t$Merged )
naked_label    = as.character( gold_t$Naked_merged )
naked_label_nr = unique(as.character( unlist( str_split( naked_label, "," )  ) ) ) # 1342

# absolute numbers

counter_function = function( intersects ){
  
    seen_d <<- as.list( rep( F, length( gold_t$CL ) ) )
    names( seen_d ) = as.character( gold_t$CL )
    
    expt_mat = as.vector( rep( 0, 10 ) )
    names( expt_mat ) = c(
        "COSMIC_only",
        "CCLE_only",
        "CELL_only",
        "COSMIC_CCLE",
        "COSMIC_CELL",
        "CCLE_CELL",
        "ALL",
        "COSMIC_MULTI",
        "CCLE_MULTI",
        "CELL_MULTI"
    )
    
    for ( i in 1:length( intersects ) ) {
      
        entry_line_raw = intersects[i]
        queries        = as.character( unlist( str_split( entry_line_raw, "," )  ) )
        
        for (entry in queries){
          
            if ( entry == ""  )
                break
          
            print(entry)
            
            if ( ! as.logical( seen_d[[ entry ]]  ) ){
            
                seen_d[[ entry ]] = T
            
            } else {
            
                entry_line_raw = str_replace( entry_line_raw, pattern = entry, ""  )
            }
            
        }
      
        cosmic_in_entry = grepl( "COSMIC",    entry_line_raw )
        ccle_in_entry   = grepl( "CCLE",      entry_line_raw )
        cell_in_entry   = grepl( "CELLMINER", entry_line_raw )
        
        cosmic_only     =     cosmic_in_entry   & ( ! ccle_in_entry ) & ( ! cell_in_entry )
        cosmic_ccle     =     cosmic_in_entry   &     ccle_in_entry   & ( ! cell_in_entry )
        cosmic_cell     =     cosmic_in_entry   & ( ! ccle_in_entry ) &     cell_in_entry 
        ccle_only       = ( ! cosmic_in_entry ) &     ccle_in_entry   & ( ! cell_in_entry )
        ccle_cell       = ( ! cosmic_in_entry ) &     ccle_in_entry   &     cell_in_entry 
        cell_only       = ( ! cosmic_in_entry ) & ( ! ccle_in_entry ) &     cell_in_entry
        all             =     cosmic_in_entry   &     ccle_in_entry   &     cell_in_entry
        
        if (cosmic_only){
          
          expt_mat["COSMIC_only"]  = expt_mat["COSMIC_only"] + 1
          expt_mat["COSMIC_MULTI"] = expt_mat["COSMIC_MULTI"] + str_count(  "COSMIC" ) - 1
          
        } else if (ccle_only) {
          
            expt_mat["CCLE_only"] = expt_mat["CCLE_only"] + 1
            expt_mat["CCLE_MULTI"] = expt_mat["CCLE_MULTI"] + str_count(  "CCLE" ) - 1
          
        } else if (cell_only) {
          
            expt_mat["CELL_only"] = expt_mat["CELL_only"] + 1
            expt_mat["CELL_MULTI"] = expt_mat["CELL_MULTI"] + str_count(  "CELLMINER" ) - 1
          
        } else if (cosmic_ccle){
          
          expt_mat["COSMIC_CCLE"] = expt_mat["COSMIC_CCLE"] + 1
          
        } else if (cosmic_cell) {
          
          expt_mat["COSMIC_CELL"] = expt_mat["COSMIC_CELL"] + 1
          
        } else if (ccle_cell) {
          
          expt_mat["CCLE_CELL"] = expt_mat["CCLE_CELL"] + 1
          
        } else if (all) {
        
          expt_mat["ALL"] = expt_mat["ALL"] + 1
          
        }
    }
    
  return(expt_mat)
}

c_vec = counter_function( intersects )

# intersect

member_mat = matrix( integer(), nrow = 0, ncol = 3 )

member_mat = rbind( member_mat, cbind( rep( 1, c_vec["COSMIC_only"]      ), rep( 0, c_vec["COSMIC_only"]      ), rep( 0, c_vec["COSMIC_only"] ) ) )
member_mat = rbind( member_mat, cbind( rep( 0, c_vec["CCLE_only"]        ), rep( 1, c_vec["CCLE_only"]        ), rep( 0, c_vec["CCLE_only"] ) ) )
member_mat = rbind( member_mat, cbind( rep( 0, c_vec["CELL_only"]        ), rep( 0, c_vec["CELL_only"]        ), rep( 1, c_vec["CELL_only"] ) ) )

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
  print.mode = 'raw'
)

# CoSMIC CLP 1027 - 1024 = 4 shared; 1238/2+264/6+20/2 = 673 shared; 351 self-confined
# CCLE 950 - 904 = 46 shared; 1238/2+264/6+2/2 = 664 shared; 240 self-confined
# CELL 72 - 60 = 12 shared; 20/2 + 264/6 + 2/2 = 55 shared; 5 self-confined