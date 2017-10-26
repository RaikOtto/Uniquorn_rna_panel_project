library("devtools")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()
library("stringr")
ref_gen = "GRCH37"

known_pairs_t = read.table("~/Uniquorn_rna_panel_project/Misc/known_relationships.tsv",sep ="\t", header =T, fill=T,stringsAsFactors = F)

contained_ccls = show_contained_cls()

library_names = Uniquorn::read_library_names(ref_gen = ref_gen)

for (library_name in library_names)
    if ( str_detect(pattern = library_name, b_file) )
        source_file <<- library_name

raw_files_path = "~/Uniquorn_data/benchmark_vcf_files/panel_raw_files/"
i_files = list.files( raw_files_path, pattern = ".vcf$", full.names = T )

all_cls = as.character( sapply(i_files, FUN = function(vec){return(
  tail(as.character(unlist(str_split(vec,"/"))),1))}))

lib_vec <<- c()
for (cls_name in all_cls){
    for (lib_name in library_names){
        if ( str_detect(cls_name, pattern = lib_name))
            lib_vec = c(lib_vec, lib_name)
    }
}

all_cls = str_replace_all(all_cls, pattern = 
    paste(c(library_names,"\\.","vcf"), collapse = "|",sep = ""),"")

identifier_plane = str_replace_all(all_cls, pattern = 
    paste(c("\\_","\\-"), collapse = "|",sep = ""),
    "")
identifier_plane = str_to_upper(identifier_plane)
# all

known_cls = unique(c(known_pairs_t$Cl1,known_pairs_t$Cl2))

for (i in 1:length(all_cls)){
  
    ccl_name = all_cls[i]
    lib_name = lib_vec[i]
  
    if (str_detect( ccl_name, pattern = "TT.CCLE" ))
        identifier_plane[i] = "TTALT1"
  
    if (str_detect( ccl_name, pattern = "TT.COSMIC" ))
        identifier_plane[i] = "TTALT1"
    
    if (str_detect( ccl_name, pattern = "T-T.COSMIC" ))
        identifier_plane[i] = "TTALT2"

    if (str_detect( ccl_name, pattern = "KMH-2.COSMIC" ))
        identifier_plane[i] = "KMH2ALT1"

    if (str_detect( ccl_name, pattern = "KM-H2.COSMIC" ))
        identifier_plane[i] = "KMH2ALT2"
}

# similarity due to identical names

identify_similar_names = function( index_ccl ){
  
  sim_map       = which( as.character( identifier_plane ) == as.character( identifier_plane[index_ccl]) )
  ident_ccls    = paste(all_cls[ sim_map ], lib_vec[ sim_map ], sep ="_" )
  similar_cls   = paste0( c(ident_ccls  ), collapse = "," )
  
  return( similar_cls )
}

all_cls = cbind( 
  paste( all_cls, lib_vec, sep ="_" ),
  identifier_plane,
  sapply( seq( length( all_cls ) ), FUN = identify_similar_names ))

### similarity due to reports

known_pairs   = paste( known_pairs_t[,1], known_pairs_t[,2], sep ="_"  )
known_pairs   = c( known_pairs, paste( known_pairs_t[,2], known_pairs_t[,1], sep ="_"  ) )
all_cls       = cbind( all_cls, rep("", length(all_cls[,1])) )

identify_related_names = function( index_ccl ){
  
    stripped_name = as.character( identifier_plane[index_ccl] )

    related_cls     = character()
    ori_related_cls = character()
  
    related_map   = which( grepl( paste0( c( "(^", stripped_name  ,"_)+"), collapse = "" ), known_pairs ) )
  
    for (map in related_map){
    
        other_cl = tail(as.character( unlist( str_split(known_pairs[map], "_") ) ), 1)
        ident_ccls = identifier_plane[ identifier_plane == other_cl]
        ident_libraries = lib_vec[ identifier_plane == other_cl]
    
        ori_names = paste(
            ident_ccls,
            ident_libraries, sep = "_"
        )
        
        for (o_name in ori_names){
            ori_related_cls = c( ori_related_cls, o_name )
        }
    }
  
    related_cls = paste0( unique( ori_related_cls ), collapse = "," )
  
    return( related_cls )
}
all_cls[ , 4] = sapply( seq( length( all_cls[,1] ) ), FUN = identify_related_names )

colnames(all_cls) = c("CL","CL_plane","Name_identical","Related")

write.table(x = all_cls,"~//Uniquorn_rna_panel_project/Misc//Goldstandard.tsv",sep="\t",quote =F , row.names=F)

### two-sided 
all_cls = as.data.frame(all_cls)
#all_cls       = cbind( all_cls, rep("", length(all_cls[,1])) )

merged_sim_rel = mapply(
    c,
    str_split( all_cls$Name_identical, "," ),
    str_split( all_cls$Related, "," ),
    SIMPLIFY=FALSE
)

clean_me = function( a ){ a=as.character(unlist(a)); a = a[a != ""]; return(a) }
merged_sim_rel = lapply( FUN = clean_me, merged_sim_rel )
merged_sim_rel = lapply( FUN = sort, merged_sim_rel )

naked_merge = lapply( FUN = str_replace_all, merged_sim_rel, pattern = "(_CCLE*)|(_COSMIC*)|(_CELLMINER*)|(_EGA*)|(_GDC*)", replacement = ""  )
naked_merge = lapply( FUN = str_replace_all, naked_merge,    pattern = "([^A-Z0-9,])*", replacement = ""  )
naked_merge = sapply( FUN = str_split, naked_merge, pattern = ","  )
naked_merge = lapply( FUN = unique, naked_merge )

merged_sim_rel = lapply( FUN = paste0, merged_sim_rel, collapse = ","  )
naked_merge    = lapply( FUN = paste0, naked_merge, collapse = ","  )
naked_merge    = str_replace( naked_merge, pattern = "CUSTOM", "" )

all_cls = cbind(all_cls, as.character( unlist( merged_sim_rel ) ) )
all_cls = cbind(all_cls, as.character( unlist( naked_merge ) ) )

colnames(all_cls) = c("CL","CL_plane","Name_identical","Related","Merged","Naked_merged")

write.table(x = all_cls,"~/Uniquorn_rna_panel_project/Misc//Goldstandard.tsv",sep="\t",quote =F , row.names=F)
