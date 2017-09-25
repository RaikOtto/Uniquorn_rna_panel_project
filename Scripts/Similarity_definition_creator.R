library("Uniquorn")
library("stringr")

known_pairs_t = read.table("~/Uniquorn_data//benchmark_vcf_files//known_relationships.tab",sep ="\t", header =T, fill=T,stringsAsFactors = F)

c = show_amount_cls_per_database()

ccle = c[grepl("CCLE",c[,1]),1]
cosmic = c[grepl("COSMIC",c[,1]),1]
cellminer = c[grepl("CELLMINER",c[,1]),1]

# ega

ega_CCLs = list.files("~/Uniquorn_data/EGA_results/ident/", pattern = ".ident")
ega_CCLs = str_replace(ega_CCLs, pattern = ".hg19.vcf.ident","")
ega_CCLs = str_to_upper(ega_CCLs)
ega_CCLs = sapply( ega_CCLs, FUN = paste, "EGA", sep ="_" )

# gdc

gdc_CCLs = list.files("~/Uniquorn_data/GDC_results/", pattern = ".ident")
gdc_CCLs = str_replace(gdc_CCLs, pattern = ".1.hg19.vcf.ident","")
gdc_CCLs = str_replace(gdc_CCLs, pattern = ".2.hg19.vcf.ident","")
gdc_CCLs = str_to_upper(gdc_CCLs)
gdc_CCLs = sapply( gdc_CCLs, FUN = paste, "GDC", sep ="_" )

# all

all_cls = c(ega_CCLs, gdc_CCLs,c[grepl("COSMIC",c[,1]),1], c[grepl("CCLE",c[,1]),1],c[grepl("CELLMINER",c[,1]),1]  )
all_cls = all_cls[all_cls != ""]
identifier = str_replace(all_cls,"_COSMIC","")
identifier = str_replace(identifier,"_CCLE","")
identifier = str_replace(identifier,"_CELLMINER","")
identifier = str_replace(identifier,"_CUSTOM","")
identifier = str_replace(identifier,"-","")
identifier = str_replace(identifier,"-","")
identifier = str_replace(identifier,"_","")

known_cls = unique(c(known_pairs_t$Cl1,known_pairs_t$Cl2))

make_name_comparable = function( replace_name ){
  
  identical_names = c("")
  replace_name = as.character( replace_name )
  replace_name = str_to_upper( replace_name )
  
  if (replace_name == "TT_CCLE"  ) 
    replace_name = "TTALT1_CCLE"
  
  if (replace_name == "TT_COSMIC"  ) 
    replace_name = "TTALT1_COSMIC"
  
  if (replace_name == "T-T_COSMIC"  ) 
    replace_name = "TTALT2_COSMIC"
  
  if (replace_name == "KMH-2_COSMIC"  ) 
    replace_name = "KMH2ALT1_COSMIC"
  
  if (replace_name == "KM-H2_COSMIC"  ) 
    replace_name = "KMH2ALT2_COSMIC"
  
  replace_name = as.character( str_replace_all( string = replace_name, pattern = "(_COSMIC)|(_CCLE)|(_CELLMINER)|(EGA)|(GDC)", "" ) )
  replace_name = as.character( str_replace_all( string = replace_name, pattern = "([^A-Z0-9])*", "" ) )
  
  return( replace_name )
}

stripped_names = sapply( all_cls, FUN = make_name_comparable )
names(stripped_names) = all_cls
# similarity due to identical names

identify_similar_names = function( index_cl ){
  
  stripped_name = stripped_names[ index_cl ]
  original_name = names(stripped_names[ index_cl ])
  
  sim_map       = which( as.character( stripped_names ) == stripped_name)
  similar_cls   = paste0( c( names( stripped_names )[ sim_map ] ), collapse = "," )
  
  return( similar_cls )
}

all_cls = cbind( all_cls, sapply( seq( length( all_cls ) ), FUN = identify_similar_names ))
dim(all_cls)
### similarity due to reports

known_pairs   = paste( known_pairs_t[,1], known_pairs_t[,2], sep ="_"  )
known_pairs   = c( known_pairs, paste( known_pairs_t[,2], known_pairs_t[,1], sep ="_"  ) )
all_cls       = cbind( all_cls, rep("", length(all_cls[,1])) )

identify_related_names = function( index_cl ){
  
  stripped_name = as.character( stripped_names[ index_cl ] )
  stripped_name = str_replace( stripped_name, pattern = "CUSTOM$", "" )
  #original_name = names(stripped_names[ index_cl ])
  related_cls     = character()
  ori_related_cls = character()
  
  related_map   = which( grepl( paste0( c( "(^",stripped_name,"_)+"), collapse = "" ), known_pairs ) )
  
  for (map in related_map){
    
    other_cl = tail(as.character( unlist( str_split(known_pairs[map], "_") ) ), 1)
    grep_other_cl = paste0( c( "^",other_cl,"$"), collapse = "" )
    
    ori_names = names(
      stripped_names[ grepl(
        grep_other_cl,
        stripped_names
      )]
    )
    
    for (o_name in ori_names){
      ori_related_cls = c( ori_related_cls, o_name )
    }
  }
  
  related_cls = paste0( unique( ori_related_cls ), collapse = "," )
  
  return( related_cls )
}
all_cls[ , 3] = sapply( seq( length( all_cls[,1] ) ), FUN = identify_related_names )

colnames(all_cls) = c("CL","Name_identical","Related")

write.table(x = all_cls,"~/Dropbox/PhD/Uniquorn_project/Pub/Goldstandard.tab",sep="\t",quote =F , row.names=F)

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

colnames(all_cls) = c("CL","Name_identical","Related","Merged","Naked_merged")

write.table(x = all_cls,"~/Dropbox/PhD/Uniquorn_project/Pub/Goldstandard.tsv",sep="\t",quote =F , row.names=F)
