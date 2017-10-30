match_query_ccl_to_database = function(
  g_query,
  ref_gen = "GRCH37",
  library_name,
  mutational_weight_inclusion_threshold
){
  
  g_mat = read_mutation_grange_objects(
    ref_gen = ref_gen,
    library_name = library_name,
    mutational_weight_inclusion_threshold = mutational_weight_inclusion_threshold
  )
  
  fo_query = GenomicRanges::findOverlaps(
      query = g_query,
      subject = g_mat,
      select = "arbitrary",
      type = "equal"
  )
  match = IRanges::subsetByOverlaps(g_mat, g_query)
  cl_ids = as.character(unlist(str_split(GenomicRanges::mcols(match)$Member_CCLs, pattern = "," )))
  cl_ids = str_replace_all( cl_ids, pattern = paste("_",library_name,sep = ""),"" )
  
  #sort(table(hit_ccls),decreasing = T)
  matching_variants = sort(table(cl_ids), decreasing = T)
  cl_ids_all = str_replace_all( GenomicRanges::mcols(g_mat)$Member_CCLs, pattern = paste("_",library_name,sep = ""),"" )
  cl_ids_all = unique(as.character(unlist(str_split(cl_ids_all,pattern = ","))))
  missing_ccls = table( cl_ids_all[ !( cl_ids_all %in% names(matching_variants) ) ] )
  missing_ccls = missing_ccls - 1
  matching_variants = c(matching_variants, missing_ccls)
  
  matching_variants_t = data.frame(
      CCL = names(matching_variants),
      Matches = as.character(matching_variants),
      Library = rep(library_name, length(matching_variants))
  )
  
  ### add stats info
  
  cl_data <<-  show_contained_cls( verbose = FALSE)
  cl_data = cl_data[cl_data$Library == library_name,]
  
  switch(as.character(mutational_weight_inclusion_threshold),
         "0" = { weight_name = "W0"  },
         "0.25" = { weight_name = "W25"  },
         "0.5" = { weight_name = "W05"  },
         "1" = { weight_name = "W1"  },
         stop(paste( "Could not recocgnize mutational weight ", 
             mutational_weight_inclusion_threshold))
  )
  
  all_muts = as.character(cl_data[[weight_name]])
  ccl_match = match( as.character(matching_variants_t$CCL), as.character(cl_data$CCL), nomatch = 0)
  
  missing_libraries = unique(as.character(matching_variants_t$Library[ccl_match == 0] ) )
  
  if(length(missing_libraries) > 0){
    
    for (mis_lib_name in  missing_libraries){
      
      print( paste("Updating database for library: ", mis_lib_name, sep ="")  )
      mis_mat = read_mutation_grange_objects(
          ref_gen = ref_gen,
          library_name = mis_lib_name,
          mutational_weight_inclusion_threshold = mutational_weight_inclusion_threshold
      )
      write_w0_and_split_w0_into_lower_weights(
          g_mat = mis_mat,
          library_name = mis_lib_name,
          ref_gen = ref_gen
      )
    }
    
    cl_data <<-  show_contained_cls( verbose = FALSE)
    cl_data = cl_data[cl_data$Library == library_name,]
    
    switch(as.character(mutational_weight_inclusion_threshold),
           "0" = { weight_name = "W0"  },
           "0.25" = { weight_name = "W25"  },
           "0.5" = { weight_name = "W05"  },
           "1" = { weight_name = "W1"  },
           stop(paste( "Could not recocgnize mutational weight ", 
                       mutational_weight_inclusion_threshold))
    )
    
    all_muts <<- as.character(cl_data[[weight_name]])
    ccl_match <<- match( as.character(matching_variants_t$CCL), as.character(cl_data$CCL), nomatch = 0)
  }
  
  matching_variants_t$All_variants = all_muts[ccl_match]
  
  return(matching_variants_t)
}
