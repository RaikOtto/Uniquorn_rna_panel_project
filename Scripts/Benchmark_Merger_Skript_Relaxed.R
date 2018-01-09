run_relaxed_merging = function( 
  input_path_identification,
  output_table_path
){
    library("stringr")

    gold_std_path     = "~/Uniquorn_rna_panel_project/Misc/Goldstandard.tsv"
    
    input_data        = read.table(input_path_identification, sep ="\t", header =T   )
    gold_std_t        = read.table( 
        gold_std_path,
        sep = "\t",
        header = T,
        fill = T,
        stringsAsFactors = F
    )
    
    identifier       = unique(gold_std_t$Naked_merged)
    identifier_match = match( gold_std_t$Naked_merged, identifier )
    
    comparisons_to_check = aggregate(
      x = gold_std_t$CL,
      by = list( gold_std_t$Naked_merged  ),
      FUN = paste0, sep = "",
      collapse = ","
    )
    
    search_term = as.character(
      input_data$Name_query
    )
    
    map_true_ident = as.character(
        sapply(
            search_term,
            FUN = grep,
            comparisons_to_check$Group.1,
            value = T
        )
    )
    
    identification_successful = aggregate( 
        as.character(
            unlist(
                input_data$True_positive
            )
        ),
        by = list(map_true_ident),
        FUN=paste,
        sep =","
    )
    identifier       = as.character(unlist( identification_successful$Group.1))
    success          = str_replace_all( pattern = "([^A-Z0-9,])", as.character(( identification_successful[,2] ) ), "" )
    expected         = str_replace_all( pattern = "([^A-Z0-9,])", as.character( aggregate( input_data$Expected,       by = list( map_true_ident ), FUN=paste, sep ="," )[,2] ), "" )
    found            = str_replace_all( pattern = "([^A-Z0-9,])", as.character( aggregate( input_data$Found,          by = list( map_true_ident ), FUN=paste, sep ="," )[,2] ), "" )
    true_positives   = str_replace_all( pattern = "([^A-Z0-9,])", as.character( aggregate( input_data$True_positive,  by = list( map_true_ident ), FUN=paste, sep ="," )[,2] ), "" )
    false_negatives  = str_replace_all( pattern = "([^A-Z0-9,])", as.character( aggregate( input_data$False_negative, by = list( map_true_ident ), FUN=paste, sep ="," )[,2] ), "" )
    false_positives  = str_replace_all( pattern = "([^A-Z0-9,])", as.character( aggregate( input_data$False_positive, by = list( map_true_ident ), FUN=paste, sep ="," )[,2] ), "" )
    
    identifier = str_replace_all(pattern = "([^A-Z0-9,])" , identifier, ""  )
    identifier = str_replace_all(pattern = "^(,)*|(,$)*" , identifier, ""  )
    identifier = str_replace_all(pattern = "(,){2,10}" , identifier, ","  )
    identifier = str_replace_all(pattern = "," , identifier, ", "  )
    identifier = str_replace_all(pattern = "COSMIC" , identifier, "_COSMIC"  )
    identifier = str_replace_all(pattern = "CELLMINER" , identifier, "_CELLMINER"  )
    identifier = str_replace_all(pattern = "CCLE" , identifier, "_CCLE"  )
    
    splitit = sapply( identifier, FUN = str_split, ", "  )
    #splitit = sapply( splitit, FUN = unique  )
    identifier = as.character(sapply( splitit, FUN = paste0, collapse = ", "  ))
    
    false_negatives = str_replace_all(pattern = "^(,)*|(,$)*" , false_negatives, ""  )
    false_negatives = str_replace_all(pattern = "(,){2,10}" , false_negatives, ","  )
    false_negatives = str_replace_all(pattern = "," , false_negatives, ", "  )
    false_negatives = str_replace_all(pattern = "COSMIC" , false_negatives, "_COSMIC"  )
    false_negatives = str_replace_all(pattern = "CELLMINER" , false_negatives, "_CELLMINER"  )
    false_negatives = str_replace_all(pattern = "CCLE" , false_negatives, "_CCLE"  )
    
    splitit = sapply( false_negatives, FUN = str_split, ", "  )
    #splitit = sapply( splitit, FUN = unique  )
    false_negatives = sapply( splitit, FUN = paste0, collapse = ", "  )
    
    false_positives = str_replace_all(pattern = "^(,)*|(,$)*" , false_positives, ""  )
    false_positives = str_replace_all(pattern = "(,){2,10}" , false_positives, ","  )
    false_positives = str_replace_all(pattern = "," , false_positives, ", "  )
    false_positives = str_replace_all(pattern = "COSMIC" , false_positives, "_COSMIC"  )
    false_positives = str_replace_all(pattern = "CELLMINER" , false_positives, "_CELLMINER"  )
    false_positives = str_replace_all(pattern = "CCLE" , false_positives, "_CCLE"  )
    
    splitit = sapply( false_positives, FUN = str_split, ", "  )
    #splitit = sapply( splitit, FUN = unique  )
    false_positives = sapply( splitit, FUN = paste0, collapse = ", "  )
    
    found = str_replace_all(pattern = "^(,)*|(,$)*" , found, ""  )
    found = str_replace_all(pattern = "(,){2,10}" , found, ","  )
    found = str_replace_all(pattern = "," , found, ", "  )
    found = str_replace_all(pattern = "COSMIC" , found, "_COSMIC"  )
    found = str_replace_all(pattern = "CELLMINER" , found, "_CELLMINER"  )
    found = str_replace_all(pattern = "CCLE" , found, "_CCLE"  )
    
    splitit = sapply( found, FUN = str_split, ", "  )
    #splitit = sapply( splitit, FUN = unique  )
    found = sapply( splitit, FUN = paste0, collapse = ", "  )
    
    expected = str_replace_all(pattern = "^(,)*|(,$)*" , expected, ""  )
    expected = str_replace_all(pattern = "(,){2,10}" , expected, ","  )
    expected = str_replace_all(pattern = "," , expected, ", "  )
    expected = str_replace_all(pattern = "COSMIC" , expected, "_COSMIC"  )
    expected = str_replace_all(pattern = "CELLMINER" , expected, "_CELLMINER"  )
    expected = str_replace_all(pattern = "CCLE" , expected, "_CCLE"  )
    
    splitit = sapply( expected, FUN = str_split, ", "  )
    #splitit = sapply( splitit, FUN = unique  )
    expected = sapply( splitit, FUN = paste0, collapse = ", "  )
    
    true_positives = str_replace_all(pattern = "^(,)*|(,$)*" , true_positives, ""  )
    true_positives = str_replace_all(pattern = "(,){2,10}" , true_positives, ","  )
    true_positives = str_replace_all(pattern = "," , true_positives, ", "  )
    true_positives = str_replace_all(pattern = "COSMIC" , true_positives, "_COSMIC"  )
    true_positives = str_replace_all(pattern = "CELLMINER" , true_positives, "_CELLMINER"  )
    true_positives = str_replace_all(pattern = "CCLE" , true_positives, "_CCLE"  )
    
    splitit = sapply( true_positives, FUN = str_split, ", "  )
    #splitit = sapply( splitit, FUN = unique  )
    true_positives = sapply( splitit, FUN = paste0, collapse = ", "  )
    
    outcome_table = as.data.frame(
      cbind(
        identifier,
        success,
        expected,
        found,
        true_positives,
        false_negatives,
        false_positives
      )
    )
    colnames(outcome_table) = c(
      "Identifier",
      "Success",
      "Expected",
      "Found",
      "True_positives",
      "False_negatives",
      "False_positives"
    )
    
    outcome_table = unique(outcome_table)
    outcome_table$Success = as.character(grepl("TRUE",outcome_table$Success))
    
    write.table( outcome_table, output_table_path,sep="\t",quote =F, row.names = F)
    
    table( outcome_table$Success  )
}
