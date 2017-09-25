library("ggplot2")
library("stringr")
library("Uniquorn")


c = show_contained_mutations()

group = rep( "", dim(c)[1] )
group[ grepl( x = c$CL, pattern = "_COSMIC") ] = "CoSMIC"
group[ grepl( x = c$CL, pattern = "_CCLE") ] = "CCLE"
group[ grepl( x = c$CL, pattern = "_CELLMINER") ] = "CellMiner"

wf_t = data.frame( cbind(c$Weight, c$CL, group) )
colnames( wf_t ) = c("Weight","CL","Panel")
wf_t$Weight = as.character( wf_t$Weight )

## matrix calculation

agg_res = function( selector, mat ){
  
  mat = mat[ mat$Panel == selector,]
  print(selector)
  
  mat$Weight[ as.double(mat$Weight) <= .25] = .25
  mat$Weight[ ( as.double(mat$Weight) <= .5 ) & (as.double(mat$Weight) > .25)] = .5
  mat$Weight[ as.double(mat$Weight) > .5] = 1
  
  interval_mat <<- rbind(
    interval_mat,
    cbind( 
      as.double( mat$Weight ),
      as.character( mat$CL ),
      rep(
        selector, 
        length(mat$Weight ) 
      ) 
    )
  )
}


interval_mat <<- matrix( as.double(), nrow = 0, ncol = 3 )
lapply( as.character( unique( wf_t$Panel  ) ), FUN = agg_res, wf_t )

interval_mat = data.frame( matrix( as.character( unlist( interval_mat ) ), ncol = 3 ) )
colnames(interval_mat) = c("Weight","CL","Group")

## matrix compact

compact_mat <<- matrix( as.double(), nrow = 0, ncol = 3 )

agg_res = function( selector, mat ){
  
  mat = mat[ mat$Panel == selector,]
  print(selector)
  
  h = hist( as.double( mat$Weight ), breaks = c(0,.25,.5,1), plot = F  )
  vals = h$counts
  
  compact_mat <<- rbind( compact_mat,
     cbind( h$breaks[-1],
            vals,
            rep(selector, length( vals ) )
     )
  )
}

lapply( as.character( unique( wf_t$Panel  ) ), FUN = agg_res, wf_t )

compact_mat = data.frame( matrix( as.character( unlist( compact_mat ) ), ncol = 3 ) )
colnames(compact_mat) = c("Weight","Count","Group")
compact_mat$Group = factor(compact_mat$Group, levels = unique(compact_mat$Group))

## cancer cell line calculation unique

cl_mat_calc_unique = function( selector, mat, groups ){
  
  mat = mat[ groups == selector,]
  mat = mat[ as.character( mat$Weight ) == "1", ]
  sum_vec = rep(1, length( mat$Fingerprint ) )
  print(selector)
  
  vals = aggregate( sum_vec, FUN = sum, by = list( as.character( mat$CL ) ) )
  
  cl_mat_unique <<- rbind( cl_mat_unique,
    cbind( vals$x,
       vals$Group.1,
       rep(selector, length( unique( vals$Group.1) ) )
    )
  )
}

cl_mat_unique <<- matrix( as.double(), nrow = 0, ncol = 3 )
lapply( as.character( unique( group  ) ), FUN = cl_mat_calc_unique, mat = c,  groups = group )

cl_mat_unique = data.frame( matrix( as.character( unlist( cl_mat_unique ) ), ncol = 3 ), stringsAsFactors = F )
colnames( cl_mat_unique ) = c("Count","CL","Group")
#cl_mat_unique$Group = factor( group, levels = unique( group))

## cancer cell line calculation all

cl_mat_calc_all = function( selector, mat, groups ){
  
  mat = mat[ groups == selector,]

  sum_vec = rep(1, length( mat$Fingerprint ) )
  print(selector)
  
  vals = aggregate( sum_vec, FUN = sum, by = list( as.character( mat$CL ) ) )
  
  cl_mat_all <<- rbind( cl_mat_all,
                           cbind( vals$x,
                                  vals$Group.1,
                                  rep(selector, length( unique( vals$Group.1) ) )
                           )
  )
}

cl_mat_all <<- matrix( as.double(), nrow = 0, ncol = 3 )
lapply( as.character( unique( group  ) ), FUN = cl_mat_calc_all, mat = c,  groups = group )

cl_mat_all = data.frame( matrix( as.character( unlist( cl_mat_all ) ), ncol = 3 ), stringsAsFactors = F )
colnames( cl_mat_all ) = c("Count","CL","Group")

### repair missing

missing_cls = as.character( cl_mat_all$CL[ which( !( cl_mat_all$CL %in% cl_mat_unique$CL )  ) ] )
missing_mat  = cbind( c(0,0), missing_cls, c("CCLE","CCLE" ) )
colnames(missing_mat) = colnames(cl_mat_unique)
cl_mat_unique = rbind( cl_mat_unique, missing_mat )

cl_mat_all = cl_mat_all[ order( as.character(  cl_mat_all$CL ) ), ]
cl_mat_unique = cl_mat_unique[ order( as.character( cl_mat_unique$CL ) ), ]

dif = as.integer( cl_mat_all$Count ) - as.integer( cl_mat_unique$Count )
dif

cl_mat_unique$Group = paste( cl_mat_all$Group, "filtered", sep = " " )
cl_mat_all$Group = paste( cl_mat_all$Group, "raw", sep = " " )

cl_mat_merge = rbind( cl_mat_all, cl_mat_unique )
library(tidyr)
cl_mat_merge = tidyr::separate(cl_mat_merge, Group, into = c("Panel", "Subtype", "Type"), sep = " ", remove = F  )

## mat unique average

cl_mat_average_all = function( selector, mat, weights ){
  
  cl_mat = mat[ mat$CL == selector,]
  
  for ( weight in weights ){
    
    sub_mat = cl_mat[ cl_mat$Weight == weight, ]
    print(c(selector,weight))
    
    weight_mat = cbind(
      length( sub_mat$Weight ),
      as.character( weight ),
      as.character( selector ),
      tail( as.character( cl_mat$Group ), 1 )
    )
    
    cl_count_mat <<- rbind( 
      cl_count_mat,
      weight_mat
    )
  }
}

cl_count_mat <<- matrix( as.character(), nrow = 0, ncol = 4 )
lapply( unique( as.character( interval_mat$CL  ) ), FUN = cl_mat_average_all, interval_mat,  unique( as.character( interval_mat$Weight ) ) )

cl_count_mat = data.frame( matrix( as.character( unlist( cl_count_mat ) ), ncol = 4 ), stringsAsFactors = F )
colnames( cl_count_mat ) = c("Count","Weight","CL","Group")

# abs_mat 

cl_mat_average_cl = function( selector, mat, weights ){
  
  cl_mat = mat[ mat$Group == selector,]
  cl_mat$Weight = as.character( cl_mat$Weight  )
  
  list_cls = unique(cl_mat$CL)
  
  for ( cl in list_cls){
    
    sub_mat_cl = cl_mat[ cl_mat$CL == cl,]
    nr_all_muts_weight = length( sub_mat_cl$Weight )
    
    for ( weight in weights ){
      
      sub_mat = sub_mat_cl[ as.double( sub_mat_cl$Weight ) == as.double( weight ), ]
      print(c(selector, as.character( cl ),weight))
      
      nr_muts_weight = length( sub_mat$Weight )
      
      weight_mat = cbind(
        as.character( nr_muts_weight ),
        as.character( weight ),
        as.character( cl ),
        as.character( selector )
      )
      
      cl_abs_count_mat <<- rbind( 
        cl_abs_count_mat,
        weight_mat
      )
    }
  }
}

interval_mat$Weight = as.character(interval_mat$Weight)
#interval_mat$Weight[ interval_mat$Weight == 3] = "1.0"
#interval_mat$Weight[ interval_mat$Weight == 2] = "0.5"
#interval_mat$Weight[ interval_mat$Weight == 1] = "0.25"

cl_abs_count_mat <<- matrix( as.character(), nrow = 0, ncol = 4 )
lapply( unique( as.character( interval_mat$Group  ) ), FUN = cl_mat_average_cl, interval_mat,  unique( as.character( interval_mat$Weight ) ) )

cl_abs_count_mat = data.frame( matrix( as.character( unlist( cl_abs_count_mat ) ), ncol = 4 ), stringsAsFactors = F )
colnames( cl_abs_count_mat ) = c("Count","Weight","CL","Group")
cl_abs_count_mat_save = cl_abs_count_mat
