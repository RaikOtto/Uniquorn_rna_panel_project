library("stringr")
library(devtools)
load_all()

output_directory = "~/Uniquorn_data/benchmark_vcf_files/raw_files/"
#output_directory = "/Users/raikotto/Dropbox/PhD/MapTor-Net/Andranik/First_draft_RNA_var/"

#package_path = "/Users/raikotto/Downloads/"
#package_path = "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/Uniquorn/"

#database_path = paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )

sim_list       = initiate_db_and_load_data( 
    ref_gen = ref_gen, 
    request_table = "sim_list" 
)
sim_list_stats = initiate_db_and_load_data(
    ref_gen = ref_gen, 
    request_table = "sim_list_stats"
)

sim_list = sim_list[ ! grepl("CUSTOM", sim_list$CL) ,]
sim_list_stats = sim_list_stats[ ! grepl("CUSTOM", sim_list_stats$CL) ,]

list_of_cls = as.character( sim_list_stats$CL )

return_chrom_start = function( vec, res ){
  
  vec = as.character( unlist( vec ) )
  
    alternative = paste0( 
      c( rep( "A", as.integer( vec[3] ) -  as.integer( vec[2] ) + 1 ) ), 
      collapse = ""
    )
  
  vec_res = matrix(
    ncol = length(labels), 
    c(
        vec[1],
        vec[2],
        ".",
        "C",
        alternative,
        ".",
        ".",
        ".",
        ".",
        "1|1"
    )
  )
  
  res = rbind( vec_res, res  )
}

labels = as.vector( c( "##fileformat=VCFv4.1\r\n##reference=GRCh37\r\n##contig=<ID=chrM,length=16571,assembly=hg19>\r\n##contig=<ID=chr1,length=249250621,assembly=hg19>\r\n##contig=<ID=chr2,length=243199373,assembly=hg19>\r\n##contig=<ID=chr3,length=198022430,assembly=hg19>\r\n##contig=<ID=chr4,length=191154276,assembly=hg19>\r\n##contig=<ID=chr5,length=180915260,assembly=hg19>\r\n##contig=<ID=chr6,length=171115067,assembly=hg19>\r\n##contig=<ID=chr7,length=159138663,assembly=hg19>\r\n##contig=<ID=chr8,length=146364022,assembly=hg19>\r\n##contig=<ID=chr9,length=141213431,assembly=hg19>\r\n##contig=<ID=chr10,length=135534747,assembly=hg19>\r\n##contig=<ID=chr11,length=135006516,assembly=hg19>\r\n##contig=<ID=chr12,length=133851895,assembly=hg19>\r\n##contig=<ID=chr13,length=115169878,assembly=hg19>\r\n##contig=<ID=chr14,length=107349540,assembly=hg19>\r\n##contig=<ID=chr15,length=102531392,assembly=hg19>\r\n##contig=<ID=chr16,length=90354753,assembly=hg19>\r\n##contig=<ID=chr17,length=81195210,assembly=hg19>\r\n##contig=<ID=chr18,length=78077248,assembly=hg19>\r\n##contig=<ID=chr19,length=59128983,assembly=hg19>\r\n##contig=<ID=chr20,length=63025520,assembly=hg19>\r\n##contig=<ID=chr21,length=48129895,assembly=hg19>\r\n##contig=<ID=chr22,length=51304566,assembly=hg19>\r\n##contig=<ID=chrX,length=155270560,assembly=hg19>\r\n##contig=<ID=chrY,length=59373566,assembly=hg19>\r\n#CHROM",
                       "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","MEMBER"))

for (CL in list_of_cls){
  
    CL_output = str_replace( CL,"/","")
    print (CL)
  
    output_file = paste0( 
    c( 
        output_directory, 
        paste( CL_output, "vcf", sep ="." ) 
    ), 
        collapse = "/"
    )
  
  #if(  ! file.exists(output_file)){
  
    res <-- matrix( ncol = length(labels), nrow = 0 )
    
    mapping = which( sim_list$CL %in% CL, arr.ind = T  )
    mutations = sim_list$Fingerprint[ mapping ]
    
    split_str = str_split( mutations, "_"  )
    
    res = sapply( split_str, FUN = return_chrom_start, res )
    res = t(res)
    colnames(res)  = labels
    write.table( x = res, file = output_file, sep ="\t", row.names =F, quote = F )
  #}
}
