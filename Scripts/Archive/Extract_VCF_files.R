library("stringr")
library(devtools)
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

output_directory = "~/Uniquorn_data/benchmark_vcf_files/raw_files_new/"

database = show_contained_cls()

ccl_list = database$CCL

# header 

head = "##fileformat=VCFv4.1\r\n##reference=GRCh37\r\n##contig=<ID=chrM,length=16571,assembly=hg19>\r\n##contig=<ID=chr1,length=249250621,assembly=hg19>\r\n##contig=<ID=chr2,length=243199373,assembly=hg19>\r\n##contig=<ID=chr3,length=198022430,assembly=hg19>\r\n##contig=<ID=chr4,length=191154276,assembly=hg19>\r\n##contig=<ID=chr5,length=180915260,assembly=hg19>\r\n##contig=<ID=chr6,length=171115067,assembly=hg19>\r\n##contig=<ID=chr7,length=159138663,assembly=hg19>\r\n##contig=<ID=chr8,length=146364022,assembly=hg19>\r\n##contig=<ID=chr9,length=141213431,assembly=hg19>\r\n##contig=<ID=chr10,length=135534747,assembly=hg19>\r\n##contig=<ID=chr11,length=135006516,assembly=hg19>\r\n##contig=<ID=chr12,length=133851895,assembly=hg19>\r\n##contig=<ID=chr13,length=115169878,assembly=hg19>\r\n##contig=<ID=chr14,length=107349540,assembly=hg19>\r\n##contig=<ID=chr15,length=102531392,assembly=hg19>\r\n##contig=<ID=chr16,length=90354753,assembly=hg19>\r\n##contig=<ID=chr17,length=81195210,assembly=hg19>\r\n##contig=<ID=chr18,length=78077248,assembly=hg19>\r\n##contig=<ID=chr19,length=59128983,assembly=hg19>\r\n##contig=<ID=chr20,length=63025520,assembly=hg19>\r\n##contig=<ID=chr21,length=48129895,assembly=hg19>\r\n##contig=<ID=chr22,length=51304566,assembly=hg19>\r\n##contig=<ID=chrX,length=155270560,assembly=hg19>\r\n##contig=<ID=chrY,length=59373566,assembly=hg19>\r\n#"
header_row =  c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","MEMBER")

###

library_names = Uniquorn:::read_library_names(ref_gen = "GRCH37")

for (library_name in library_names){
    
    print(library_name)
    ccl_list = Uniquorn::show_contained_ccls(ref_gen = "GRCH37")

    g_mat = read_mutation_grange_objects(
        library_name = library_name,
        ref_gen = ref_gen,
        mutational_weight_inclusion_threshold = 0
    )
    
    member_ccls = as.character(ccl_list$CCL[ccl_list$Library == library_name])
    member_info = as.character( g_mat$Member_CCLs )
    member_info = as.character(sapply(member_info, FUN = function(vec){return(str_replace_all(vec,pattern = "\\[|\\]",""))}))
    
    for (ccl_name in member_ccls){
      
        print(ccl_name)
        ccl_name = as.character( str_replace_all(ccl_name,pattern = "\\[|\\]","") )
        
        vcf_file_name = paste0(
          c("~/Uniquorn_data/benchmark_vcf_files/raw_files_new/", ccl_name, ".",library_name,".vcf"),
          collapse= ""
        )
        
        if (file.exists(vcf_file_name))
            next()
      
        index = str_detect(member_info, pattern = ccl_name)
        var_info = g_mat[index,]
        as.data.frame(var_info)
        
        gr = g_mat[index,]
        
        chroms = paste0( "chr", as.character( seqnames(var_info) ) )
        positions = as.data.frame(ranges( gr))
        starts = positions$start
        length = positions$width
        id = rep(".", length(chroms))
        ref = rep("C", length(chroms))
        member_col = rep("1|1", length(chroms))
        vars   = sapply( length, FUN = function(vec){return(paste(rep("A", vec),collapse= "", sep =""))})
        
        vcf_file = matrix( cbind(
            chroms,
            starts,
            id,
            ref ,
            vars,
            id,
            id,
            id,
            id,
            member_col
        ) , nrow = length(chroms))
        vcf_file = matrix(rbind(header_row, vcf_file), nrow = length(chroms) + 1)
        vcf_file[1,1] = paste(head, vcf_file[1,1], sep ="")
        
        write.table(vcf_file, file = vcf_file_name, sep = "\t", quote =F, row.names = F, col.names = F)
        
    }
}
