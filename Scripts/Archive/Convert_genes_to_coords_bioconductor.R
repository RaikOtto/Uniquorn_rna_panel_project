library(stringr)
source("https://bioconductor.org/biocLite.R")
biocLite("Homo.sapiens")
library(Homo.sapiens)
genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)

t = read.table("~/Uniquorn_rna_panel_project/Misc/Hotspot_panel_v2_Themro_fisher_cancer_panel_49_genes.hg19.bed")
mycoords.list = as.character( t[,1] )

library("biomaRt")
listMarts(host="www.ensembl.org")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org")
results = getBM(
    values = mycoords.list,
    attributes = c(
      "hgnc_symbol",
      "ensembl_gene_id"
    ),
    filters = "hgnc_symbol",
    mart = ensembl
)

gb <- getBM(
  attributes=c(
    #"hgnc_symbol",
    "ensembl_gene_id",
    "chromosome_name",
    'ensembl_exon_id',
    "exon_chrom_start",
    "exon_chrom_end"
  ),
  filters = "ensembl_gene_id",
  values = results[,2],
  mart = ensembl
)
unique(gb[,1])

bed_file = cbind( 
    paste0( "chr",gb$chromosome_name ),
    gb$exon_chrom_start,
    gb$exon_chrom_end,
    gb$ensembl_gene_id
)
write.table(bed_file,"~/Uniquorn_rna_panel_project/Misc/Hotspot_panel_v2_Thermo_fisher_cancer_panel_49_genes.hg19.bed", sep ="\t", row.names = F, col.names = F, quote = F)
