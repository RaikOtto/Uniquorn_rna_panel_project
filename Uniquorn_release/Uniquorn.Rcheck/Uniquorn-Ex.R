pkgname <- "Uniquorn"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Uniquorn')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("add_custom_vcf_to_database")
### * add_custom_vcf_to_database

flush(stderr()); flush(stdout())

### Name: add_custom_vcf_to_database
### Title: Adds a custom vcf file to the three existing cancer cell line
###   panels
### Aliases: add_custom_vcf_to_database

### ** Examples

HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn");
add_custom_vcf_to_database( 
vcf_input_files = HT29_vcf_file,
library = "",
ref_gen = "GRCH37",
test_mode = TRUE,
n_threads = 1)



cleanEx()
nameEx("identify_vcf_file")
### * identify_vcf_file

flush(stderr()); flush(stdout())

### Name: identify_vcf_file
### Title: identify_VCF_file
### Aliases: identify_vcf_file

### ** Examples

HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn");

identification = identify_vcf_file( HT29_vcf_file )



cleanEx()
nameEx("initiate_canonical_databases")
### * initiate_canonical_databases

flush(stderr()); flush(stdout())

### Name: initiate_canonical_databases
### Title: initiate_canonical_databases
### Aliases: initiate_canonical_databases

### ** Examples

initiate_canonical_databases(
cosmic_file = "CosmicCLP_MutantExport.tsv",
ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
ref_gen = "GRCH37")



cleanEx()
nameEx("remove_custom_vcf_from_database")
### * remove_custom_vcf_from_database

flush(stderr()); flush(stdout())

### Name: remove_custom_vcf_from_database
### Title: Removes a cancer cell line training fingerprint (vcf file) from
###   the database. The names of all training sets can be seen by using the
###   function 'show_contained_cls'.
### Aliases: remove_custom_vcf_from_database

### ** Examples

remove_custom_vcf_from_database( 
name_cl = "HT29_CELLMINER", 
ref_gen = "GRCH37",
test_mode = TRUE )



cleanEx()
nameEx("show_contained_cls")
### * show_contained_cls

flush(stderr()); flush(stdout())

### Name: show_contained_cls
### Title: show_contained_cls
### Aliases: show_contained_cls

### ** Examples

contained_cls = show_contained_cls( 
ref_gen = "GRCH37")



cleanEx()
nameEx("show_contained_mutations")
### * show_contained_mutations

flush(stderr()); flush(stdout())

### Name: show_contained_mutations
### Title: show_contained_mutations
### Aliases: show_contained_mutations

### ** Examples

contained_cls = show_contained_mutations( ref_gen = "GRCH37" )



cleanEx()
nameEx("show_contained_mutations_for_cl")
### * show_contained_mutations_for_cl

flush(stderr()); flush(stdout())

### Name: show_contained_mutations_for_cl
### Title: show_contained_mutations_for_cl
### Aliases: show_contained_mutations_for_cl

### ** Examples

SK_OV_3_CELLMINER_mutations = show_contained_mutations_for_cl(
name_cl = "SK_OV_3_CELLMINER_mutations",
ref_gen = "GRCH37")



cleanEx()
nameEx("show_which_cls_contain_mutation")
### * show_which_cls_contain_mutation

flush(stderr()); flush(stdout())

### Name: show_which_cls_contain_mutation
### Title: show_which_cls_contain_mutation
### Aliases: show_which_cls_contain_mutation

### ** Examples

Cls_containing_mutations = show_which_cls_contain_mutation( 
mutation_name = "10_103354427_103354427", 
ref_gen = "GRCH37")



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
