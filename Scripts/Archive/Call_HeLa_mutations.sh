export PERL5LIB=/gpfs01/sw/prod_23082013/software_depot/x86_64/RedHat/6.3/vcftools/0.1.11/perl/
#VCF_TOOLS='/gpfs01/sw/prod_23082013/software_depot/x86_64/RedHat/6.3/samtools/0.1.19/bcftools/vcfutils.pl'

PROGRAM_FOLDER="/Users/raikotto/Downloads/"
FWD_READS="/Users/raikotto/Downloads/SRR1232554_1.fastq"
BKWD_READS="/Users/raikotto/Downloads/SRR1232554_2.fastq"
OUTPUT_FILE="/Users/raikotto/Downloads/Hela_rna_seq.sam"

# just check to parse the input right by deleting the last slash

#OUTPUT_FILE=$OUTPUT_FOLDER_SAMPLE'/'$SAMPLE'_celline_ident_tmp'
#RESULTS_FILE=$OUTPUT_FOLDER_SAMPLE'/'$SAMPLE'_celline_ident.txt'
REF_GEN_HG19_BOWTIE2='/Users/raikotto/Downloads/ref_gen/hg19/full_genome/' # change reference genome here
REF_GEN_HG19_FASTA='/Users/raikotto/Downloads/ref_gen/hg19/full_genome/hg19.fa'

# map the reads and create the pileup # mapping step

bowtie2 -p 8  --local -x $REF_GEN_HG19_BOWTIE2 -1 $FWD_READS -2 $BKW_READS -S $OUTPUT_FILE

# variant calling samtools

samtools view -bS $OUTPUT_FILE'.sam' > $OUTPUT_FILE'.bam'; # rewrite readable to compressed format
samtools sort $OUTPUT_FILE'.bam' $OUTPUT_FILE'_sorted'; # sort the file, note that .bam is added on default
samtools index $OUTPUT_FILE'_sorted.bam' # create lookup index
samtools mpileup -uf $REF_GEN_HG19_FASTA $OUTPUT_FILE'_sorted.bam' > $OUTPUT_FILE'.mpileup'; 

#cat $OUTPUT_FILE'.mpileup' | bcftools view -vcg - > $OUTPUT_FILE'_samtools.vcf'
#grep -v '##' $OUTPUT_FILE'_samtools.vcf' > $OUTPUT_FILE'_samtools.vcf_tmp'
#mv -f $OUTPUT_FILE'_samtools.vcf_tmp' $OUTPUT_FILE'_samtools.vcf'

picard AddOrReplaceReadGroups.jar I=$inputbam O=$bam RGLB="library" RGPL="DNA-seq" RGPU="" RGSM="HELA"

picard AddOrReplaceReadGroups RGSM=1804 RGLB=library RGPL=DNA-seq RGPU=default INPUT=Hela_dna_seq_sorted.bam OUTPUT=Hela_dna_seq_sorted.bam

# clean up
#rm $OUTPUT_FILE'.sam' $OUTPUT_FILE'.bam' $OUTPUT_FILE'_sorted.bam' $OUTPUT_FILE'_sorted.bam.bai' $OUTPUT_FILE'.mpileup' $OUTPUT_FILE'_samtools.vcf'

### GATK

picard AddOrReplaceReadGroups RGSM= default RGLB= library RGPL= DNA-seq RGPU= default INPUT= Hela_dna_seq_sorted.bam OUTPUT= output.bam CREATE_INDEX= True SORT_ORDER= coordinate VALIDATION_STRINGENCY= LENIENT
java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R ucsc.hg19.fasta -I Ovcar8.sorted.annot.bam -o gatk_test.vcf
# alias picard='/usr/libexec/java_home --exec java -jar ~/Downloads/picard-tools-2.0.1/picard.jar '
samtools faidx hg19.fa
picard CreateSequenceDictionary R= hg19.fa O= hg19.dict
