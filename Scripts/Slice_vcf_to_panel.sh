PROGRAM_FOLDER=$HOME"/Uniquorn_data/benchmark_vcf_files/raw_files/"
PANEL_BED_FILE=$HOME"/Downloads/trusight_cancer_manifest_a.bed"
OUTPUT_PATH=$HOME"/Uniquorn_data/benchmark_vcf_files/panel_raw_files/"

INPUT_FILES=`find $PROGRAM_FOLDER -name "*.vcf" | egrep -e "*.vcf$" | egrep -e 'CCLE' `

for VCF_FILE in $INPUT_FILES; do

    OUT_FILE_NAME=`basename $VCF_FILE`
    OUTPUT_FILE=$OUTPUT_PATH"/"$OUT_FILE_NAME
    echo $OUTPUT_FILE
    
    #egrep -e "#" $VCF_FILE > $OUTPUT_FILE"_TMP"
    #cat $VCF_FILE | egrep -v '#' | awk '{print "chr"$0}' >> $OUTPUT_FILE"_TMP"
    #mv $OUTPUT_FILE"_TMP" $VCF_FILE
    
    bedtools intersect -a $VCF_FILE -b $PANEL_BED_FILE >> $OUTPUT_FILE
done
