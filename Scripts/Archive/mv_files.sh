#!/bin/bash


for FILE_NAME in /Users/raikotto/Uniquorn_data/vcfHG19/*.hg19.vcf_uniquorn_ident.tab; do

    echo $FILE_NAME

    NEW_FILE_NAME=${FILE_NAME//.hg19.vcf_uniquorn_ident.tab/.vcf_uniquorn_ident.tab}
    echo $NEW_FILE_NAME
    
    mv $FILE_NAME $NEW_FILE_NAME

done