#!/bin/bash

cd /home/mehmet/Documents/projects/rnaseq
BAM_INPUT_DIR="bam_files"
MAIN_OUTPUT_DIR="stringtie_results"

mkdir -p $MAIN_OUTPUT_DIR

for file in $BAM_INPUT_DIR/Rep*D*.bam; do
    d_num=${file#*D}         # Remove everything before and including 'D'
    d_num=Day${d_num%%[^0-9]*} # Keep only the number immediately after 'D'

    mkdir -p "$BAM_INPUT_DIR/$d_num"  # Create directory if it doesn't exist
    mkdir -p "$MAIN_OUTPUT_DIR/$d_num" # Create dictroy for StringTie output
    mv "$file" "$BAM_INPUT_DIR/$d_num/"  # Move the file
    # add outher output direcetories like plots as well
done

echo "Files sorted into respective Day directories!"

###	This may become hand in the future if you want to make a similar program that will make KO vs WT comparasion
