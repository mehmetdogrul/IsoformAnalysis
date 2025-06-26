#!/bin/bash

cd /home/mehmet/Documents/projects/rnaseq
BAM_INPUT_DIR="bam_files"
MAIN_OUTPUT_DIR="stringtie_results"

mkdir -p $MAIN_OUTPUT_DIR

for file in $BAM_INPUT_DIR/Rep*D*.bam; do
    d_num=${file#*D}         # Remove everything before and including 'D'
    d_num=Day${d_num%%[^0-9]*} # Keep only the number immediately after 'D'

    mkdir -p "$BAM_INPUT_DIR/$d_num"  # Create directory if it doesn't exist
    mv "$file" "$BAM_INPUT_DIR/$d_num/"  # Move the file

    DAY_DIR=($BAM_INPUT_DIR/Day*/)
    DIR_COUNT={#DAY_DIR[@]}

	### FIX THIS PART. THIS IS BROKEN. WRITE A FOR LOOP. YOU HAVE A TEMPLETE IN YOUR TABLET
    while [ $DIR_COUNT -ge 2 ]; do
        echo "Creating ${DAY_DIR[0]} vs ${DAY_DIR[1]} output directory"

        DAY_PAIR="${DAY_DIR[DIR_COUNT-2]} vs ${DAY_DIR[DIR_COUNT-1]}"

	mkdir -p "$MAIN_OUTPUT_DIR/$DAY_PAIR"

	((DIR_COUNT--))

    done

done

echo "Files sorted into respective Day directories!"

