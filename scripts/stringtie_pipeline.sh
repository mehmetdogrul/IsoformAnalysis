#!/bin/bash

cd /home/mehmet/Documents/projects/rnaseq
# Set variables
THREADS=8
REF_GTF="annotations/reference.gtf"
MERGED_GTF="merged.gtf"
MERGE_LIST="mergelist.txt"
OUTPUT_DIR="stringtie_results/$1"
BAM_DIR="bam_files/$1"

# Step 1: Run StringTie for each BAM file
echo "Running StringTie for each BAM file..."
for BAM in $BAM_DIR/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    mkdir -p "$OUTPUT_DIR/$SAMPLE"
    stringtie "$BAM" -p $THREADS -G $REF_GTF -o "$OUTPUT_DIR/$SAMPLE/${SAMPLE}.gtf"
done

# Step 2: Merge transcripts
echo "Merging transcripts..."
ls $OUTPUT_DIR/*/*.gtf > $OUTPUT_DIR/$MERGE_LIST
stringtie --merge -p $THREADS -G $REF_GTF -o $OUTPUT_DIR/$MERGED_GTF $OUTPUT_DIR/$MERGE_LIST

# Step 3: Re-estimate transcript abundances
echo "Re-estimating transcript abundances in .ctab format"
for BAM in $BAM_DIR/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    stringtie "$BAM" -p $THREADS -G $OUTPUT_DIR/$MERGED_GTF -o "$OUTPUT_DIR/$SAMPLE/${SAMPLE}_final.gtf" -e -B
done

# Step 4: Generate sample information file
echo "Creating sample information file..."
echo "sampleID,condition,day" > $OUTPUT_DIR/sample_info.csv
for BAM in $BAM_DIR/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    # Assign condition based on filename
    if [[ "$SAMPLE" == *KO* ]]; then
        CONDITION="KO"
    elif [[ "$SAMPLE" == *WT* ]]; then
        CONDITION="WT"
    else
        CONDITION="Unknown"
    fi
    echo "$SAMPLE,$CONDITION,$1" >> $OUTPUT_DIR/sample_info.csv
done

echo "Pipeline completed! Processed data is in the '$OUTPUT_DIR' directory."

