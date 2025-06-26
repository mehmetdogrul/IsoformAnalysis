#!/bin/bash

echo "sample_sorter.sh is running..."
source sample_sorter.sh

cd /home/mehmet/Documents/projects/rnaseq # this is already the pwd

for directory in $BAM_INPUT_DIR/Day*;do
    DAY=$(basename $directory)
    echo "StringTie is processing $DAY files."
    cd scripts
    source stringtie_pipeline.sh $DAY
done

# add r script here.

