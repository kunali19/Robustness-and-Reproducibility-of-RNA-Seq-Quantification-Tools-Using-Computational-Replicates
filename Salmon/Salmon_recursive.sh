#!/bin/bash
# Define the path to the kallisto index
INDEX="/project/fmohebbi_1178/RNA-seq_project/Kunali_project/Salmon/salmon_index"
# Define the output directory
OUTPUT_DIR="/project/fmohebbi_1178/RNA-seq_project/Kunali_project/Salmon/Samples/Results"
# Define the list of file prefixes
FILES=("ERR188071" "ERR188218" "ERR188294" "ERR188423" "ERR188434" "ERR188444")
# Define the list of substitutions
SUBSTITUTIONS=("g" "rv" "s1" "s2" "s3")
# Iterate over each file prefix
for FILE_PREFIX in "${FILES[@]}"; do
    # Iterate over each substitution
    for SUBSTITUTION in "${SUBSTITUTIONS[@]}"; do
        # Run salmon quant with appropriate substitutions
        Salmon quant -i "$INDEX" \
                       -o "$OUTPUT_DIR/${FILE_PREFIX}_${SUBSTITUTION}" \
                       "${FILE_PREFIX}_${SUBSTITUTION}_1.fastq" \
                       "${FILE_PREFIX}_${SUBSTITUTION}_2.fastq"
    done
done