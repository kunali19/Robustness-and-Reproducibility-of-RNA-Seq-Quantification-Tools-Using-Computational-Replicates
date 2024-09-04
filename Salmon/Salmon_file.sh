import os
import pandas as pd
# List of file suffixes
file_suffixes = ['_g', '_rv', '_v1', '_v2', '_v3']
# Loop over each file suffix
for suffix in file_suffixes:
    # Load the abundance file
    file_path = "/project/fmohebbi_1178/RNA-seq_project/Kunali_project/Salmon/New_Samples/Results/SRR896807" + suffix + "/quant.sf"
    df = pd.read_csv(file_path, sep='\t')
    # Calculate total mapped reads
    total_mapped_reads = 1000000
    # Convert TPM to gene counts
    df['gene_counts'] = (df['TPM'] * total_mapped_reads) / 10**6
    # Save the modified file
    output_file_path = "/project/fmohebbi_1178/RNA-seq_project/Kunali_project/Salmon/New_Samples/Results/SRR896807" + suffix + "/abundance_gene_counts.tsv"
    df.to_csv(output_file_path, sep='\t', index=False)
