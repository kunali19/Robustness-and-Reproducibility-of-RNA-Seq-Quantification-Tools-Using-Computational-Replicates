import pandas as pd
# Load the abundance file
file_path = "/project/fmohebbi_1178/RNA-seq_project/Kunali_project/Kallisto/Samples/Results/ERR188071_g/abundance.tsv"
df = pd.read_csv(file_path, sep='\t')
# Calculate total mapped reads
total_mapped_reads = df['est_counts'].sum()
# Convert TPM to gene counts
df['gene_counts'] = (df['tpm'] * total_mapped_reads) / 10**6
# Save the modified file
output_file_path = "/project/fmohebbi_1178/RNA-seq_project/Kunali_project/Kallisto/Samples/Results/ERR188071_g/abundance_gene_counts.tsv"
df.to_csv(output_file_path, sep='\t', index=False)
