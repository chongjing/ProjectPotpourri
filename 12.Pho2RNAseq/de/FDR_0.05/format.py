import csv
import sys

## this script is to add MSU_geneID and RAP_geneID from /rds-d7/user/cx264/hpc-work/project/0.ref/1.Rice_Nipponbare/NIP-T2T.gene_MSU_RAP.bed to edgeR results .csv file
## Run: python reformat.py /rds-d7/user/cx264/hpc-work/project/0.ref/1.Rice_Nipponbare/NIP-T2T.gene_MSU_RAP.bed 1.4.QLF.DE.results.WTLP-myc_vs_mock.FINAL.csv 1.5.QLF.DE.results.WTLP-myc_vs_mock.FINAL.csv
## Chongjing Xia, 20240303

# The script takes three arguments: the two input file names and the output file name
first_file = sys.argv[1]
second_file = sys.argv[2]
output_file = sys.argv[3]

# Read the first file and store the geneID and its corresponding data in a dictionary
with open(first_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # Skip the header row
    gene_dict = {rows[0]: rows[1:3] for rows in reader}  # Store geneID and its corresponding data

# Open the second file and the output file
with open(second_file, 'r') as f, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(f)
    writer = csv.writer(outfile)
    header = next(reader)  # Read the header row
    header.extend(['MSUloc', 'RAP_Id'])  # Extend the header row with new column names
    writer.writerow(header)  # Write the header row to the output file

    # For each row in the second file, if the geneID is in the dictionary, get the data and write to the output file
    for row in reader:
        if row[1] in gene_dict:
            row.extend(gene_dict[row[1]])
            writer.writerow(row)

