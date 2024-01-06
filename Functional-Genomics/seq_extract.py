from Bio import SeqIO

def extract_sequences_seqio(input_file, output_file, genome_name):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            print(f"Checking: {record.description}")
            if f'genome_name:{genome_name}' in record.description:
                found = True
                SeqIO.write(record, outfile, "fasta")
    
    return found

input_file = '/Users/ebh/Desktop/GO/Filtered_019_MLST_Hits.fa'
output_file = '/Users/ebh/Desktop/GO/GSL17_019_Hits.fa'
genome_name = 'GSL17_019'

if not extract_sequences_seqio(input_file, output_file, genome_name):
    print(f"No sequences found for genome name {genome_name}")