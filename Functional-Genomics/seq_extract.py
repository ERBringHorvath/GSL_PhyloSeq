from Bio import SeqIO

##Function to extract sequencing matching search key
def extract_sequences_seqio(input_file, output_file, genome_name):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            print(f"Checking: {record.description}")
            if f'genome_name:{genome_name}' in record.description:
                found = True
                SeqIO.write(record, outfile, "fasta")
    
    return found

##Define input and output files
input_file = '/path/to/file'
output_file = '/path/to/file'
genome_name = 'search_key' ##Define seach key

if not extract_sequences_seqio(input_file, output_file, genome_name):
    print(f"No sequences found for genome name {genome_name}")
