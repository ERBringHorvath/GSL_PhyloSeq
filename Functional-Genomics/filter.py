from Bio import SeqIO

##Function to filter sequences
def filter_sequences_by_length(input_file, output_file, min_length):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequences = SeqIO.parse(infile, "fasta")
        filtered_sequences = (record for record in sequences if len(record.seq) >= min_length)

        count = SeqIO.write(filtered_sequences, outfile, "fasta")

    return count

##Load files and define output file
input_file = '/path/to/file'  
output_file = '/path/to/file'  
min_length = 40  # Minimum length requirement

num_sequences_written = filter_sequences_by_length(input_file, output_file, min_length)

print(f"Number of sequences longer than {min_length} bp written to {output_file}: {num_sequences_written}")
