from Bio import SeqIO

def filter_sequences_by_length(input_file, output_file, min_length):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequences = SeqIO.parse(infile, "fasta")
        filtered_sequences = (record for record in sequences if len(record.seq) >= min_length)

        count = SeqIO.write(filtered_sequences, outfile, "fasta")

    return count

input_file = '/Users/ebh/Desktop/GO/GSL17-019.fa'  # Replace with the path to your input FASTA file
output_file = '/Users/ebh/Desktop/GO/Filtered_GSL17_019.fa'  # Replace with the path to your desired output file
min_length = 200  # Minimum length requirement

num_sequences_written = filter_sequences_by_length(input_file, output_file, min_length)

print(f"Number of sequences longer than {min_length} bp written to {output_file}: {num_sequences_written}")
