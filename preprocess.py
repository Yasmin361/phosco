from Bio import SeqIO

def split_multifasta(input_file, max_sequences=250):
    # Initialize variables
    chunk_count = 0
    sequence_chunk = []

    # Read the multifasta file
    with open(input_file, "r") as fasta_file:
        file_name = input_file.split(".fasta")[0].split("/")[-1]
        print(file_name)
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_chunk.append(record)
            if len(sequence_chunk) == max_sequences:
                chunk_count += 1
                output_file = f"preprocessed_full/{file_name}_chunk_{chunk_count}.fasta"
                with open(output_file, "w") as output:
                    SeqIO.write(sequence_chunk, output, "fasta")
                sequence_chunk = []
        # Write remaining sequences
        if sequence_chunk:
            chunk_count += 1
            output_file = f"preprocessed_full/{file_name}_chunk_{chunk_count}.fasta"
            with open(output_file, "w") as output:
                SeqIO.write(sequence_chunk, output, "fasta")


# Example usage
input_file = "input/uniprotkb_organism_id_9606_AND_keyword_2025_02_09.fasta"  # Replace with your actual file path
split_multifasta(input_file)
