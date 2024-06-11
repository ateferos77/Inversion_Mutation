

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import random


def apply_inversion(sequence, start, end):
    """Apply inversion mutation to a DNA sequence."""
    inverted_region = sequence[start:end].reverse_complement()
    sequence = sequence[:start] + inverted_region + sequence[end:]
    return sequence

def generate_inversions(genome_file, mutation_length, mutation_percentage, output_file):
    """Generate inversions in the genome sequence."""
    genome_record = SeqIO.read(genome_file, "fasta")
    genome_length = len(genome_record.seq)
    num_mutations = int(genome_length * mutation_percentage / 100)
    modified_regions = []
    
    for _ in range(num_mutations):
        attempt = 0
        max_attempts = 1000  # Set a limit to avoid infinite loops
        while attempt < max_attempts:
            start = random.randint(0, genome_length - mutation_length)
            end = start + mutation_length
            if not any(start < modified_end and end > modified_start for modified_start, modified_end in modified_regions):
                break
            attempt += 1
        else:
            print("Max attempts reached, cannot find non-overlapping region.")
            break
        
        genome_record.seq = apply_inversion(genome_record.seq, start, end)
        modified_regions.append((start, end))
        feature = SeqFeature(FeatureLocation(start, end), type="inversion")
        genome_record.features.append(feature)

    with open(output_file, "w") as output_handle:
        SeqIO.write(genome_record, output_handle, "fasta")

    with open(output_file + "_modifications.txt", "w") as modification_handle:
        modification_handle.write(f"Number of mutations applied: {num_mutations}\n")
        modification_handle.write("Modified regions:\n")
        for start, end in modified_regions:
            modification_handle.write(f"{start}-{end}\n")

if __name__ == "__main__":
    genome_file = input("Enter the path to the genome file: ")
    mutation_length = int(input("Enter the mutation length: "))
    mutation_percentage = float(input("Enter the mutation percentage: "))
    output_file = input("Enter the output file path: ")

    generate_inversions(genome_file, mutation_length, mutation_percentage, output_file)
