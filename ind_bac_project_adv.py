from Bio.Align import PairwiseAligner
from Bio import SeqIO


def get_bacdata_dict(file_path):
    # Dictionary comprehension to map IDs to sequence strings
    bacdata_dict = {
        record.description: str(record.seq) 
        for record in SeqIO.parse(file_path, "fasta")
    }
    return bacdata_dict

def read_single_sequence(file_path):
    """
    Reads a FASTA file containing exactly one record 
    and returns the sequence as a string.
    """
    try:
        record = SeqIO.read(file_path, "fasta")
        return str(record.seq)
    except ValueError:
        # This triggers if the file has 0 or >1 records
        print("Error: File must contain exactly one sequence.")
        return None
    
# Usage
bacdata_dict = get_bacdata_dict("seq.fasta")
for species, sequence in bacdata_dict.items():
    print(f"Species: {species}")
    print(f"Sequence: {sequence}")
    print("-" * 20) # Separator for readability


unknown_species = read_single_sequence("unknown.fasta")
if sequence:
    print(f"Sequence Length: {len(sequence)}")
    print(f"Sequence: {sequence}")

aligner = PairwiseAligner()
aligner.mode = 'global'

aligner.match_score = 1 
aligner.mismatch_score = 0 
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5 

best_match = None
highest_identity = 0

for species, sequence in bacdata_dict.items():
    alignments = aligner.align(unknown_species, sequence)
    best_alignment = alignments[0]
   
    matches = sum(
        1 for a, b in zip(best_alignment.target, best_alignment.query)
        if a == b
    )
    identity = (matches / len(best_alignment.target)) * 100
    print(f"{species}: {identity:.2f}% identity")

    if identity > highest_identity:
        highest_identity = identity
        best_match = species
print(f"\nUnknown organism is most similar to {best_match} ({highest_identity:.2f}% identity)")