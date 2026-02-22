bacdata_dict = {
    "Escherichia_coli": "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
    "Bacillus_subtilis": "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
    "Staphylococcus_aureus": "ATATATAFCATATATCGATATATATATATATTATTATTATTA"
}


unknown_species = "ATATATAFCATATATCGATATATATATATATTATTATTATTA"


def sequence_identity(seq1, seq2):
    length = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches /length) * 100 if length > 0 else 0

best_match = None
highest_identity = 0
for species, sequence in bacdata_dict.items():
    identity = sequence_identity(unknown_species, sequence)
    if identity > highest_identity:
        highest_identity = identity
        best_match = species
print(f"The unknown species is most likely: {best_match} with {highest_identity:.2f}% identity.")