#!/usr/bin/env python3
"""Generate all possible DNA sequences coding for a given protein."""

import argparse
from typing import List
from Bio.SeqIO.FastaIO import SimpleFastaParser

CODON_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M': ['ATG'],
    'N': ['AAT', 'AAC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA']
}

VALID_AMINO_ACIDS = set(CODON_TABLE.keys())


def parse_protein_sequence(fasta_path: str) -> str:
    """Read and validate a protein sequence from a FASTA file."""
    with open(fasta_path, encoding='utf-8') as handle:
        for _, sequence in SimpleFastaParser(handle):
            protein_sequence = sequence.upper()
            invalid_residues = set(protein_sequence) - VALID_AMINO_ACIDS
            if invalid_residues:
                raise ValueError(
                    f"Invalid amino acid(s): {', '.join(invalid_residues)}")
            return protein_sequence
    raise ValueError("No sequences found in FASTA file")


def generate_coding_sequences(protein_sequence: str) -> List[str]:
    """Generate all possible DNA sequences for a protein sequence."""
    def build_combinations(current: str, remaining: str, sequences: List[str]) -> None:
        if not remaining:
            sequences.append(current)
            return

        amino_acid = remaining[0]
        for codon in CODON_TABLE[amino_acid]:
            build_combinations(current + codon, remaining[1:], sequences)

    sequences: List[str] = []
    build_combinations("", protein_sequence, sequences)
    return sequences


def save_dna_sequences(sequences: List[str], output_path: str) -> None:
    """Write DNA sequences to a FASTA file."""
    with open(output_path, 'w', encoding='utf-8') as f:
        for idx, seq in enumerate(sequences, 1):
            f.write(f">variant_{idx}\n{seq}\n")


def main() -> None:
    """Process protein sequence and generate coding DNA sequences."""
    parser = argparse.ArgumentParser(
        description='Generate all possible DNA sequences for a protein')
    parser.add_argument('fasta_file', help='Input protein FASTA file')
    parser.add_argument('-o', '--output', default='dna_variants.fasta',
                        help='Output FASTA file (default: dna_variants.fasta)')
    args = parser.parse_args()

    protein_sequence = parse_protein_sequence(args.fasta_file)
    dna_sequences = generate_coding_sequences(protein_sequence)
    save_dna_sequences(dna_sequences, args.output)

    print(f"Generated {len(dna_sequences)} possible DNA sequences")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
