#!/usr/bin/env python3
"""Generate all possible DNA sequences coding for a given protein.

This module provides functionality to generate all possible DNA sequences that could
code for a given protein sequence, taking into account the redundancy of the genetic code.

Typical usage example:
    python ch1_1.py input.fasta -o output.fasta
"""

import argparse
import sys
from typing import List, Dict, Optional

from Bio.Data import CodonTable
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Constants
MAX_PROTEIN_LENGTH = 100  # Practical limit to avoid memory issues
CODON_TABLE = None  # Will be initialized in main()


def create_codon_table() -> Dict[str, List[str]]:
    """Create a dictionary mapping amino acids to their possible codons.

    Returns:
        Dict[str, List[str]]: Mapping of amino acids to their possible codons.
    """
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    codon_table = {}
    for codon, amino in standard_table.forward_table.items():
        if amino not in codon_table:
            codon_table[amino] = []
        codon_table[amino].append(codon)
    # Add stop codons
    codon_table['*'] = standard_table.stop_codons
    return codon_table


def save_dna_sequences(sequences: List[str], output_path: str) -> None:
    """Write DNA sequences to a FASTA file.

    Args:
        sequences: List of DNA sequences to write.
        output_path: Path to the output FASTA file.

    Raises:
        IOError: If there's an error writing to the output file.
    """
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            for idx, seq in enumerate(sequences, 1):
                f.write(f">variant_{idx}\n{seq}\n")
    except IOError as e:
        sys.exit(f"Error writing to output file: {e}")


def parse_protein_sequence(fasta_path: str) -> Optional[str]:
    """Read and validate a protein sequence from a FASTA file.

    Args:
        fasta_path: Path to the input FASTA file.

    Returns:
        The protein sequence if valid.

    Raises:
        ValueError: If the sequence contains invalid amino acids or is too long.
        IOError: If there's an error reading the input file.
    """
    try:
        with open(fasta_path, encoding='utf-8') as handle:
            for _, sequence in SimpleFastaParser(handle):
                protein_sequence = sequence.upper()
                
                if len(protein_sequence) > MAX_PROTEIN_LENGTH:
                    raise ValueError(
                        f"Protein sequence too long. Maximum length is {MAX_PROTEIN_LENGTH}")
                
                invalid_residues = set(protein_sequence) - VALID_AMINO_ACIDS
                if invalid_residues:
                    raise ValueError(
                        f"Invalid amino acid(s): {', '.join(invalid_residues)}")
                
                return protein_sequence
            raise ValueError("No sequences found in FASTA file")
    except IOError as e:
        sys.exit(f"Error reading input file: {e}")


def generate_coding_sequences(protein_sequence: str) -> List[str]:
    """Generate all possible DNA sequences for a protein sequence.

    Args:
        protein_sequence: The protein sequence to generate DNA sequences for.

    Returns:
        List of possible DNA sequences.
    """
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


def main() -> None:
    """Process protein sequence and generate coding DNA sequences."""
    global CODON_TABLE, VALID_AMINO_ACIDS
    
    try:
        parser = argparse.ArgumentParser(
            description='Generate all possible DNA sequences for a protein')
        parser.add_argument('fasta_file', help='Input protein FASTA file')
        parser.add_argument('-o', '--output', default='dna_variants.fasta',
                          help='Output FASTA file (default: dna_variants.fasta)')
        args = parser.parse_args()

        # Initialize global constants
        CODON_TABLE = create_codon_table()
        VALID_AMINO_ACIDS = set(CODON_TABLE.keys())

        protein_sequence = parse_protein_sequence(args.fasta_file)
        if not protein_sequence:
            sys.exit("No valid protein sequence found")

        dna_sequences = generate_coding_sequences(protein_sequence)
        save_dna_sequences(dna_sequences, args.output)

        print(f"Generated {len(dna_sequences)} possible DNA sequences")
        print(f"Results written to {args.output}")

    except Exception as e:
        sys.exit(f"Error: {str(e)}")


if __name__ == "__main__":
    main()
