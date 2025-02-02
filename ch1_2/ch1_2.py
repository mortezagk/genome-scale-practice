#!/usr/bin/env python3
"""Generate DNA sequences from protein sequences using codon usage frequencies."""

import random
import sys
import logging
from pathlib import Path
from typing import Dict, DefaultDict
from collections import defaultdict

from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Configure logging only for errors and debugging
logging.basicConfig(
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

CodonFrequencies = Dict[str, float]
CodonProbabilities = Dict[str, Dict[str, float]]


def get_codon_tables() -> Dict[str, str]:
    """Get standard genetic code tables in DNA format.
    
    Returns:
        Dict mapping codons to amino acids
    """
    standard_table = CodonTable.standard_dna_table
    dna_code = {**standard_table.forward_table}
    for codon in standard_table.stop_codons:
        dna_code[codon] = '*'
    return dna_code


def normalize_codon_frequencies(codon_freqs: CodonFrequencies) -> CodonProbabilities:
    """Normalize codon frequencies into probabilities for each amino acid.
    
    Args:
        codon_freqs: Dictionary of codon frequencies
        
    Returns:
        Nested dictionary mapping amino acids to their codon probabilities
    """
    aa_codons: DefaultDict[str, Dict[str, float]] = defaultdict(dict)
    for codon, freq in codon_freqs.items():
        aa = DNA_CODONS[codon]
        aa_codons[aa][codon] = freq
    
    return {
        aa: {codon: freq/sum(codons.values()) 
             for codon, freq in codons.items()}
        for aa, codons in aa_codons.items()
    }


def generate_dna_sequence(protein_seq: str, 
                         codon_probs: CodonProbabilities) -> str:
    """Generate random DNA sequence based on protein sequence and codon probabilities."""
    # Check for invalid amino acids first
    invalid_aas = [aa for aa in protein_seq if aa not in codon_probs]
    if invalid_aas:
        raise ValueError(f"Invalid amino acid(s) in sequence: {', '.join(invalid_aas)}")
    
    return ''.join(
        random.choices(
            list(codon_probs[aa].keys()),
            weights=list(codon_probs[aa].values())
        )[0]
        for aa in protein_seq
    )


def create_complete_codon_freqs(min_freq: int = 50, 
                              max_freq: int = 1000) -> CodonFrequencies:
    """Create a dictionary with frequencies for standard genetic code codons."""
    return {codon: random.randint(min_freq, max_freq) 
            for codon in DNA_CODONS}


def read_protein_sequence(fasta_path: Path) -> SeqRecord:
    """Read protein sequence from FASTA file."""
    try:
        return next(SeqIO.parse(fasta_path, "fasta"))
    except StopIteration:
        raise ValueError("No sequences found in FASTA file")
    except Exception as e:
        raise ValueError(f"Error reading FASTA file: {e}")


def main() -> None:
    """Main execution function."""
    if len(sys.argv) != 2:
        print("Usage: python ch1_2.py <protein_fasta_file>")
        sys.exit(1)

    try:
        fasta_path = Path(sys.argv[1])
        protein_record = read_protein_sequence(fasta_path)
        protein_sequence = str(protein_record.seq)

        example_freqs = create_complete_codon_freqs()
        codon_probs = normalize_codon_frequencies(example_freqs)
        dna_sequence = generate_dna_sequence(protein_sequence, codon_probs)

        print("\nInput protein sequence:")
        print(f"Description: {protein_record.description}")
        print(f"Length: {len(protein_sequence)} amino acids")
        print(f"Sequence: {protein_sequence}")
        
        print("\nGenerated DNA sequence:")
        print(f"Length: {len(dna_sequence)} nucleotides")
        print(f"Sequence: {dna_sequence}")

    except Exception as e:
        logger.error("Error: %s", str(e))
        sys.exit(1)


# Global constants
DNA_CODONS = get_codon_tables()


if __name__ == "__main__":
    main()