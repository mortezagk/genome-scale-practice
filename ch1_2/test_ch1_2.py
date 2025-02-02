import pytest
from pathlib import Path
from unittest.mock import patch
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from ch1_2 import (
    get_codon_tables,
    normalize_codon_frequencies,
    generate_dna_sequence,
    create_complete_codon_freqs,
    read_protein_sequence,
    DNA_CODONS
)

@pytest.fixture
def sample_protein_fasta(tmp_path):
    """Create a sample protein FASTA file."""
    fasta_path = tmp_path / "test_protein.fasta"
    with open(fasta_path, "w") as f:
        f.write(">test_protein\nMA*\n")
    return fasta_path

@pytest.fixture
def mock_codon_freqs():
    """Create sample codon frequencies."""
    return {
        'ATG': 100,  # M
        'GCT': 200,  # A
        'GCC': 100,  # A
        'TAA': 50,   # *
        'TAG': 50,   # *
    }

def test_get_codon_tables():
    """Test creation of codon translation table."""
    codons = get_codon_tables()
    assert 'ATG' in codons
    assert codons['ATG'] == 'M'
    assert 'TAA' in codons
    assert codons['TAA'] == '*'
    assert len(codons) == 64  # Standard genetic code size

def test_normalize_codon_frequencies(mock_codon_freqs):
    """Test normalization of codon frequencies to probabilities."""
    probs = normalize_codon_frequencies(mock_codon_freqs)
    
    # Check Methionine probabilities
    assert 'M' in probs
    assert probs['M']['ATG'] == 1.0
    
    # Check Alanine probabilities
    assert 'A' in probs
    assert abs(probs['A']['GCT'] - 0.666667) < 0.0001
    assert abs(probs['A']['GCC'] - 0.333333) < 0.0001
    
    # Check stop codon probabilities
    assert '*' in probs
    assert probs['*']['TAA'] == 0.5
    assert probs['*']['TAG'] == 0.5

def test_generate_dna_sequence():
    """Test DNA sequence generation with fixed probabilities."""
    protein_seq = "MA*"
    probs = {
        'M': {'ATG': 1.0},
        'A': {'GCT': 1.0},
        '*': {'TAA': 1.0}
    }
    
    dna_seq = generate_dna_sequence(protein_seq, probs)
    assert dna_seq == "ATGGCTTAA"

def test_generate_dna_sequence_invalid_aa():
    """Test DNA sequence generation with invalid amino acid."""
    protein_seq = "MXA"  # X is invalid
    probs = {'M': {'ATG': 1.0}, 'A': {'GCT': 1.0}}
    
    with pytest.raises(ValueError) as exc:
        generate_dna_sequence(protein_seq, probs)
    assert "Invalid amino acid(s)" in str(exc.value)

@patch('random.randint')
def test_create_complete_codon_freqs(mock_randint):
    """Test creation of complete codon frequencies."""
    mock_randint.return_value = 100
    freqs = create_complete_codon_freqs()
    
    assert len(freqs) == 64
    assert all(freq == 100 for freq in freqs.values())
    assert all(codon in DNA_CODONS for codon in freqs)

def test_read_protein_sequence(sample_protein_fasta):
    """Test reading protein sequence from FASTA."""
    record = read_protein_sequence(sample_protein_fasta)
    
    assert isinstance(record, SeqRecord)
    assert str(record.seq) == "MA*"
    assert record.id == "test_protein"

def test_read_protein_sequence_empty(tmp_path):
    """Test reading from empty FASTA file."""
    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.touch()
    
    with pytest.raises(ValueError) as exc:
        read_protein_sequence(empty_fasta)
    assert "No sequences found" in str(exc.value)

def test_read_protein_sequence_nonexistent():
    """Test reading from non-existent file."""
    with pytest.raises(ValueError) as exc:
        read_protein_sequence(Path("nonexistent.fasta"))
    assert "Error reading FASTA file" in str(exc.value)

@patch('random.choices')
def test_generate_dna_sequence_random(mock_choices):
    """Test DNA sequence generation with mocked random choices."""
    mock_choices.side_effect = lambda k, weights: [list(k)[0]]
    
    protein_seq = "MA"
    probs = {
        'M': {'ATG': 1.0},
        'A': {'GCT': 0.5, 'GCC': 0.5}
    }
    
    dna_seq = generate_dna_sequence(protein_seq, probs)
    assert dna_seq == "ATGGCT"
    assert mock_choices.call_count == 2
