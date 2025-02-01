import os
import pytest
from ch1_1 import (
    parse_protein_sequence,
    generate_coding_sequences,
    save_dna_sequences,
    VALID_AMINO_ACIDS
)

@pytest.fixture
def temp_fasta(tmp_path):
    """Create a temporary FASTA file with a test sequence."""
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(">test_protein\nMA*\n")
    return fasta_path

@pytest.fixture
def temp_output(tmp_path):
    """Provide a temporary output file path."""
    return tmp_path / "output.fasta"

def test_parse_protein_sequence(temp_fasta):
    """Test reading and validating protein sequence."""
    sequence = parse_protein_sequence(temp_fasta)
    assert sequence == "MA*"
    assert all(aa in VALID_AMINO_ACIDS for aa in sequence)

def test_parse_invalid_sequence(tmp_path):
    """Test handling of invalid amino acids."""
    invalid_fasta = tmp_path / "invalid.fasta"
    with open(invalid_fasta, "w", encoding="utf-8") as f:
        f.write(">test_protein\nMAX\n")
    
    with pytest.raises(ValueError) as exc:
        parse_protein_sequence(invalid_fasta)
    assert "Invalid amino acid(s)" in str(exc.value)

def test_generate_coding_sequences():
    """Test DNA sequence generation."""
    # Test with Methionine (M) which has only one codon
    sequences = generate_coding_sequences("M")
    assert sequences == ["ATG"]

    # Test with Alanine (A) which has four codons
    sequences = generate_coding_sequences("A")
    assert sorted(sequences) == sorted(["GCT", "GCC", "GCA", "GCG"])

    # Test with sequence "MA" (should have 4 combinations)
    sequences = generate_coding_sequences("MA")
    assert len(sequences) == 4
    assert all(seq.startswith("ATG") for seq in sequences)

def test_save_dna_sequences(temp_output):
    """Test saving DNA sequences to FASTA file."""
    test_sequences = ["ATGGCT", "ATGGCC"]
    save_dna_sequences(test_sequences, temp_output)

    # Verify file contents
    with open(temp_output, encoding="utf-8") as f:
        content = f.read().splitlines()
    
    assert len(content) == 4  # 2 sequences * 2 lines each
    assert content[0] == ">variant_1"
    assert content[1] == "ATGGCT"
    assert content[2] == ">variant_2"
    assert content[3] == "ATGGCC"

def test_empty_fasta(tmp_path):
    """Test handling of empty FASTA file."""
    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.touch()
    
    with pytest.raises(ValueError) as exc:
        parse_protein_sequence(empty_fasta)
    assert "No sequences found" in str(exc.value)

def test_nonexistent_file():
    """Test handling of non-existent input file."""
    with pytest.raises(FileNotFoundError):
        parse_protein_sequence("nonexistent.fasta")
