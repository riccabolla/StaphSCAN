import urllib.request
import gzip
import shutil
import subprocess
import pytest
import pandas as pd
from pathlib import Path 

# reference genome test1
@pytest.fixture(scope="module")
def s_aureus_fasta(tmp_path_factory):
    """Downloads and unzips the reference genome before the test runs."""
    tmp_dir = tmp_path_factory.mktemp("data")
    gz_path = tmp_dir / "GCF_000013425.1.fna.gz"
    fasta_path = tmp_dir / "GCF_000013425.1.fna"

    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz"
    urllib.request.urlretrieve(url, gz_path)

    with gzip.open(gz_path, 'rb') as f_in:
        with open(fasta_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return fasta_path

# mixed genome test2
@pytest.fixture(scope="module")
def local_fasta_mixed():
    """Locates fasta file to test"""
    current_dir = Path(__file__).parent
    fasta_mixed_path = current_dir / "data" / "test_mixed.fasta"
    
    assert fasta_mixed_path.exists(), f"Could not find local test file at {fasta_mixed_path}"
    return fasta_mixed_path

# pass genome test3
@pytest.fixture(scope="module")
def local_fasta_pass():
    """Locates fasta file to test"""
    current_dir = Path(__file__).parent
    fasta_pass_path = current_dir / "data" / "test_pass.fasta"
    
    assert fasta_pass_path.exists(), f"Could not find local test file at {fasta_pass_path}"
    return fasta_pass_path


# Reference genome
def test_reference_genome(s_aureus_fasta, tmp_path):
    """Runs test on the downloaded NCBI reference."""
    out_dir = tmp_path / "staphscan_output"
    
    cmd = ["staphscan", "-i", str(s_aureus_fasta), "-o", str(out_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"Pipeline failed! Error log:\n{result.stderr}"
    
    summary_file = out_dir / "staphscan_summary.tsv"
    assert summary_file.exists(), f"Expected output file not found at {summary_file}"
    
    df = pd.read_csv(summary_file, sep='\t')
    assert not df.empty, "staphscan_summary.tsv is empty"

    result_row = df.iloc[0]
    
    assert result_row['Species'] == "S. aureus"
    assert result_row['QC'] == "FAILED (Ambiguous_Bases)"
    assert 2600000 <= result_row['Total_size'] <= 3100000


# Contaminated genome
def test_mixed_genome(local_fasta_mixed, tmp_path):
    """Runs test on custom fasta"""
    out_dir = tmp_path / "staphscan_local_output"
    
    cmd = ["staphscan", "-i", str(local_fasta_mixed), "-o", str(out_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"Pipeline failed! Error log:\n{result.stderr}"
    
    summary_file = out_dir / "staphscan_summary.tsv"
    assert summary_file.exists(), f"Expected output file not found at {summary_file}"
    
    df = pd.read_csv(summary_file, sep='\t')
    assert not df.empty, "staphscan_summary.tsv is empty"

    result_row = df.iloc[0]
    
    assert result_row['Species'] == "S. aureus"  
    assert result_row['QC'] == "FAILED (Too long,Mixed)"                 
    assert result_row['Total_size'] > 5500000              
    assert float(result_row['Mash_distance']) >= 0.02 

# pass genome    
def test_pass_genome(local_fasta_pass, tmp_path):
    """Runs test on custom fasta"""
    out_dir = tmp_path / "staphscan_local_output"
    
    cmd = ["staphscan", "-i", str(local_fasta_pass), "-o", str(out_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"Pipeline failed! Error log:\n{result.stderr}"
    
    summary_file = out_dir / "staphscan_summary.tsv"
    assert summary_file.exists(), f"Expected output file not found at {summary_file}"
    
    df = pd.read_csv(summary_file, sep='\t')
    assert not df.empty, "staphscan_summary.tsv is empty"

    result_row = df.iloc[0]
    
    assert result_row['Species'] == "S. aureus"  
    assert result_row['QC'] == "PASS"                 
    assert result_row['Total_size'] == 2845438              
    assert float(result_row['Mash_distance']) == 0.00785531
    assert result_row["ST"] == "ST1"
    assert result_row["res_score"] == 2
    assert result_row["clfAB"] == "Complete"
