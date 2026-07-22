import urllib.request
import gzip
import shutil
import subprocess
import pytest
import pandas as pd
from pathlib import Path

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

def test_full_staphscan_pipeline(s_aureus_fasta, tmp_path):
    """Runs the full staphscan CLI tool and verifies the TSV output."""
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