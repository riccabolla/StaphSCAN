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
    # tmp_path is a built-in pytest fixture that creates a unique temp folder for the output
    out_dir = tmp_path / "staphscan_output"
    
    # 1. Execute the full pipeline exactly as a user would in the terminal
    cmd = ["staphscan", "-i", str(s_aureus_fasta), "-o", str(out_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # 2. Verify the command ran successfully (exit code 0)
    assert result.returncode == 0, f"Pipeline failed! Error log:\n{result.stderr}"
    
    # 3. Verify the pipeline generated the expected output file
    summary_file = out_dir / "staphscan_summary.tsv"
    assert summary_file.exists(), f"Expected output file not found at {summary_file}"
    
    # 4. Parse the output file and validate the biological results
    df = pd.read_csv(summary_file, sep='\t')
    
    assert not df.empty, "staphscan_summary.tsv is empty"
    
    # Extract the first row (since we only passed one input)
    result_row = df.iloc[0]
    
    # Assert pipeline-level expectations
    assert result_row['Species'] == "S. aureus"
    assert result_row['QC'] == "FAILED (Ambiguous_Bases)"
    assert 2600000 <= result_row['Total_size'] <= 3100000