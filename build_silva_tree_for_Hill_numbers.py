#!/usr/bin/env python3
"""
Build phylogenetic tree from SILVA database 138.2
This script extracts sequences based on accession numbers and builds a phylogenetic tree
"""

import os
import sys
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests

# Configuration
WORK_DIR = Path("/Users/paulaeterovick/silva_phylogeny")
SILVA_URL = "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
ACCESSION_FILE = Path("/Users/paulaeterovick/silva_accessions.txt")
OUTPUT_TREE_FILE = WORK_DIR / "phylogenetic_tree.tree"
OUTPUT_NEWICK_FILE = WORK_DIR / "phylogenetic_tree.newick"

def setup_workdir():
    """Create working directory"""
    WORK_DIR.mkdir(exist_ok=True)
    print(f"Working directory: {WORK_DIR}")

def read_accessions():
    """Read SILVA accession numbers from file"""
    print(f"Reading accession numbers from {ACCESSION_FILE}")
    with open(ACCESSION_FILE, 'r') as f:
        accessions = set(line.strip() for line in f if line.strip())
    print(f"Found {len(accessions)} accession numbers")
    return accessions

def download_silva_database():
    """Download SILVA database if not present"""
    silva_fasta_gz = WORK_DIR / "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
    silva_fasta = WORK_DIR / "SILVA_138.2_SSURef_NR99_tax_silva.fasta"
    
    if silva_fasta.exists():
        print(f"SILVA database already exists at {silva_fasta}")
        return silva_fasta
    
    if not silva_fasta_gz.exists():
        print(f"Downloading SILVA database from {SILVA_URL}")
        print("This may take a while (file is ~400MB)...")
        subprocess.run(["wget", "-O", str(silva_fasta_gz), SILVA_URL], check=True)
    
    print("Decompressing SILVA database...")
    subprocess.run(["gunzip", "-k", str(silva_fasta_gz)], check=True)
    
    return silva_fasta

def extract_sequences(silva_fasta, accessions):
    """Extract sequences from SILVA database matching accession numbers"""
    output_fasta = WORK_DIR / "silva_sequences.fasta"
    
    print(f"Extracting sequences from {silva_fasta}")
    print(f"Looking for {len(accessions)} sequences...")
    
    found_sequences = []
    found_accessions = set()
    sequence_name_counts = {}
    
    # Parse SILVA fasta file
    # SILVA format: >accession.version taxonomy
    for record in SeqIO.parse(silva_fasta, "fasta"):
        # Extract accession from header (format: >AB000278.1.1498 Bacteria;...)
        silva_id = record.id.split('.')[0]  # Get base accession without version
        
        if silva_id in accessions:
            found_accessions.add(silva_id)
            # Clean up the ID for tree building and handle duplicates
            if silva_id in sequence_name_counts:
                sequence_name_counts[silva_id] += 1
                record.id = f"{silva_id}_{sequence_name_counts[silva_id]}"
            else:
                sequence_name_counts[silva_id] = 1
                record.id = silva_id
            record.description = ""
            found_sequences.append(record)
            
            if len(found_sequences) % 100 == 0:
                print(f"  Found {len(found_sequences)} sequences...")
    
    print(f"\nFound {len(found_sequences)} out of {len(accessions)} requested sequences")
    
    # Report duplicates
    duplicates = {k: v for k, v in sequence_name_counts.items() if v > 1}
    if duplicates:
        print(f"Found {len(duplicates)} accessions with duplicates (renamed with suffixes)")
        for acc, count in sorted(duplicates.items()):
            print(f"  {acc}: {count} copies")
    
    # Report missing accessions
    missing = accessions - found_accessions
    if missing:
        print(f"Missing {len(missing)} accessions")
        missing_file = WORK_DIR / "missing_accessions.txt"
        with open(missing_file, 'w') as f:
            for acc in sorted(missing):
                f.write(f"{acc}\n")
        print(f"Missing accessions written to {missing_file}")
    
    # Write extracted sequences
    print(f"Writing extracted sequences to {output_fasta}")
    SeqIO.write(found_sequences, output_fasta, "fasta")
    
    return output_fasta, len(found_sequences)

def align_sequences(input_fasta):
    """Align sequences using MAFFT"""
    output_aligned = WORK_DIR / "silva_aligned.fasta"
    
    print(f"Aligning sequences with MAFFT...")
    print("This may take a while for large datasets...")
    
    # Use MAFFT auto mode
    cmd = ["mafft", "--auto", "--thread", "-1", str(input_fasta)]
    
    with open(output_aligned, 'w') as outf:
        subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, check=True)
    
    print(f"Aligned sequences written to {output_aligned}")
    return output_aligned
def build_tree(aligned_fasta):
    """Build phylogenetic tree using FastTree"""
    print(f"Building phylogenetic tree with FastTree...")
    print("Using GTR+CAT model for nucleotide sequences...")
    
    # FastTree for nucleotide sequences
    # Use full path to 'fasttree' executable if available to avoid command not found errors
    # Also capture stderr to analyze errors
    cmd = ["fasttree", "-gtr", "-nt", str(aligned_fasta)]
    
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"FastTree error output:\n{result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    
    with open(OUTPUT_TREE_FILE, 'w') as outf:
        outf.write(result.stdout)
    
    print(f"Tree built successfully: {OUTPUT_TREE_FILE}")
    
    # Create copy with .newick extension
    import shutil
    shutil.copy(OUTPUT_TREE_FILE, OUTPUT_NEWICK_FILE)
    print(f"Tree also saved as: {OUTPUT_NEWICK_FILE}")
    
    return OUTPUT_TREE_FILE

def check_dependencies():
    """Check if required tools are installed"""
    required = {
        'mafft': 'Multiple sequence alignment tool',
        'FastTree': 'Phylogenetic tree building tool (try: conda install -c bioconda fasttree)',
    }
    
    missing = []
    for tool, description in required.items():
        result = subprocess.run(['which', tool], capture_output=True)
        if result.returncode != 0:
            missing.append(f"{tool}: {description}")
    
    if missing:
        print("\nERROR: Missing required tools:")
        for item in missing:
            print(f"  - {item}")
        print("\nTo install with conda:")
        print("  conda install -c bioconda mafft fasttree")
        print("\nOr with homebrew:")
        print("  brew install mafft fasttree")
        sys.exit(1)

def main():
    print("=" * 80)
    print("SILVA Phylogenetic Tree Builder")
    print("=" * 80)
    
    # Check dependencies
    check_dependencies()
    
    # Setup
    setup_workdir()
    
    # Read accessions
    accessions = read_accessions()
    
    # Download SILVA database
    silva_fasta = download_silva_database()
    
    # Extract sequences
    extracted_fasta, num_seqs = extract_sequences(silva_fasta, accessions)
    
    if num_seqs == 0:
        print("\nERROR: No sequences found. Cannot build tree.")
        sys.exit(1)
    
    if num_seqs < 4:
        print(f"\nWARNING: Only {num_seqs} sequences found. Need at least 4 for meaningful tree.")
    
    # Align sequences
    aligned_fasta = align_sequences(extracted_fasta)
    
    # Build tree
    tree_file = build_tree(aligned_fasta)
    
    print("\n" + "=" * 80)
    print("SUCCESS!")
    print("=" * 80)
    print(f"Phylogenetic tree files created:")
    print(f"  - {OUTPUT_TREE_FILE}")
    print(f"  - {OUTPUT_NEWICK_FILE}")
    print(f"\nSequences: {num_seqs}")
    print(f"Working directory: {WORK_DIR}")
    print("=" * 80)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
