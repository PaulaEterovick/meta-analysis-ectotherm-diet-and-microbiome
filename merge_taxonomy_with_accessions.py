#!/usr/bin/env python3
"""
Merge QIIME2 Taxonomy with SILVA Accession Numbers

This script takes the taxonomy.tsv and search-results.tsv files produced by
classify_with_silva_accessions.sh and merges them into a comprehensive 
mapping file that contains:
- Feature ID (OTU)
- Full Taxonomy
- Confidence Score
- SILVA Accession Number
- Percent Identity
- E-value
- Bit-score

Usage:
    python merge_taxonomy_with_accessions.py

Output:
    silva_taxonomy_mapping.tsv - Comprehensive mapping file
"""

import pandas as pd
import sys
import os
from pathlib import Path

# File paths - modify these if your files are in different locations
TAXONOMY_FILE = "taxonomyV3_V4_444_with_accessions_export/taxonomy.tsv"
BLAST_RESULTS_FILE = "silva_blast_results_export/blast6.tsv"
OUTPUT_FILE = "silva_taxonomy_mapping.tsv"

def load_taxonomy(file_path):
    """Load taxonomy data from QIIME2 export"""
    print(f"Loading taxonomy from {file_path}...")
    
    try:
        # Read taxonomy data
        taxonomy_df = pd.read_csv(file_path, sep='\t')
        
        # Check if there's a header in the expected format
        if not all(col in taxonomy_df.columns for col in ['Feature ID', 'Taxon', 'Confidence']):
            # If not, try again with no header
            taxonomy_df = pd.read_csv(file_path, sep='\t', 
                                      names=['Feature ID', 'Taxon', 'Confidence'])
            
        print(f"Loaded {len(taxonomy_df)} taxonomy records")
        return taxonomy_df
    
    except Exception as e:
        print(f"Error loading taxonomy: {e}")
        sys.exit(1)

def load_blast_results(file_path):
    """Load BLAST results from QIIME2 export"""
    print(f"Loading BLAST results from {file_path}...")
    
    try:
        # BLAST results columns are typically:
        # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        
        # Read BLAST data
        blast_df = pd.read_csv(file_path, sep='\t', names=names)
        
        # Keep only the best hit for each query sequence
        best_hits = blast_df.sort_values('bitscore', ascending=False).drop_duplicates('qseqid')
        
        print(f"Loaded {len(best_hits)} best BLAST hits out of {len(blast_df)} total hits")
        return best_hits
    
    except Exception as e:
        print(f"Error loading BLAST results: {e}")
        sys.exit(1)

def merge_data(taxonomy_df, blast_df):
    """Merge taxonomy and BLAST data"""
    print("Merging taxonomy and BLAST data...")
    
    # Rename columns for consistency
    blast_df = blast_df.rename(columns={
        'qseqid': 'Feature ID',
        'sseqid': 'SILVA_Accession',
        'pident': 'Percent_Identity',
        'evalue': 'E_value',
        'bitscore': 'Bit_Score'
    })
    
    # Select relevant columns
    blast_df = blast_df[['Feature ID', 'SILVA_Accession', 'Percent_Identity', 'E_value', 'Bit_Score']]
    
    # Merge the dataframes
    merged_df = taxonomy_df.merge(blast_df, on='Feature ID', how='left')
    
    # Clean up SILVA accession (remove version numbers if present)
    if 'SILVA_Accession' in merged_df.columns:
        # Extract the accession number pattern (e.g., "AB123456" from "AB123456.1")
        merged_df['SILVA_Accession'] = merged_df['SILVA_Accession'].str.split('.').str[0]
    
    # Count unmatched features
    unmatched = merged_df['SILVA_Accession'].isna().sum()
    print(f"Merged data: {len(merged_df)} features, {unmatched} without SILVA accessions")
    
    return merged_df

def format_output(merged_df):
    """Format the merged data for output"""
    # Reorder columns for better readability
    columns = ['Feature ID', 'Taxon', 'Confidence', 'SILVA_Accession', 
               'Percent_Identity', 'E_value', 'Bit_Score']
    
    # Ensure all expected columns exist
    for col in columns:
        if col not in merged_df.columns:
            merged_df[col] = None
    
    return merged_df[columns]

def main():
    print("\n=== QIIME2 Taxonomy and SILVA Accessions Merger ===\n")
    
    # Check if input files exist
    for file_path in [TAXONOMY_FILE, BLAST_RESULTS_FILE]:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}")
            print("Please run classify_with_silva_accessions.sh first!")
            sys.exit(1)
    
    # Load data
    taxonomy_df = load_taxonomy(TAXONOMY_FILE)
    blast_df = load_blast_results(BLAST_RESULTS_FILE)
    
    # Merge data
    merged_df = merge_data(taxonomy_df, blast_df)
    
    # Format output
    output_df = format_output(merged_df)
    
    # Save the merged data
    print(f"Saving merged data to {OUTPUT_FILE}...")
    output_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    
    print("\n✅ Merge completed successfully!")
    print(f"Output file: {OUTPUT_FILE}")
    print("\nThe merged file contains:")
    print("  - Feature ID (OTU)")
    print("  - Full Taxonomy")
    print("  - Confidence Score")
    print("  - SILVA Accession Number")
    print("  - Percent Identity")
    print("  - E-value")
    print("  - Bit-score")

if __name__ == "__main__":
    main()
