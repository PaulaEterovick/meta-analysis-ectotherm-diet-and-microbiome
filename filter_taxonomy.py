#!/usr/bin/env python3
"""
Filter silva_taxonomy_mappingV3_V4_444.tsv to keep only Feature IDs 
that are present as #OTU ID in sample_countsV3_V4_444.xlsx
"""

import pandas as pd
import sys
import os

def main():
    # File paths
    excel_file = '/Volumes/TOSHIBA EXT/metanalysis/clphyloseq_files/sample_countsV3_V4_444.xlsx'
    taxonomy_file = '/Volumes/TOSHIBA EXT/metanalysis/clphyloseq_files/silva_taxonomy_mappingV3_V4_444.tsv'
    output_file = '/Volumes/TOSHIBA EXT/metanalysis/clphyloseq_files/silva_taxonomy_mappingV3_V4_444_filtered.tsv'
    
    # Check if files exist
    if not os.path.exists(excel_file):
        print(f"Error: Excel file not found: {excel_file}")
        sys.exit(1)
    
    if not os.path.exists(taxonomy_file):
        print(f"Error: Taxonomy file not found: {taxonomy_file}")
        sys.exit(1)
    
    print("Reading OTU IDs from Excel file...")
    # Read the Excel file to get OTU IDs (header is at row 2, 0-indexed)
    df_excel = pd.read_excel(excel_file, header=2)
    
    # Extract OTU IDs from the first column (#OTU ID)
    otu_ids = set(df_excel.iloc[:, 0].tolist())
    
    print(f"Found {len(otu_ids)} unique OTU IDs in Excel file")
    
    print("Reading taxonomy mapping file...")
    # Read the taxonomy file
    df_taxonomy = pd.read_csv(taxonomy_file, sep='\t')
    
    print(f"Original taxonomy file has {len(df_taxonomy)} rows")
    print(f"Columns: {list(df_taxonomy.columns)}")
    
    # Filter taxonomy file to keep only rows where Feature ID is in the OTU IDs set
    df_filtered = df_taxonomy[df_taxonomy['Feature ID'].isin(otu_ids)]
    
    print(f"Filtered taxonomy file has {len(df_filtered)} rows")
    print(f"Removed {len(df_taxonomy) - len(df_filtered)} rows")
    
    # Check how many OTU IDs from Excel were matched
    matched_ids = set(df_filtered['Feature ID'].tolist())
    unmatched_ids = otu_ids - matched_ids
    
    print(f"Matched {len(matched_ids)} out of {len(otu_ids)} OTU IDs from Excel file")
    if unmatched_ids:
        print(f"Warning: {len(unmatched_ids)} OTU IDs from Excel file were not found in taxonomy file")
        if len(unmatched_ids) <= 10:
            print("Unmatched IDs:", list(unmatched_ids))
    
    # Save the filtered file
    print(f"Saving filtered taxonomy file to: {output_file}")
    df_filtered.to_csv(output_file, sep='\t', index=False)
    
    print("Filtering completed successfully!")
    print(f"Original file: {len(df_taxonomy)} rows")
    print(f"Filtered file: {len(df_filtered)} rows")

if __name__ == "__main__":
    main()
