#!/usr/bin/env python3
"""
Merge rows with the same OTU in both sample_countsV3_V4_444 and taxonomyV3_V4_444 sheets 
of phyloseqV3_V4_444.xlsx. For sample counts, sum the values. For taxonomy, 
take the most complete/specific taxonomy information available.
"""

import pandas as pd
import sys
import os

def merge_taxonomy_for_otu(otu_rows):
    """
    Merge taxonomy information for duplicate OTU rows.
    Takes the most complete information available for each taxonomic level.
    """
    merged_row = otu_rows.iloc[0].copy()  # Start with first row
    
    # Taxonomy columns in order of specificity
    tax_columns = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    for col in tax_columns:
        if col in otu_rows.columns:
            # Get all non-null, non-empty values for this column
            values = otu_rows[col].dropna()
            values = values[values != '']
            values = values[~values.str.startswith(f'{col.lower()[0]}__$')]  # Remove empty taxonomy markers like 'g__'
            
            if len(values) > 0:
                # Take the most frequent value, or first if tie
                value_counts = values.value_counts()
                if len(value_counts) > 0:
                    merged_row[col] = value_counts.index[0]
    
    return merged_row

def main():
    # File paths
    input_file = '/pathway/phyloseqV3_V4_444.xlsx'
    output_file = '/pathway/phyloseqV3_V4_444_merged_complete.xlsx'
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)
    
    print("Reading phyloseqV3_V4_444.xlsx file...")
    
    # Read all sheets
    xl = pd.ExcelFile(input_file)
    all_sheets = {}
    
    for sheet_name in xl.sheet_names:
        print(f"Reading sheet: {sheet_name}")
        all_sheets[sheet_name] = pd.read_excel(input_file, sheet_name=sheet_name)
    
    # Process sample_countsV3_V4_444 sheet
    if 'sample_countsV3_V4_444' not in all_sheets:
        print("Error: sample_countsV3_V4_444 sheet not found")
        sys.exit(1)
    
    df_counts = all_sheets['sample_countsV3_V4_444']
    
    print(f"\\nProcessing sample_countsV3_V4_444 sheet:")
    print(f"Original: {df_counts.shape[0]} rows, {df_counts.shape[1]} columns")
    print(f"Unique OTUs before merging: {df_counts['otu'].nunique()}")
    
    # Check for duplicates
    duplicates_counts = df_counts[df_counts.duplicated(subset=['otu'], keep=False)]
    print(f"Found {len(duplicates_counts)} duplicate rows affecting {duplicates_counts['otu'].nunique()} unique OTUs")
    
    # Get the sample columns (all columns except 'otu')
    sample_columns = [col for col in df_counts.columns if col != 'otu']
    
    # Group by OTU and sum the counts for each sample
    print("Merging duplicate OTUs by summing counts...")
    df_counts_merged = df_counts.groupby('otu', as_index=False)[sample_columns].sum()
    
    print(f"Merged: {df_counts_merged.shape[0]} rows, {df_counts_merged.shape[1]} columns")
    print(f"Removed {df_counts.shape[0] - df_counts_merged.shape[0]} duplicate rows")
    
    # Process taxonomyV3_V4_444 sheet
    if 'taxonomyV3_V4_444' in all_sheets:
        df_tax = all_sheets['taxonomyV3_V4_444']
        
        print(f"\\nProcessing taxonomyV3_V4_444 sheet:")
        print(f"Original: {df_tax.shape[0]} rows, {df_tax.shape[1]} columns")
        print(f"Unique OTUs before merging: {df_tax['otu'].nunique()}")
        
        # Check for duplicates
        duplicates_tax = df_tax[df_tax.duplicated(subset=['otu'], keep=False)]
        print(f"Found {len(duplicates_tax)} duplicate rows affecting {duplicates_tax['otu'].nunique()} unique OTUs")
        
        if len(duplicates_tax) > 0:
            print("Merging duplicate taxonomy entries...")
            
            # Group by OTU and merge taxonomy
            merged_tax_rows = []
            
            for otu in df_tax['otu'].unique():
                otu_rows = df_tax[df_tax['otu'] == otu]
                if len(otu_rows) > 1:
                    # Multiple rows for this OTU - merge them
                    merged_row = merge_taxonomy_for_otu(otu_rows)
                    merged_tax_rows.append(merged_row)
                else:
                    # Single row - keep as is
                    merged_tax_rows.append(otu_rows.iloc[0])
            
            df_tax_merged = pd.DataFrame(merged_tax_rows)
            
            print(f"Merged: {df_tax_merged.shape[0]} rows, {df_tax_merged.shape[1]} columns")
            print(f"Removed {df_tax.shape[0] - df_tax_merged.shape[0]} duplicate rows")
            
            all_sheets['taxonomyV3_V4_444'] = df_tax_merged
        else:
            print("No duplicates found in taxonomy sheet")
    
    # Update the sheets dictionary with merged data
    all_sheets['sample_countsV3_V4_444'] = df_counts_merged
    
    # Verify that sample_countsV3_V4_444 and taxonomyV3_V4_444 have matching OTU sets
    if 'taxonomyV3_V4_444' in all_sheets:
        counts_otus = set(all_sheets['sample_countsV3_V4_444']['otu'])
        tax_otus = set(all_sheets['taxonomyV3_V4_444']['otu'])
        
        print(f"\\nOTU consistency check:")
        print(f"OTUs in sample_countsV3_V4_444: {len(counts_otus)}")
        print(f"OTUs in taxonomyV3_V4_444: {len(tax_otus)}")
        print(f"Common OTUs: {len(counts_otus.intersection(tax_otus))}")
        
        if counts_otus != tax_otus:
            print("Warning: OTU sets don't match between sheets")
            missing_in_tax = counts_otus - tax_otus
            missing_in_counts = tax_otus - counts_otus
            if missing_in_tax:
                print(f"  Missing in taxonomy: {len(missing_in_tax)} OTUs")
            if missing_in_counts:
                print(f"  Missing in counts: {len(missing_in_counts)} OTUs")
        else:
            print("✓ OTU sets match perfectly between sheets")
    
    # Save all sheets to the new Excel file
    print(f"\\nSaving merged data to: {output_file}")
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for sheet_name, sheet_df in all_sheets.items():
            print(f"Writing sheet: {sheet_name}")
            sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    print("\\nMerging completed successfully!")
    print(f"Original file: {input_file}")
    print(f"Merged file: {output_file}")

if __name__ == "__main__":
    main()
