#!/bin/bash

# QIIME2 Classification with SILVA Accession Numbers
# This script modifies your original classification command to output 
# both taxonomy classifications AND the specific SILVA accession numbers
# that each sequence was matched against.

echo "Starting QIIME2 classification with SILVA accession tracking..."


# Input files
REP_SEQS="rep-seqs-deblurV3_V4_444.qza"
SILVA_SEQS="silva-138.2-ssu-nr99-seqs-derep-uniq.qza"
SILVA_TAX="silva-138.2-ssu-nr99-tax-derep-uniq.qza"

# Output files
CLASSIFICATION_OUT="bespoke-taxonomyV3_V4_444-with-accessions.qza"
BLAST_RESULTS_OUT="silva-blast-results-V3_V4_444.qza"

# Check if input files exist
echo "Checking input files..."
for file in "$REP_SEQS" "$SILVA_SEQS" "$SILVA_TAX"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file $file not found!"
        echo "Please ensure all files are in the current directory."
        exit 1
    fi
    echo "✓ Found: $file"
done
echo ""

# Run QIIME2 BLAST-based classification
echo "Running QIIME2 BLAST-based classification..."
echo "This will take longer than sklearn but provides SILVA accession numbers."
echo ""

qiime feature-classifier classify-consensus-blast \
    --i-query "$REP_SEQS" \
    --i-reference-reads "$SILVA_SEQS" \
    --i-reference-taxonomy "$SILVA_TAX" \
    --p-maxaccepts 10 \
    --p-perc-identity 0.8 \
    --p-query-cov 0.8 \
    --p-min-consensus 0.51 \
    --p-num-threads 4 \
    --o-classification "$CLASSIFICATION_OUT" \
    --o-search-results "$BLAST_RESULTS_OUT" \
    --verbose

# Check if command succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Classification completed successfully!"
    echo ""
    echo "Output files created:"
    echo "  1. $CLASSIFICATION_OUT - Taxonomy classifications"
    echo "  2. $BLAST_RESULTS_OUT - BLAST results with SILVA accessions"
    echo ""
    
    # Export the results for analysis
    echo "Exporting results for analysis..."
    
    # Export taxonomy
    qiime tools export \
        --input-path "$CLASSIFICATION_OUT" \
        --output-path "taxonomyV3_V4_444_with_accessions_export"
    
    # Export BLAST results  
    qiime tools export \
        --input-path "$BLAST_RESULTS_OUT" \
        --output-path "silva_blast_results_export"
    
    echo ""
    echo "✅ Export completed!"
    echo ""
    echo "Exported files:"
    echo "  - taxonomyV3_V4_444_with_accessions_export/taxonomy.tsv"
    echo "  - silva_blast_results_export/search-results.tsv"
    echo ""
    echo "The search-results.tsv file contains:"
    echo "  - Query sequence ID (your OTU)"
    echo "  - Subject sequence ID (SILVA accession)" 
    echo "  - Percent identity"
    echo "  - Alignment length"
    echo "  - E-value"
    echo "  - Bit score"
    echo ""
    
else
    echo "❌ Classification failed!"
    echo "Please check the error messages above."
    exit 1
fi

echo "================================================================"
echo "ANALYSIS COMPLETE!"
echo "================================================================"
