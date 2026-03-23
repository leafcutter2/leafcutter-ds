leafcutter-gtf-to-exons annotation/chr10_minimal.gtf.gz annotation/chr10_minimal.exons.txt
leafcutter-ds -0 foobar --exon_file annotation/chr10_minimal.exons.txt -o ./test leafcutter2.junction_counts.gz groups.txt
leafcutter-ds -0 Cells-EBV-transformedlymphocytes --exon_file annotation/chr10_minimal.exons.txt -o ./test leafcutter2.junction_counts.gz groups.txt

