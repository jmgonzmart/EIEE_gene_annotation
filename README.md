# EIEE_gene_annotation

This repository contains several scripts used for the EIEE-related gene annotation comparison in the manuscript "Systematic re-annotation of 191 genes associated with early-onset epilepsy unmasks de novo variants linked to Dravet syndrome in novel SCN1A exons" by Steward et al., submitted to Nature Genomic Medicine.

a) Novel exon sequence after gene re-annotation
The compare_annotation.pl script counts novel exons, novel introns, shifted splice junctions and novel genomic coverage after the re-annotation of the EIEE-related genes. The comparison is made between the GENCODE v20 and GENCODE v28 releases.

perl compare_annotation.pl -in EIEE_gene_list.txt -sp human -out annotation_comparison_output.txt

b) Variants in the novel sequences after gene re-annotation
The compare_variation.pl script counts the variants overlapping the annotation of the EIEE-related genes in the GENCODE v20 and GENCODE v28 releases and which of them are uniqu to the novel annotation. The variants are extracted from the Ensembl variation dataset "All phenotype-associated - short variants (SNPs and indels)".

perl compare_variation.pl -in EIEE_gene_list.txt -sp human -out variation_comparison_output.txt


Dependencies
- Ensembl API (recommended v92)
- Getopt::Long
