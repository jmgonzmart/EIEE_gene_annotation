# EIEE gene annotation

This repository contains several scripts used for the EIEE-related gene annotation comparison in the publication "Re-annotation of 191 developmental and epileptic encephalopathy-associated genes unmasks de novo variants in SCN1A" by Steward et al., NPJ Genom Med. 2019; 4: 31.

a) Novel exon sequence after gene re-annotation
The compare_annotation.pl script counts novel exons, novel introns, shifted splice junctions and novel genomic coverage after the re-annotation of the EIEE-related genes. The comparison is made between the GENCODE v20 and GENCODE v28 releases.

perl compare_annotation.pl -in EIEE_gene_list.txt -sp human -out annotation_comparison_output.txt

b) Variants in the novel sequences after gene re-annotation
The compare_variation.pl script counts the variants overlapping the annotation of the EIEE-related genes in the GENCODE v20 and GENCODE v28 releases and which of them are uniqu to the novel annotation. The variants are extracted from the Ensembl variation dataset "All phenotype-associated - short variants (SNPs and indels)".

perl compare_variation.pl -in EIEE_gene_list.txt -sp human -out variation_comparison_output.txt


## Dependencies
- Ensembl Perl API (recommended v92)
- Getopt::Long
