#!/bin/bash
# Take the list of DEG
grep -f ../results/motif_finding/ST6512vsOKYL029_DEP.txt \
        ../data/ref_genome/genomic.gtf \
        > ../results/motif_finding/filtered_gtf.txt

# Only retain Gene coordinates
awk '$3 == "gene"' ../results/motif_finding/filtered_gtf.txt \
    > ../results/motif_finding/filtered_gtd_gene_only.txt

# If genes are on the positive strand, take start-1000
# If genes are on the negative strand, do end+1000
awk -F'\t' \
    '{if ($7 == "+") $4 = $4 - 1000; else if ($7 == "-") $5 = $5 + 1000; print}' \
    OFS='\t' \
    ../results/motif_finding/filtered_gtd_gene_only.txt > ../results/motif_finding/gene_and_promoter.txt

# Only extract promoter location
# If gene is on the + strand take first 1000 bp
# If gene is on the - only take last 1000 bp
awk -F'\t' \
    '{if ($7 == "+") $5 = $4 + 1000; else if ($7 == "-") $4 = $5 - 1000; print}' \
    OFS='\t' \
    ../results/motif_finding/gene_and_promoter.txt > ../results/motif_finding/promoter_location.txt


#bedtools to extract sequences
bedtools getfasta -fi ../data/ref_genome/GCA_000002525.1_ASM252v1_genomic.fna \
        -bed ../results/motif_finding/promoter_location.txt \
        > ../results/motif_finding/DNA_seq.txt

#checked the orientation of promoters vs genes on benchling and + and - are correctly assined/calculated